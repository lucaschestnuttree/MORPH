from __future__ import annotations

import os
import re
import sys
import functools
from copy import deepcopy
from datetime import datetime
from io import StringIO
from collections import defaultdict
from dataclasses import dataclass
from string import ascii_uppercase
from typing import Any, Dict, List, Tuple, Optional, Set

import numpy as np
import pandas as pd

from Bio.PDB import PDBParser, PDBIO, Select, PPBuilder
from Bio.PDB import NeighborSearch
from Bio.PDB.Selection import unfold_entities

from pymut import *  # mutate(), load_rotamers(), is_backbone(), ...


# Ensure stdout/stderr are line-buffered for live logging in SLURM
try:
    sys.stdout.reconfigure(line_buffering=True)
    sys.stderr.reconfigure(line_buffering=True)
except Exception:
    pass

fprint = functools.partial(print, flush=True)

# ============================================================
# KNOBS
# ============================================================

PYMUT_MUTATION_TYPE: str = "first"
PYMUT_FALLBACK_MUTATION_TYPE: str = "best"
MAX_FALLBACK_REDO: int = 10000

# Distance-cutoff clash detector (Å), SAME CHAIN ONLY
CLASH_NEIGHBOR_CUTOFF_A: float = 8.0
CLASH_HARD_CUTOFF_A: float = 3.0

DEFAULT_BFACTOR = 42.0
BACKBONE_NAMES = {"N", "CA", "C", "O", "OXT"}

# For ancestry categories we assign chains starting at U, V, W, X, Y, Z...
CHAIN_LABELS = list(ascii_uppercase[19:])  # U..Z


# ============================================================
# Residue name maps
# ============================================================

AA1_TO_AA3: Dict[str, str] = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE",
    "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}

AA3_NORMALIZE: Dict[str, str] = {
    "HID": "HIS", "HIE": "HIS", "HIP": "HIS",
    "CYX": "CYS", "CYN": "CYS", "CSS": "CYS", "CYM": "CYS",
    "MSE": "MET",
}


# ============================================================
# Data structures
# ============================================================

@dataclass(frozen=True)
class Mutation:
    wt3: str
    pos: int
    mut3: str


class NoHydrogensSelect(Select):
    def accept_atom(self, atom):
        return atom.element != "H"


class AtomRecordsOnlySelect(Select):
    """Write only non-hydrogen atoms from standard (ATOM) residues."""

    def accept_residue(self, residue):
        hetflag = residue.id[0]
        return not hetflag.strip()

    def accept_atom(self, atom):
        return atom.element != "H"


# ============================================================
# TSV processing
# ============================================================

def filter_tsv_file(file_path: str) -> pd.DataFrame:
    """
    Filters rows in a TSV file based on specific criteria, including allele_count,
    protein_variant, allele_frequency, intron, uniprot_id, and am_pathogenicity.
    """
    df = pd.read_csv(file_path, sep="\t")

    if "intron" not in df.columns:
        df["intron"] = np.nan

    df_filtered = df[
        (df["allele_count"].notna())
        & (df["allele_count"] != 0.0)
        & (df["intron"].isna() | df["intron"].isin([None, ""]))
        & (df["uniprot_id"] != "splice")
        & (~df["am_class"].isin(["new_stop", "synonymous"]))
        & (~df["protein_variant"].str.contains(r"[A-Z]1[A-Z]", na=False))
    ]

    return df_filtered


def extract_position(variant: str) -> int:
    match = re.search(r"\d+", variant)
    if match:
        return int(match.group())
    raise ValueError(f"Invalid variant format: {variant}")


def map_frequency_to_bin(frequency: float) -> int:
    frequency_bins = [
        (1e-9, 1e-8, 0),
        (1e-8, 1e-7, 1),
        (1e-7, 1e-6, 2),
        (1e-6, 1e-5, 3),
        (1e-5, 1e-4, 4),
        (1e-4, 1e-3, 5),
        (1e-3, 1e-2, 6),
        (1e-2, 0.51, 7),
        (0.51, 1.0, 8),
    ]
    for lower, upper, value in frequency_bins:
        if lower <= frequency < upper:
            return value
    return 0


def calculate_top_frequencies_per_category(
    df: pd.DataFrame, categories: List[str]
) -> Dict[str, Dict[int, Dict[str, Any]]]:
    """
    For each category, per position keep the MAX frequency and the last_character (mut AA1).
    """
    category_position_data = {category: {} for category in categories}

    for category in categories:
        column_name = (
            "allele_frequency" if category == "allele_frequency" else f"{category}_allele_frequency"
        )
        if column_name not in df.columns:
            print(f"Warning: Column '{column_name}' not found. Skipping category '{category}'.")
            continue

        category_data = defaultdict(lambda: {"frequency": 0.0, "last_character": None})

        for _, row in df.iterrows():
            if pd.isna(row[column_name]) or float(row[column_name]) <= 0.0:
                continue
            try:
                position = extract_position(row["protein_variant"])
                frequency = float(row[column_name])
                last_character = str(row["protein_variant"])[-1]

                if frequency > category_data[position]["frequency"]:
                    category_data[position]["frequency"] = frequency
                    category_data[position]["last_character"] = last_character
            except ValueError as e:
                print(f"Skipping invalid row: {e}")

        category_position_data[category] = dict(category_data)

    return category_position_data


def chainA_frequency(top_category_data: Dict[str, Dict[int, Dict[str, Any]]]) -> Dict[int, float]:
    """
    Max frequency across all categories per position (for Chain A B-factors).
    """
    position_max: Dict[int, float] = {}
    for _, positions in top_category_data.items():
        for pos, data in positions.items():
            f = float(data.get("frequency", 0.0) or 0.0)
            position_max[pos] = max(position_max.get(pos, 0.0), f)
    return position_max


# ============================================================
# PDB helpers
# ============================================================

def parse_structure(pdb_path: str, name: str = "pdb"):
    parser = PDBParser(QUIET=True)
    return parser.get_structure(name, pdb_path)


def get_chain_resname_map(pdb_path: str, chain_id: str) -> Dict[int, str]:
    structure = parse_structure(pdb_path, "resmap")
    model = structure[0]
    if chain_id not in [c.id for c in model]:
        raise ValueError(f"Chain {chain_id} not found in {pdb_path}")
    chain = model[chain_id]

    res_map: Dict[int, str] = {}
    for res in chain:
        hetflag, resseq, icode = res.id
        if hetflag.strip():
            continue
        if icode != " ":
            continue
        name = res.resname.strip().upper()
        name = AA3_NORMALIZE.get(name, name)
        res_map[resseq] = name
    return res_map


def free_reserved_chain_ids(model, reserved_chain_ids: Set[str]) -> Dict[str, str]:
    existing = {c.id for c in model}
    mapping: Dict[str, str] = {}
    candidates = (
        list("BCDEFGHIJKLMNOPQRSTUVWXYZ")
        + list("0123456789")
        + list("abcdefghijklmnopqrstuvwxyz")
    )

    def pick_unused() -> str:
        for cid in candidates:
            if cid not in existing and cid not in reserved_chain_ids:
                return cid
        raise RuntimeError("Could not find an unused chain id.")

    for old in sorted(reserved_chain_ids):
        if old in existing:
            new = pick_unused()
            model[old].id = new
            mapping[old] = new
            existing.remove(old)
            existing.add(new)
    return mapping


def normalize_pdb_chains(
    pdb_in: str,
    pdb_out: str,
    reserved: Set[str],
    normalize_residues: bool = True,
) -> Dict[str, str]:
    structure = parse_structure(pdb_in, "norm")
    model = structure[0]

    mapping = free_reserved_chain_ids(model, reserved_chain_ids=reserved)

    if normalize_residues:
        for chain in model:
            for res in chain:
                hetflag, _, icode = res.id
                if hetflag.strip():
                    continue
                if icode != " ":
                    continue
                old = res.resname.strip().upper()
                new = AA3_NORMALIZE.get(old)
                if new and new != old:
                    res.resname = new

    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_out, preserve_atom_numbering=True)
    return mapping


def extract_chainA_as_chain(
    pdb_in: str,
    pdb_out: str,
    chain_dst: str,
    chain_src: str = "A",
) -> None:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("src", pdb_in)
    model = structure[0]
    if chain_src not in [c.id for c in model]:
        raise ValueError(f"Chain {chain_src} not found in {pdb_in}")

    src_chain = model[chain_src]
    new_chain = deepcopy(src_chain)
    new_chain.id = chain_dst

    one_structure = deepcopy(structure)
    one_model = one_structure[0]
    for ch in list(one_model):
        one_model.detach_child(ch.id)
    one_model.add(new_chain)

    io = PDBIO()
    io.set_structure(one_structure)
    # The category-chain template must contain only records that were ATOM
    # records in chain A. HETATM residues from chain A are intentionally omitted.
    io.save(
        pdb_out,
        select=AtomRecordsOnlySelect(),
        preserve_atom_numbering=True,
    )


# ============================================================
# FAST clash helpers
# ============================================================

def build_same_chain_neighborsearch(structure, chain_id: str):
    model = structure[0]
    chain = model[chain_id]
    atoms_same_chain = list(unfold_entities(chain, "A"))
    return chain, NeighborSearch(atoms_same_chain)


def clash_and_min_distance_same_chain(
    ns: NeighborSearch,
    residue,
    *,
    neighbor_cutoff: float = CLASH_NEIGHBOR_CUTOFF_A,
    hard_cutoff: float = CLASH_HARD_CUTOFF_A,
) -> Tuple[bool, float]:
    min_dist = float("inf")
    has_clash = False

    for atom in residue:
        if is_backbone(atom) or atom.get_id().strip() in BACKBONE_NAMES:
            continue
        if getattr(atom, "element", None) == "H":
            continue

        for close_atom in ns.search(atom.coord, neighbor_cutoff):
            close_res = close_atom.get_parent()
            if close_res is residue:
                continue
            if getattr(close_atom, "element", None) == "H":
                continue

            dist = float(np.linalg.norm(close_atom.coord - atom.coord))
            if dist < min_dist:
                min_dist = dist
            if dist < hard_cutoff:
                has_clash = True

    return has_clash, min_dist


# ============================================================
# Category mutation planning (ancestry-specific)
# ============================================================

def plan_mutations_for_category(
    category_data: Dict[int, Dict[str, Any]],
    chainA_res_map: Dict[int, str],
) -> List[Mutation]:
    """
    For each position with a max-frequency variant in this category,
    build Mutation(wt3,pos,mut3) if WT != MUT.
    """
    muts: List[Mutation] = []
    for pos, data in category_data.items():
        aa1 = data.get("last_character")
        if not aa1 or not isinstance(aa1, str):
            continue
        aa1 = aa1.upper()
        mut3 = AA1_TO_AA3.get(aa1)
        if mut3 is None:
            continue

        wt3 = chainA_res_map.get(int(pos))
        if wt3 is None:
            continue

        mut3 = AA3_NORMALIZE.get(mut3, mut3)
        wt3 = AA3_NORMALIZE.get(wt3, wt3)

        if wt3 == mut3:
            continue

        muts.append(Mutation(wt3=wt3, pos=int(pos), mut3=mut3))
    return muts


# ============================================================
# pyMUT mutation + re-rotamerization
# ============================================================

def apply_mutations_with_pymut(
    pdb_in: str,
    pdb_out: str,
    chain_id: str,
    mutations: List[Mutation],
) -> None:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("mut", pdb_in)

    try:
        rotamer_lib = load_rotamers()
    except Exception as e:
        rotamer_lib = None
        fprint(f"Warning: could not preload rotamer library: {e}")

    mutated_positions: List[int] = []
    desired_res_by_pos: Dict[int, str] = {}

    # Phase 1: apply mutations
    for idx, m in enumerate(mutations, start=1):
        pos = int(m.pos)
        new_res = AA3_NORMALIZE.get(m.mut3, m.mut3)
        desired_res_by_pos[pos] = new_res
        try:
            mutate(
                structure,
                chain_id,
                pos,
                new_res,
                mutation_type=PYMUT_MUTATION_TYPE,
                rotamer_lib=rotamer_lib,
            )
            mutated_positions.append(pos)
        except Exception as e:
            fprint(f"Warning: pyMUT failed {chain_id}:{pos} -> {new_res} using {PYMUT_MUTATION_TYPE}: {e}")

        if idx % 50 == 0 or idx == len(mutations):
            fprint(f" - pyMUT phase1 {chain_id}: {idx}/{len(mutations)} mutations applied")

    if not mutated_positions:
        io = PDBIO()
        io.set_structure(structure)
        io.save(pdb_out)
        return

    # Phase 2: clash detect + fallback re-rotamerize worst
    try:
        chain, ns = build_same_chain_neighborsearch(structure, chain_id)
    except Exception as e:
        fprint(f"Warning: could not build NeighborSearch for chain {chain_id}: {e}")
        io = PDBIO()
        io.set_structure(structure)
        io.save(pdb_out)
        return

    clashing: List[Tuple[float, int]] = []
    for pos in mutated_positions:
        try:
            residue = chain[pos]
        except Exception:
            continue

        has_clash, min_dist = clash_and_min_distance_same_chain(
            ns,
            residue,
            neighbor_cutoff=CLASH_NEIGHBOR_CUTOFF_A,
            hard_cutoff=CLASH_HARD_CUTOFF_A,
        )
        if has_clash:
            if min_dist == float("inf"):
                min_dist = 0.0
            clashing.append((float(min_dist), int(pos)))

    if not clashing:
        fprint(f" - pyMUT phase2 {chain_id}: no clashes found (hard<{CLASH_HARD_CUTOFF_A:.2f}Å)")
        io = PDBIO()
        io.set_structure(structure)
        io.save(pdb_out)
        return

    clashing.sort(key=lambda t: t[0])
    worst = clashing[:MAX_FALLBACK_REDO]

    n_redo = 0
    for pre_min_dist, pos in worst:
        new_res = desired_res_by_pos[pos]
        try:
            mutate(
                structure,
                chain_id,
                pos,
                new_res,
                mutation_type=PYMUT_FALLBACK_MUTATION_TYPE,
                rotamer_lib=rotamer_lib,
            )
            # rebuild NS after coordinate changes
            chain, ns = build_same_chain_neighborsearch(structure, chain_id)
            n_redo += 1
        except Exception as e:
            fprint(f"Warning: pyMUT fallback failed {chain_id}:{pos} -> {new_res} using {PYMUT_FALLBACK_MUTATION_TYPE}: {e}")

    fprint(
        f" - pyMUT phase2 {chain_id}: clashing sites={len(clashing)}, "
        f"remutated={n_redo}/{min(len(clashing), MAX_FALLBACK_REDO)} "
        f"(hard<{CLASH_HARD_CUTOFF_A:.2f}Å; neighbor<{CLASH_NEIGHBOR_CUTOFF_A:.1f}Å)"
    )

    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_out)



# ============================================================
# B-factor annotation + merge (ancestry version)
# ============================================================

def set_chain_bfactors_from_category(chain, category_data: Dict[int, Dict[str, Any]]) -> None:
    """
    For a category chain (U/V/W/X/Y/...), B-factor comes from that category's frequency at residue.
    """
    for residue in chain:
        resnum = residue.id[1]
        if resnum in category_data:
            freq = float(category_data[resnum].get("frequency", 0.0) or 0.0)
            b = float(map_frequency_to_bin(freq))
        else:
            b = float(DEFAULT_BFACTOR)
        for atom in residue:
            atom.bfactor = b


def set_chainA_bfactors_from_max_category_freq(
    structure,
    top_category_data: Dict[str, Dict[int, Dict[str, Any]]],
    default_b: float = DEFAULT_BFACTOR,
) -> None:
    """
    Chain A B-factor = max bin across categories per residue position.
    """
    model = structure[0]
    if "A" not in [c.id for c in model]:
        raise ValueError("Chain A not found in merged structure.")

    # pos -> max allele frequency across all categories
    pos_maxfreq = chainA_frequency(top_category_data)

    for residue in model["A"]:
        resnum = residue.id[1]
        if resnum in pos_maxfreq:
            b = float(map_frequency_to_bin(float(pos_maxfreq[resnum] or 0.0)))
        else:
            b = float(default_b)

        for atom in residue:
            atom.bfactor = b
            
def set_chainA_bfactors_from_max_category_freq(
    structure,
    top_category_data: Dict[str, Dict[int, Dict[str, Any]]],
    default_b: float = DEFAULT_BFACTOR,
) -> None:
    """
    Chain A B-factor = max bin across categories per residue position.
    """
    model = structure[0]
    if "A" not in [c.id for c in model]:
        raise ValueError("Chain A not found in merged structure.")

    # pos -> max allele frequency across all categories
    pos_maxfreq = chainA_frequency(top_category_data)

    for residue in model["A"]:
        resnum = residue.id[1]
        if resnum in pos_maxfreq:
            b = float(map_frequency_to_bin(float(pos_maxfreq[resnum] or 0.0)))
        else:
            b = float(default_b)

        for atom in residue:
            atom.bfactor = b


def merge_normalized_original_with_variant_chains(
    normalized_original_pdb: str,
    variant_chain_pdbs: List[Tuple[str, str, str]],  # (pdb_path, chain_label, category_name)
    top_category_data: Dict[str, Dict[int, Dict[str, Any]]],
) -> object:
    """
    Final merge logic:
      1) start from normalized original (keeps ALL original chains)
      2) overwrite Chain A B-factors = max bin across categories per position
      3) add/replace variant chains U/V/W/X/Y/... with category-specific B-factors
    """
    original_structure = parse_structure(normalized_original_pdb, "original_norm")
    model = original_structure[0]

    # (2) Set chain A b-factors from max across categories
    set_chainA_bfactors_from_max_category_freq(original_structure, top_category_data, DEFAULT_BFACTOR)

    # (3) Add/replace category chains
    for pdb_path, chain_label, category in variant_chain_pdbs:
        s = parse_structure(pdb_path, chain_label)
        m = s[0]
        if chain_label not in [c.id for c in m]:
            raise ValueError(f"Expected chain {chain_label} not found in {pdb_path}")

        new_chain = deepcopy(m[chain_label])

        # category-chain b-factors from that category's frequency
        set_chain_bfactors_from_category(new_chain, top_category_data.get(category, {}))

        # Give copied/mutated chains unique atom serials above every atom in the
        # normalized original. This lets original CONECT records continue to
        # reference the original ATOM/HETATM serial numbers.
        existing_serials = [
            int(getattr(atom, "serial_number", 0) or 0)
            for atom in original_structure.get_atoms()
        ]
        next_serial = max(existing_serials, default=0) + 1
        for atom in new_chain.get_atoms():
            atom.serial_number = next_serial
            next_serial += 1

        if chain_label in [c.id for c in model]:
            model.detach_child(chain_label)
        model.add(new_chain)

    return original_structure


# ============================================================
# ONE function that builds ALL remarks (category mapping + sequences)
# ============================================================

def build_metadata_lines_with_category_mapping_and_sequences(
    original_pdb: str,
    structure,  # merged Biopython Structure
    *,
    categories: List[str],
    chain_labels: List[str],
    chains_to_report: Optional[Set[str]] = None,
) -> List[str]:
    """
    Builds metadata lines by:
      - taking non-coordinate records from original_pdb
      - inserting ONE REMARK   1 block with:
          * category mapping + disclaimers
          * sequences for requested chains
    Inserted before the first existing REMARK   1 line (or at end if none).
    """
    if chains_to_report is None:
        chains_to_report = set(["A"] + list(chain_labels))

    with open(original_pdb, "r") as f:
        lines = f.readlines()

    drop_prefixes = ("ATOM", "HETATM", "ANISOU", "MODEL", "END", "ENDMDL", "TER", "CONECT", "MASTER")
    metadata = [ln for ln in lines if not ln.startswith(drop_prefixes)]

    today_date = datetime.now().strftime("%d-%b-%y").upper()
    metadata = [
        f"HEADER    GNOMAD CATEGORY MAPPING                          {today_date}\n"
        if ln.startswith("HEADER") else ln
        for ln in metadata
    ]

    # ---------- build mapping/disclaimer lines (ALL REMARK 1) ----------
    remark_1_lines = [
        "REMARK   1\n",
        "REMARK   1 REFERENCE\n",
        "REMARK   1  AUTH   LUCAS CHATAIGNER & HANNA JULIENNE\n",
        "REMARK   1  TITL   GNOMAD CATEGORY (ANCESTRY/SEX) STRUCTURE MAPPING\n",
        "REMARK   1\n",
        "REMARK   1 DISCLAIMERS\n",
        "REMARK   1 PROVIDED HERE IS A HEAVILY ANNOTATED FILE FOR VISUALIZATION PURPOSES\n",
        "REMARK   1\n",
        "REMARK   1 B-FACTORS ARE ARBITRARILY ASSIGNED ACCORDING TO A SCHEME USEFUL.\n",
        "REMARK   1 FOR VISUALIZATION OF GNOMAD ANCESTRY SPECIFIC ALLELE FREQUENCY\n",
        "REMARK   1 CHAINS U/V/W/X/Y/... REPRESENT MOST COMMON VARIANTS PER CATEGORY AT EACH POSITION.\n",
        "REMARK   1 THESE CHAINS ARE NOT REAL PROTEINS AND CONTAIN UNREALISTIC STRUCTURAL ELEMENTS..\n",
        "REMARK   1\n",
    ]
    for cl, cat in zip(chain_labels, categories):
        remark_1_lines.append(f"REMARK   1 CHAIN {cl} CORRESPONDS TO CATEGORY: {cat.upper()}\n")

    # ---------- sequences (ALL REMARK 1) ----------
    ppb = PPBuilder()
    chain_sequences: Dict[str, str] = {}
    model = structure[0]

    for chain in model:
        if chain.id not in chains_to_report:
            continue
        peptides = ppb.build_peptides(chain)
        seq = "".join(str(p.get_sequence()) for p in peptides)
        if seq:
            chain_sequences[chain.id] = seq

    seq_lines = [
        "REMARK   1\n",
        "REMARK   1 POLYPEPTIDE SEQUENCES (ONE-LETTER) FOR CHAINS:\n",
    ]
    for cid in sorted(chain_sequences.keys()):
        seq_lines.append(f"REMARK   1  CHAIN {cid}: {chain_sequences[cid]}\n")

    # Insert before first existing REMARK   1 (or at end if none)
    first_remark_1_index = next(
        (i for i, ln in enumerate(metadata) if ln.startswith("REMARK   1")),
        len(metadata),
    )

    combined = (
        metadata[:first_remark_1_index]
        + remark_1_lines
        + seq_lines
        + metadata[first_remark_1_index:]
    )
    return combined

def _remap_chain_id_in_coordinate_line(line: str, chain_mapping: Dict[str, str]) -> str:
    """Apply normalized chain-id changes to an original coordinate record."""
    if len(line) > 21:
        old_chain = line[21]
        if old_chain in chain_mapping:
            line = line[:21] + chain_mapping[old_chain] + line[22:]
    return line


def read_original_hetatm_and_conect_records(
    original_pdb: str,
    chain_mapping: Optional[Dict[str, str]] = None,
) -> Tuple[List[str], List[str]]:
    """Return original HETATM and CONECT records, preserving their text/order."""
    chain_mapping = chain_mapping or {}
    hetatm_lines: List[str] = []
    conect_lines: List[str] = []

    with open(original_pdb, "r") as handle:
        for line in handle:
            if line.startswith("HETATM"):
                line = _remap_chain_id_in_coordinate_line(line, chain_mapping)
                hetatm_lines.append(line if line.endswith("\n") else line + "\n")
            elif line.startswith("CONECT"):
                conect_lines.append(line if line.endswith("\n") else line + "\n")

    return hetatm_lines, conect_lines


def write_final_pdb(
    structure,
    metadata_lines: List[str],
    original_pdb: str,
    out_pdb: str,
    original_chain_mapping: Optional[Dict[str, str]] = None,
) -> None:
    """
    Write generated protein ATOM records, then restore original HETATM/CONECT.

    Category chains therefore contain only ATOM records copied from chain A,
    while original hetero atoms and connectivity records appear exactly once.
    """
    pdb_io = PDBIO()
    pdb_io.set_structure(structure)

    atom_buffer = StringIO()
    pdb_io.save(
        atom_buffer,
        select=AtomRecordsOnlySelect(),
        write_end=False,
        preserve_atom_numbering=True,
    )

    hetatm_lines, conect_lines = read_original_hetatm_and_conect_records(
        original_pdb,
        chain_mapping=original_chain_mapping,
    )

    atom_lines = atom_buffer.getvalue().splitlines(keepends=True)

    # PDBIO creates TER serials from the preceding ATOM serial. When an
    # original HETATM occupies that next serial, this can create a duplicate.
    # Move generated TER serials above every ATOM/HETATM serial in the output.
    coordinate_serials: List[int] = []
    for line in atom_lines + hetatm_lines:
        if line.startswith(("ATOM  ", "HETATM")):
            try:
                coordinate_serials.append(int(line[6:11]))
            except ValueError:
                pass
    next_ter_serial = max(coordinate_serials, default=0) + 1
    for index, line in enumerate(atom_lines):
        if line.startswith("TER"):
            atom_lines[index] = line[:6] + f"{next_ter_serial:5d}" + line[11:]
            next_ter_serial += 1

    with open(out_pdb, "w") as out:
        out.writelines(metadata_lines)
        out.writelines(atom_lines)
        out.writelines(hetatm_lines)
        out.writelines(conect_lines)
        out.write("END\n")

# ============================================================
# Main pipeline (ancestry categories)
# ============================================================

def run_pipeline_categories(input_tsv: str, input_pdb: str, output_dir: str, categories: List[str]) -> str:
    os.makedirs(output_dir, exist_ok=True)
    start_time = datetime.now()

    def log(msg: str) -> None:
        now = datetime.now()
        elapsed = now - start_time
        fprint(f"[{now.strftime('%H:%M:%S')}] {msg} (t+{elapsed})")

    uniprot_id = os.path.basename(input_tsv).split("_")[0]
    final_out = os.path.join(output_dir, f"{uniprot_id}_category_gnomad.pdb")

    temp_paths: List[str] = []
    chain_labels: List[str] = CHAIN_LABELS[: len(categories)]
    reserved = set(chain_labels)

    try:
        log("1) TSV -> per-category top variant per position")
        df = filter_tsv_file(input_tsv)
        top_category_data = calculate_top_frequencies_per_category(df, categories)
        log(f" - completed TSV processing (categories: {len(categories)})")

        log("2) WT residue map from chain A (original PDB)")
        chainA_res_map = get_chain_resname_map(input_pdb, "A")

        log("3) Normalize original PDB (free reserved chain ids + normalize residue names)")
        normalized_pdb = os.path.join(output_dir, f"{uniprot_id}__tmp_normalized_original.pdb")
        mapping = normalize_pdb_chains(
            pdb_in=input_pdb,
            pdb_out=normalized_pdb,
            reserved=reserved,
            normalize_residues=True,
        )
        temp_paths.append(normalized_pdb)
        if mapping:
            log(f" - renamed original chains to free reserved ids: {mapping}")

        log("4) Create + mutate category chains using pyMUT (single-chain temp PDBs)")
        variant_chain_pdbs: List[Tuple[str, str, str]] = []  # (pdb_path, chain_label, category)

        for category, chain_label in zip(categories, chain_labels):
            chain_start = datetime.now()
            log(f"--- Category '{category}' -> Chain '{chain_label}' ---")

            # base file contains ONLY the one copied chain (chain_label)
            base = os.path.join(output_dir, f"{uniprot_id}__tmp_base_{chain_label}.pdb")
            extract_chainA_as_chain(
                pdb_in=normalized_pdb,
                pdb_out=base,
                chain_dst=chain_label,
                chain_src="A",
            )
            temp_paths.append(base)

            category_data = top_category_data.get(category, {})
            muts = plan_mutations_for_category(category_data, chainA_res_map)

            if not muts:
                log(f" - {chain_label}: no mutations (WT copy)")
                variant_chain_pdbs.append((base, chain_label, category))
                continue

            mut_out = os.path.join(output_dir, f"{uniprot_id}__tmp_mut_{chain_label}.pdb")
            apply_mutations_with_pymut(
                pdb_in=base,
                pdb_out=mut_out,
                chain_id=chain_label,
                mutations=muts,
            )
            temp_paths.append(mut_out)

            log(
                f" - {chain_label}: mutated {len(muts)} sites "
                f"(pyMUT={PYMUT_MUTATION_TYPE}; fallback={PYMUT_FALLBACK_MUTATION_TYPE}; "
                f"hard<{CLASH_HARD_CUTOFF_A:.2f}Å; neighbor<{CLASH_NEIGHBOR_CUTOFF_A:.1f}Å) "
                f"in {datetime.now() - chain_start}"
            )
            variant_chain_pdbs.append((mut_out, chain_label, category))

        log("5) Merge: start from normalized original (ALL original chains) + set chain A bfactors + add category chains")
        # IMPORTANT: merge function must:
        #  - keep all original chains from normalized_pdb
        #  - set chain A B-factors = max-bin across categories per residue
        #  - add/replace U/V/W/X/Y/... chains with category B-factors
        structure = merge_normalized_original_with_variant_chains(
            normalized_original_pdb=normalized_pdb,
            variant_chain_pdbs=variant_chain_pdbs,
            top_category_data=top_category_data,
        )

        log("6) Build ALL metadata/remarks (mapping + sequences) + write final PDB")
        metadata = build_metadata_lines_with_category_mapping_and_sequences(
            original_pdb=input_pdb,
            structure=structure,
            categories=categories,
            chain_labels=chain_labels,
            chains_to_report=set(["A"] + chain_labels),
        )
        write_final_pdb(
            structure,
            metadata,
            original_pdb=input_pdb,
            out_pdb=final_out,
            original_chain_mapping=mapping,
        )

    finally:
        for p in temp_paths:
            try:
                os.remove(p)
            except OSError:
                pass


# ============================================================
# CLI
# ============================================================

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python pymut_gnomad_ancestry_2_structure.py <input_tsv_file> <input_pdb_file> <output_directory> <categories>")
        print("Example: python pymut_gnomad_ancestry_2_structure.py input.tsv input.pdb output/ nfe,fin,amr")
        sys.exit(1)

    input_file = sys.argv[1]
    pdb_file_path = sys.argv[2]
    output_dir = sys.argv[3]
    categories_arg = sys.argv[4]
    categories = [c.strip() for c in categories_arg.split(",") if c.strip()]

    if not categories:
        raise SystemExit("Error: categories list is empty.")

    if len(categories) > len(CHAIN_LABELS):
        raise SystemExit(f"Error: too many categories ({len(categories)}). Max supported: {len(CHAIN_LABELS)}")

    run_pipeline_categories(input_file, pdb_file_path, output_dir, categories)