from __future__ import annotations

import os
import re
import sys
import functools
from copy import deepcopy
from datetime import datetime
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional, Set

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
# HARD-CODED KNOBS
# ============================================================

PYMUT_MUTATION_TYPE: str = "first"
PYMUT_FALLBACK_MUTATION_TYPE: str = "best"

MAX_FALLBACK_REDO: int = 10000  # max number of residues to re-rotamerize

# Distance-cutoff clash detector (Å), SAME CHAIN ONLY
CLASH_NEIGHBOR_CUTOFF_A: float = 8.0   # search radius (performance)
CLASH_HARD_CUTOFF_A: float = 3       # "clash" if any atom pair closer than this

# Other
DEFAULT_BFACTOR = 42.0
CHAIN_ORDER = ["V", "W", "X"]
SEQ_REMARK_CHAINS = {"A", "V", "W", "X"}
BACKBONE_NAMES = {"N", "CA", "C", "O", "OXT"}

# ============================================================
# Constants (3 lines)
# ============================================================

AA1_TO_AA3: Dict[str, str] = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE",
    "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}
AA3_NORMALIZE: Dict[str, str] = {"HID": "HIS", "HIE": "HIS", "HIP": "HIS", "CYX": "CYS", "CYN": "CYS", "CSS": "CYS", "CYM": "CYS", "MSE": "MET"}

# ============================================================
# Data structures
# ============================================================

@dataclass(frozen=True)
class Mutation:
    wt3: str
    pos: int
    mut3: str


@dataclass(frozen=True)
class Top3Entry:
    freq_bins: List[Optional[int]]  # len 3
    mut_aa1: List[Optional[str]]    # len 3


Top3Data = Dict[int, Top3Entry]


class NoHydrogensSelect(Select):
    def accept_atom(self, atom):
        return atom.element != "H"


# ============================================================
# TSV processing
# ============================================================

def filter_tsv_file(file_path: str) -> pd.DataFrame:
    df = pd.read_csv(file_path, sep="\t")
    if "intron" not in df.columns:
        df["intron"] = np.nan
    min_frequency, max_frequency = 1e-9, 1.0
    return df[
        (df["allele_count"].notna())
        & (df["allele_count"] != 0.0)
        & (df["allele_frequency"] >= min_frequency)
        & (df["allele_frequency"] <= max_frequency)
        & (~df["am_class"].isin(["new_stop", "synonymous", "splice"]))
        & (~df["protein_variant"].str.contains(r"[A-Z]1[A-Z]", na=False))
    ]


def extract_position(variant: str) -> int:
    m = re.search(r"\d+", variant)
    if not m:
        raise ValueError(f"Invalid variant format: {variant}")
    return int(m.group())


def map_frequency_to_bin(frequency: float) -> int:
    bins = [
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
    for lo, hi, val in bins:
        if lo <= frequency < hi:
            return val
    return 0


def calculate_top_3_positions(df: pd.DataFrame) -> Top3Data:
    position_values: Dict[int, List[Tuple[float, str]]] = defaultdict(list)
    for _, row in df.iterrows():
        try:
            pos = extract_position(row["protein_variant"])
            af = float(row["allele_frequency"])
            if not pd.isna(af):
                position_values[pos].append((af, row["protein_variant"]))
        except ValueError as e:
            fprint(f"Skipping invalid row: {e}")

    top3: Top3Data = {}
    for pos, values in position_values.items():
        sorted_values = sorted(values, key=lambda x: x[0], reverse=True)[:3]
        freq_bins = [map_frequency_to_bin(freq) for freq, _ in sorted_values]
        mut_aa1 = [variant[-1] for _, variant in sorted_values]
        freq_bins += [None] * (3 - len(freq_bins))
        mut_aa1 += [None] * (3 - len(mut_aa1))
        top3[pos] = Top3Entry(freq_bins=freq_bins, mut_aa1=mut_aa1)
    return top3


# ============================================================
# Chain handling (Bio.PDB)
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
    io.save(pdb_out)
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
    io.save(pdb_out)


# ============================================================
# Mutation planning
# ============================================================

def plan_mutations_for_rank(top3: Top3Data, chainA_res_map: Dict[int, str], rank: int) -> List[Mutation]:
    muts: List[Mutation] = []
    for pos, entry in top3.items():
        aa1 = entry.mut_aa1[rank]
        if aa1 is None:
            continue
        wt3 = chainA_res_map.get(pos)
        if wt3 is None:
            continue
        mut3 = AA1_TO_AA3.get(aa1.upper())
        if mut3 is None:
            continue
        if wt3 == mut3:
            continue
        muts.append(Mutation(wt3=wt3, pos=pos, mut3=mut3))
    return muts


# ============================================================
# FAST clash helpers (fixes issues 1 & 2)
# ============================================================

def build_same_chain_neighborsearch(structure, chain_id: str):
    """
    Build NeighborSearch ONCE for the given chain.
    Returns (chain, NeighborSearch) so residue lookup is cheap.
    """
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
    """
    One pass over neighbors:
      - has_clash: True if ANY sidechain-heavy atom is within hard_cutoff Å of ANY heavy atom
        in other residues of the SAME chain.
      - min_dist: minimum such heavy-atom distance found (Å), for ranking; inf if none found.
    """
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
# pyMUT mutation + re-rotamerization (distance cutoff)
# ============================================================

def apply_mutations_with_pymut(
    pdb_in: str,
    pdb_out: str,
    chain_id: str,
    mutations: List[Mutation],
) -> None:
    """
    Phase 1: mutate all sites using PYMUT_MUTATION_TYPE.
    Phase 2: build NeighborSearch ONCE; one-pass score (has_clash + min_dist) for each mutated residue;
             re-rotamerize up to MAX_FALLBACK_REDO worst (smallest min_dist) using PYMUT_FALLBACK_MUTATION_TYPE.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("mut", pdb_in)

    try:
        rotamer_lib = load_rotamers()
    except Exception as e:
        rotamer_lib = None
        fprint(f"Warning: could not preload rotamer library (fallback to per-mutation load): {e}")

    mutated_positions: List[int] = []
    desired_res_by_pos: Dict[int, str] = {}

    # -------------------------------
    # Phase 1
    # -------------------------------
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

    # -------------------------------
    # Phase 2 (FAST): build NS once + one-pass score
    # -------------------------------
    try:
        chain, ns = build_same_chain_neighborsearch(structure, chain_id)
    except Exception as e:
        fprint(f"Warning: could not build NeighborSearch for chain {chain_id}: {e}")
        io = PDBIO()
        io.set_structure(structure)
        io.save(pdb_out)
        return

    clashing: List[Tuple[float, int]] = []  # (min_dist, pos)
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

    clashing.sort(key=lambda t: t[0])  # smallest distance = worst
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
            # Rebuild NeighborSearch after coordinate changes (cheap; done only for re-mutated subset)
            chain, ns = build_same_chain_neighborsearch(structure, chain_id)

            try:
                residue = chain[pos]
            except Exception:
                residue = None

            if residue is None:
                fprint(f" -> {chain_id}:{pos} fallback {PYMUT_FALLBACK_MUTATION_TYPE} (pre min={pre_min_dist:.2f}Å, post=NA)")
            else:
                _, post_min_dist = clash_and_min_distance_same_chain(
                    ns, residue,
                    neighbor_cutoff=CLASH_NEIGHBOR_CUTOFF_A,
                    hard_cutoff=CLASH_HARD_CUTOFF_A,
                )
                if post_min_dist == float("inf"):
                    fprint(f" -> {chain_id}:{pos} fallback {PYMUT_FALLBACK_MUTATION_TYPE} (pre min={pre_min_dist:.2f}Å, post=NA)")
                else:
                    fprint(
                        f" -> {chain_id}:{pos} fallback {PYMUT_FALLBACK_MUTATION_TYPE} "
                        f"(pre min={pre_min_dist:.2f}Å, post min={post_min_dist:.2f}Å)"
                    )

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
# B-factor annotation + merge
# ============================================================

def annotate_chain_bfactors_use_V(chain, top3: Top3Data) -> Dict[int, float]:
    """
    Assign B-factors to this chain using *V's* bin values (rank-1) for every residue.
    This makes A/V/W/X all share the same bins.
    """
    bmap: Dict[int, float] = {}
    for residue in chain:
        resnum = residue.id[1]

        if resnum in top3:
            # Always use rank-1 bin (chain V's bin)
            bin_val = top3[resnum].freq_bins[0]
            b = float(bin_val) if bin_val is not None else DEFAULT_BFACTOR
        else:
            b = DEFAULT_BFACTOR

        for atom in residue:
            atom.bfactor = b
        bmap[resnum] = b

    return bmap

def merge_normalized_original_with_variant_chains(
    normalized_original_pdb: str,
    variant_chain_pdbs: List[Tuple[str, str]],
    top3: Top3Data,
) -> Tuple[object, Dict[int, float]]:
    """
    Merge normalized original + V/W/X chains,
    and set B-factors on *all* of A/V/W/X using V's bin values (rank-1).
    """
    original_structure = parse_structure(normalized_original_pdb, "original_norm")
    model = original_structure[0]

    # Add/replace variant chains, but B-factor them all with V bins
    for pdb_path, chain_label in variant_chain_pdbs:
        s = parse_structure(pdb_path, chain_label)
        m = s[0]
        if chain_label not in [c.id for c in m]:
            raise ValueError(f"Expected chain {chain_label} not found in {pdb_path}")

        new_chain = deepcopy(m[chain_label])

        # Overwrite B-factors using V-bin mapping
        annotate_chain_bfactors_use_V(new_chain, top3)

        if chain_label in [c.id for c in model]:
            model.detach_child(chain_label)
        model.add(new_chain)

    # ALSO overwrite chain A B-factors using V bins
    if "A" in [c.id for c in model]:
        annotate_chain_bfactors_use_V(model["A"], top3)

    # Keep return signature stable; build a V-map (optional)
    v_bfactor_map: Dict[int, float] = {}
    if "A" in [c.id for c in model]:
        # map from chain A now equals V bins, so it's fine as "v map"
        for residue in model["A"]:
            v_bfactor_map[residue.id[1]] = float(residue.child_list[0].bfactor) if residue.child_list else DEFAULT_BFACTOR

    return original_structure, v_bfactor_map

# ============================================================
# Metadata + writing
# ============================================================

def build_metadata_lines_from_structure(original_pdb: str, structure) -> List[str]:
    with open(original_pdb, "r") as f:
        lines = f.readlines()

    drop_prefixes = ("ATOM", "HETATM", "ANISOU", "MODEL", "END", "ENDMDL", "TER", "CONECT", "MASTER")
    metadata = [ln for ln in lines if not ln.startswith(drop_prefixes)]

    today_date = datetime.now().strftime("%d-%b-%y").upper()
    metadata = [f"HEADER GNOMAD ALLELE FREQUENCY {today_date}\n" if ln.startswith("HEADER") else ln for ln in metadata]

    remark_1_lines = [
        "REMARK 1\n", "REMARK 1 REFERENCE\n", 
        "REMARK 1 AUTH LUCAS CHATAIGNER & HANNA JULIENNE\n",
        "REMARK 1 TITL GNOMAD ALLELE FREQUENCY STRUCTURE MAPPING\n", 
        "REMARK 1\n", 
        "REMARK 1 DISCLAIMERS\n",
        "REMARK 1 PROVIDED HERE IS A HEAVILY ANNOTATED FILE FOR VISUALIZATION PURPOSES\n", 
        "REMARK 1\n",
        "REMARK 1 B-FACTORS ARE ARBITRARILY ASSIGNED ACCORDING TO A SCHEME USEFUL\n",
        "REMARK 1 FOR VISUALIZATION OF GNOMAD ALLELE FREQUENCY.\n",
        "REMARK 1\n",
        "REMARK 1 CHAINS V, W, AND X, ARE NOT REAL PROTEINS AND CONTAIN UNREALISTIC STRUCTURAL ELEMENTS..\n",
        "REMARK 1 THEY ARE USED TO REPRESENT MOST COMMON VARIANTS AT EACH POSITION.\n",
        "REMARK 1 THE MOST FREQUENT VARIANT IS SHOWN ON CHAIN V, WHILE 2nd AND 3rd ON CHAIN W AND X\n",
        "REMARK 1\n",
    ]

    ppb = PPBuilder()
    seq_lines = ["REMARK 1 POLYPEPTIDE SEQUENCES FOR CHAINS A, V, W, X:\n"]

    model = structure[0]
    for chain in model:
        if chain.id not in SEQ_REMARK_CHAINS:
            continue
        peptides = ppb.build_peptides(chain)
        seq = "".join(str(p.get_sequence()) for p in peptides)
        seq_lines.append(f"REMARK 1 CHAIN {chain.id}: {seq}\n")

    first_remark_1_index = next((i for i, ln in enumerate(metadata) if ln.startswith("REMARK 1")), len(metadata))
    return metadata[:first_remark_1_index] + remark_1_lines + seq_lines + metadata[first_remark_1_index:]


def write_final_pdb(structure, metadata_lines: List[str], out_pdb: str) -> None:
    io = PDBIO()
    io.set_structure(structure)
    with open(out_pdb, "w") as out:
        out.writelines(metadata_lines)
        io.save(out, select=NoHydrogensSelect())


# ============================================================
# Pipeline 
# ============================================================

def run_pipeline(input_tsv: str, input_pdb: str, output_dir: str) -> str:
    os.makedirs(output_dir, exist_ok=True)
    start_time = datetime.now()

    def log(msg: str) -> None:
        now = datetime.now()
        elapsed = now - start_time
        fprint(f"[{now.strftime('%H:%M:%S')}] {msg} (t+{elapsed})")

    uniprot_id = os.path.basename(input_tsv).split("_")[0]
    final_out = os.path.join(output_dir, f"{uniprot_id}_AF_gnomad.pdb")
    temp_paths: List[str] = []

    try:
        log("1) TSV -> top3")
        df = filter_tsv_file(input_tsv)
        top3 = calculate_top_3_positions(df)
        log(f" - completed TSV processing and top3 calc (positions: {len(top3)})")

        log("2) WT residue map from chain A (original PDB)")
        chainA_res_map = get_chain_resname_map(input_pdb, "A")

        log("3) Normalize original PDB (free V/W/X + normalize residue names)")
        normalized_pdb = os.path.join(output_dir, f"{uniprot_id}__tmp_normalized_original.pdb")
        chain_mapping = normalize_pdb_chains(
            pdb_in=input_pdb,
            pdb_out=normalized_pdb,
            reserved=set(CHAIN_ORDER),
            normalize_residues=True,
        )
        temp_paths.append(normalized_pdb)
        if chain_mapping:
            log(f" - renamed original chains: {chain_mapping}")

        log("4) Create + mutate V/W/X chains (fast clash detection)")
        variant_chain_pdbs: List[Tuple[str, str]] = []

        for rank, chain_label in enumerate(CHAIN_ORDER):
            chain_start = datetime.now()

            base = os.path.join(output_dir, f"{uniprot_id}__tmp_base_{chain_label}.pdb")
            extract_chainA_as_chain(
                pdb_in=normalized_pdb,
                pdb_out=base,
                chain_dst=chain_label,
                chain_src="A",
            )
            temp_paths.append(base)

            muts = plan_mutations_for_rank(top3, chainA_res_map, rank)

            if not muts:
                log(f" - {chain_label}: no mutations (WT copy)")
                variant_chain_pdbs.append((base, chain_label))
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
            variant_chain_pdbs.append((mut_out, chain_label))

        log("5) Merge normalized original + add V/W/X + annotate B-factors (use V bins for all chains)")
        # NOTE: this relies on your UPDATED merge function, which sets B-factors on A/V/W/X using V's bins
        structure, _ = merge_normalized_original_with_variant_chains(
            normalized_original_pdb=normalized_pdb,
            variant_chain_pdbs=variant_chain_pdbs,
            top3=top3,
        )

        log("6) Build metadata + write final PDB")
        metadata = build_metadata_lines_from_structure(input_pdb, structure)
        write_final_pdb(structure, metadata, final_out)

        log(f"Done. Final file: {final_out}")
        return final_out

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
    if len(sys.argv) != 4:
        fprint("Usage: python pymut_gnomad_af_2_structure.py <input_tsv_file> <input_pdb_file> <output_directory>")
        sys.exit(1)

    run_pipeline(sys.argv[1], sys.argv[2], sys.argv[3])