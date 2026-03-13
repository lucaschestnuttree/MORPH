# Standard library imports
import os
import re
import sys
import tempfile 
import numpy as np
import pandas as pd
from datetime import datetime
from collections import defaultdict
from Bio.PDB import PDBParser, PDBIO, Select, PPBuilder

# PyRosetta imports
from pyrosetta import init, pose_from_pdb
from pyrosetta.toolbox import mutate_residue
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.protocols.backrub import BackrubMover
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory, fa_rep
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
from pyrosetta.rosetta.core.pack.task.residue_selector import ClashBasedShellSelector
from pyrosetta.rosetta.core.pack.task.operation import RestrictToRepacking, IncludeCurrent

def filter_tsv_file(file_path: str) -> pd.DataFrame:
    """
    Filters rows in a TSV file based on specific criteria, including allele_count,
    protein_variant, allele_frequency, intron, uniprot_id, and am_pathogenicity.
    """
    # Load the TSV file
    df = pd.read_csv(file_path, sep='\t')

    # Add missing columns if necessary
    if 'intron' not in df.columns:
        df['intron'] = np.nan

    # Define min and max frequency for filtering
    min_frequency, max_frequency = 1e-9, 1.0

    # Apply filtering criteria
    df_filtered = df[
        (df['allele_count'].notna()) & 
        (df['allele_count'] != 0.0) & 
        (df['allele_frequency'] >= min_frequency) & 
        (df['allele_frequency'] <= max_frequency) &
        (df['uniprot_id'] != 'splice') &
        (~df['am_class'].isin(['new_stop', 'synonymous'])) & 
        (~df['protein_variant'].str.contains(r'[A-Z]1[A-Z]', na=False))
    ]

    return df_filtered

def extract_position(variant: str) -> int:
    """
    Extracts integer position from a protein variant string.
    """
    match = re.search(r'\d+', variant)
    if match:
        return int(match.group())
    else:
        raise ValueError(f"Invalid variant format: {variant}")

def map_frequency_to_bin(frequency: float) -> int:
    """
    Maps allele frequency to an integer from 0 to 8 based on defined bins.
    """
    frequency_bins = [
        (1e-9, 1e-8, 0),
        (1e-8, 1e-7, 1),
        (1e-7, 1e-6, 2),
        (1e-6, 1e-5, 3),
        (1e-5, 1e-4, 4),
        (1e-4, 1e-3, 5),
        (1e-3, 1e-2, 6),
        (1e-2, 0.51, 7),
        (0.51, 1.0, 8)
    ]

    for lower, upper, value in frequency_bins:
        if lower <= frequency < upper:
            return value
    return 0  # Default to 0 if frequency does not match any bin

def calculate_top_3_positions(df: pd.DataFrame) -> dict:
    """Identifies the top 3 highest allele frequencies and their last characters for each position."""
    position_values = defaultdict(list)

    # Gather frequencies and variants by position
    for _, row in df.iterrows():
        try:
            position = extract_position(row['protein_variant'])
            allele_frequency = float(row['allele_frequency'])
            if not pd.isna(allele_frequency):
                position_values[position].append((allele_frequency, row['protein_variant']))
        except ValueError as e:
            print(f"Skipping invalid row: {e}")

    top_3_data = {}
    for position, values in position_values.items():
        # Sort frequencies descending
        sorted_values = sorted(values, key=lambda x: x[0], reverse=True)[:3]

        # Extract the top 3 frequencies and last characters
        frequencies = [map_frequency_to_bin(freq) for freq, _ in sorted_values]
        last_characters = [variant[-1] for _, variant in sorted_values]

        # Pad to ensure three entries
        frequencies += [None] * (3 - len(frequencies))
        last_characters += [None] * (3 - len(last_characters))

        top_3_data[position] = {
            "frequencies": frequencies,
            "last_characters": last_characters,
        }

    return top_3_data

def mutate_residues(pose, mutated_residues: dict) -> None:
    """Mutates residues in a Pose object according to `mutated_residues`."""
    print("Performing residue mutations...")

    for position, target_residue in mutated_residues.items():
        try:
            rosetta_res_index = pose.pdb_info().pdb2pose('A', position)
            if rosetta_res_index == 0:
                print(f"Residue {position} not found in chain A. Skipping...")
                continue
            mutate_residue(pose, rosetta_res_index, target_residue)
        except Exception as e:
            print(f"Error mutating residue at position {position}: {e}")

def repack_rotamers(pose, clash_threshold: float = 20.0) -> None:
    """
    Repacks all rotamers and minimizes side chains for all residues in memory, 
    with a customizable clash penalty (fa_rep weight).

    Args:
        pose: The Rosetta Pose object to modify.
        clash_threshold: The weight for fa_rep in the scoring function (default: 20.0).
    """
    print("Repacking all rotamers...")

    # Task factory for repacking
    task_factory = TaskFactory()
    task_factory.push_back(RestrictToRepacking())
    task_factory.push_back(IncludeCurrent())

    pack_mover = PackRotamersMover()
    pack_mover.task_factory(task_factory)

    # Create a score function with customizable clash penalties
    scorefxn = ScoreFunctionFactory.create_score_function("ref2015")
    scorefxn.set_weight(fa_rep, clash_threshold)  # Use the provided clash threshold for fa_rep
    pack_mover.score_function(scorefxn)
    pack_mover.apply(pose)  # Apply the repacking

    # Minimize the pose side chains
    print("Minimizing side chains...")
    move_map = MoveMap()
    move_map.set_chi(True)  # Allow side-chain flexibility
    move_map.set_bb(False)  # Keep backbone rigid

    min_mover = MinMover()
    min_mover.movemap(move_map)
    min_mover.score_function(scorefxn)
    min_mover.min_type("dfpmin")  # Use quasi-Newton optimization method
    min_mover.apply(pose)

    print("Repacking and minimization complete.")

def minimize_clashes(pose, clash_threshold: float = 20.0) -> None:
    """
    Minimizes clashes in a Pose object using a customizable clash penalty.

    Args:
        pose: The Rosetta Pose object to operate on.
        clash_threshold: The weight for fa_rep in the scoring function (default: 20.0).
    """
    print("Identifying and minimizing clashes in memory...")

    # Scoring function with customizable clash penalties
    scorefxn = ScoreFunctionFactory.create_score_function("ref2015")
    scorefxn.set_weight(fa_rep, clash_threshold)  # Use the provided clash threshold for fa_rep

    # Define a ClashBasedShellSelector
    selector = ClashBasedShellSelector()
    selector.set_scorefxn(scorefxn)

    print("Identifying clash shell...")
    selected_residues = selector.apply(pose)  # Returns a boolean mask
    clash_residues = [i + 1 for i, selected in enumerate(selected_residues) if selected]

    print(f"Clashing residues identified: {clash_residues}")

    # Perform local adjustments for clashing residues
    if clash_residues:
        print("Applying BackrubMover for clashes...")
        backrub_mover = BackrubMover()
        move_map = MoveMap()

        for res_idx in clash_residues:
            move_map.set_bb(res_idx, True)  # Allow backbone movement
            move_map.set_chi(res_idx, True)  # Allow side-chain movement

        backrub_mover.set_movemap(move_map)
        backrub_mover.apply(pose)

        # Perform minimization
        print("Performing local minimization...")
        min_mover = MinMover()
        min_mover.movemap(move_map)
        min_mover.score_function(scorefxn)
        min_mover.min_type("dfpmin")  # Use quasi-Newton optimization method
        min_mover.apply(pose)
        print("Clash minimization and refinement complete.")
    else:
        print("No clashing residues identified. Skipping minimization.")

def merge_and_save_with_bfactors(
    pdb_file_path: str,            # Input file path for the original Chain A
    poses_with_chains: list,       # List of (Pose, chain_label)
    output_dir: str,
    uniprot_id: str,
    top_3_data: dict
) -> None:
    """
    Merges chains V, W, and X into the structure, applies B-factors to Chain A,
    updates metadata, removes hydrogens, and saves the final PDB file.
    """

    class NoHydrogensSelect(Select):
        """Exclude hydrogen atoms during PDB saving."""
        def accept_atom(self, atom):
            # Return False for hydrogen atoms, True otherwise
            return not atom.element == "H"
    
    # Load the original structure
    parser = PDBParser(QUIET=True)
    pdbio = PDBIO()

    # Load Chain A directly from original PDB file
    print("Loading Chain A directly from the original PDB file...")
    original_structure = parser.get_structure("original", pdb_file_path)
    merged_model = original_structure[0]  # First (and assumed only) model

    # Map of B-factors for residues in Chain V
    v_bfactor_map = {}

    # Process Chains V, W, and X, assign B-factors, and merge them
    print("Merging chains V, W, and X into the structure...")
    for pose, chain_label in poses_with_chains:
        with tempfile.NamedTemporaryFile(
            mode="w+", suffix=".pdb", prefix=f"Chain_{chain_label}_", delete=False
        ) as temp_file:
            pose.dump_pdb(temp_file.name)
            temp_file.seek(0)  # Reset pointer
            repacked_structure = parser.get_structure(chain_label, temp_file.name)

        # Iterate over model and chain
        for model in repacked_structure:
            for chain in model:
                chain.id = chain_label  # Relabel chain
                print(f"Updating B-factors for Chain {chain_label}...")

                # Assign B-factors to the chain
                chain_index = ["V", "W", "X"].index(chain_label)
                for residue in chain:
                    residue_number = residue.id[1]  # Extract residue number
                    if residue_number in top_3_data:
                        frequency_bin = top_3_data[residue_number]["frequencies"][chain_index]
                        bfactor = frequency_bin if frequency_bin is not None else 42.0
                        for atom in residue:
                            atom.bfactor = bfactor
                        # Save B-factors for Chain V to apply to Chain A later
                        if chain_label == "V":
                            v_bfactor_map[residue_number] = bfactor
                    else:
                        for atom in residue:
                            atom.bfactor = 42.0  # Default fallback B-factor

                merged_model.add(chain)  # Add chain to merged model

    # Update Chain A with B-factors from Chain V
    print("Applying B-factors from Chain V to Chain A...")
    for chain in merged_model:
        if chain.id == "A":  # Modify Chain A
            for residue in chain:
                residue_number = residue.id[1]
                if residue_number in v_bfactor_map:
                    for atom in residue:
                        atom.bfactor = v_bfactor_map[residue_number]
                else:
                    for atom in residue:
                        atom.bfactor = 42.0  # Default fallback B-factor

    # Preserve metadata and structure
    print("Reading original metadata from the PDB file...")
    with open(pdb_file_path, 'r') as pdb_file:
        lines = pdb_file.readlines()

    # Separate metadata and structure
    metadata_lines = [line for line in lines if not line.startswith(("ATOM", "HETATM", "ANISOU", "MODEL", "END"))]

    # Update HEADER with today's date
    today_date = datetime.now().strftime("%d-%b-%y").upper()
    metadata_lines = [
        f"HEADER GNOMAD ALLELE FREQUENCY                                    {today_date}\n"
        if line.startswith("HEADER") else line
        for line in metadata_lines
    ]

    # Add REMARK 2 metadata lines
    remark_2_lines = [
        "REMARK   2\n",
        "REMARK   2 REFERENCE 2\n",
        "REMARK   2  AUTH   LUCAS CHATAIGNER, OLIVIER SPERANDIO & HANNA JULIENNE\n",
        "REMARK   2  TITL   GNOMAD ALLELE FREQUENCY STRUCTURE MAPPING REVEALS AMAZINGNESS\n",
        "REMARK   2\n",
        "REMARK   2 DISCLAIMERS 2\n",
        "REMARK   2 PROVIDED HERE IS A HEAVILY ANNOTATED FILE FOR VISUALIZATION PURPOSES\n",
        "REMARK   2\n",  
        "REMARK   2 CHAIN A, B-FACTORS ARE ARBITRARILY ASSIGNED ACCORDING TO A SCHEME USEFUL\n"
        "REMARK   2 FOR VISUALIZATION OF GNOMAD ALLELE FREQUENCY.\n",
        "REMARK   2\n",
        "REMARK   2 CHAINS V, W, AND X, ARE NOT REAL PROTEINS, AND MAY CONTAIN UNREALISTIC STRUCTURAL ARTIFACTS.\n",
        "REMARK   2 THEY ARE USED TO REPRESENT MOST COMMON VARIANTS AT EACH POSITION WHERE VARIATIONS ARE REPORTED.\n",
        "REMARK   2 PER VARIABLE POSITION, THE MOST FREQUENT VARIANT IS SHOWN ON CHAIN V, WHILE SECOND\n", 
        "REMARK   2 AND THIRD MOST FREQUENT VARIANTS ARE SHOWN ON CHAIN W AND X\n",
    ]

    # Insert REMARK 2 before REMARK 1
    first_remark_1_index = next(
        (i for i, line in enumerate(metadata_lines) if line.startswith("REMARK   1")),
        len(metadata_lines)  # Default to end if REMARK 1 not found
    )
    metadata_lines = metadata_lines[:first_remark_1_index] + remark_2_lines + metadata_lines[first_remark_1_index:]

    # Combine metadata and the updated structure and save it directly
    print("Saving final output file...")
    combined_output_path = os.path.join(output_dir, f"{uniprot_id}_gnomad_AF-F1-model_v4.pdb")

    # Save directly to the final file
    with open(combined_output_path, 'w') as output_file:
        output_file.writelines(metadata_lines)
        pdbio.set_structure(original_structure)
        pdbio.save(output_file, select=NoHydrogensSelect())  # Write atomic data
        
    # --- Cleanup: remove HETATM records from chains other than A ---
    with open(combined_output_path, "r") as f:
        lines = f.readlines()

    with open(combined_output_path, "w") as f:
        for line in lines:
            if line.startswith("HETATM"):
                chain_id = line[21]  # chain ID is column 22 in PDB format
                if chain_id != "A":
                    continue  # skip HETATM from non-A chains
            f.write(line)

def extend_remark_with_sequences(
    pdb_file_path: str,
    output_file_path: str
) -> None:
    """
    Extends REMARK 2 section with polypeptide sequences for all chains as one-letter codes.

    Args:
        pdb_file_path: Path to the input PDB file.
        output_file_path: Path to the output PDB file with extended REMARK 2.
    """
    print(f"Extending REMARK 2 with sequences for all chains...")

    # Parse the PDB file and extract polypeptide sequences
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("model", pdb_file_path)
    ppb = PPBuilder()

    # Collect sequences for each chain
    chain_sequences = {}
    for chain in structure[0]:  # Assume first model
        chain_id = chain.id
        peptides = ppb.build_peptides(chain)
        sequence = "".join(str(peptide.get_sequence()) for peptide in peptides)
        chain_sequences[chain_id] = sequence
        print(f"Extracted sequence for chain {chain_id}: {sequence}")

    # Read the original PDB file
    with open(pdb_file_path, 'r') as pdb_file:
        pdb_lines = pdb_file.readlines()

    # Extend REMARK 2 with sequences
    remark_2_lines = [
        "REMARK   2\n",
        "REMARK   2 POLYPEPTIDE SEQUENCES FOR ALL CHAINS:\n",
    ]
    for chain_id, sequence in chain_sequences.items():
        remark_2_lines.append(f"REMARK   2  CHAIN {chain_id}: {sequence}\n")

    # Find where to insert the extended REMARK 2
    first_remark_1_index = next(
        (i for i, line in enumerate(pdb_lines) if line.startswith("REMARK   1")),
        len(pdb_lines)  # Default to end if REMARK 1 not found
    )
    updated_pdb_lines = (
        pdb_lines[:first_remark_1_index] +
        remark_2_lines +
        pdb_lines[first_remark_1_index:]
    )

    # Write the updated lines to the output file
    with open(output_file_path, 'w') as output_file:
        output_file.writelines(updated_pdb_lines)

    print(f"Extended REMARK 2 with sequences saved to {output_file_path}.")

def process_file(input_file: str, pdb_file_path: str, output_dir: str) -> None:
    """
    Processes a TSV file and PDB file to create mutated, repacked, and minimized chains 
    and merges them into a final structure.
    """
    uniprot_id = os.path.basename(input_file).split('_')[0]

    print("Processing TSV file...")
    df_filtered = filter_tsv_file(input_file)

    print("Calculating top 3 allele frequencies by position...")
    top_3_data = calculate_top_3_positions(df_filtered)

    print("Initializing PyRosetta...")
    init(extra_options="-ex1 -ex2aro -use_input_sc -ignore_unrecognized_res")
    original_pose = pose_from_pdb(pdb_file_path)

    # Prepare in-memory processing for chains V, W, X
    poses_with_chains = []
    chain_labels = ['V', 'W', 'X']

    for i, chain_label in enumerate(chain_labels):
        mutated_residues = {
            pos: residues["last_characters"][i]
            for pos, residues in top_3_data.items()
            if i < len(residues["last_characters"]) and residues["last_characters"][i] is not None
        }

        if mutated_residues:
            print(f"Processing chain {chain_label}...")

            # Clone Pose, mutate residues, repack, and minimize
            chain_pose = original_pose.clone()
            mutate_residues(chain_pose, mutated_residues)
            repack_rotamers(chain_pose)
            minimize_clashes(chain_pose)

            # Add processed pose
            poses_with_chains.append((chain_pose, chain_label))
        else:
            print(f"No mutation data available for chain {chain_label}, skipping...")

    print("Merging and saving final structure...")
    merged_file_path = os.path.join(output_dir, f"{uniprot_id}_gnomad_AF-F1-model_v4.pdb")
    merge_and_save_with_bfactors(
        pdb_file_path, poses_with_chains, output_dir, uniprot_id, top_3_data
    )

    print("Extending REMARK 2 with sequences...")
    extend_remark_with_sequences(
        pdb_file_path=merged_file_path,
        output_file_path=merged_file_path,  # Overwrite the same file
    )

    print(f"Process complete. Final file saved at {merged_file_path}.")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <input_tsv_file> <input_pdb_file> <output_directory>")
        sys.exit(1)

    input_file = sys.argv[1]
    pdb_file_path = sys.argv[2]
    output_dir = sys.argv[3]

    os.makedirs(output_dir, exist_ok=True)
    process_file(input_file, pdb_file_path, output_dir)