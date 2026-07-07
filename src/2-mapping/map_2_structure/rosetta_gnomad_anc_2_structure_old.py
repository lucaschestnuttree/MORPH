# Standard library imports
import os
import re
import sys
import tempfile 
import numpy as np
import pandas as pd
from datetime import datetime
from string import ascii_uppercase
from collections import defaultdict
from typing import Any, Dict, List, Tuple
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
        (df['intron'].isna() | df['intron'].isin([None, ''])) & 
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

def calculate_top_frequencies_per_category(
    df: pd.DataFrame, categories: List[str]
) -> Dict[str, Dict[int, Dict[str, Any]]]:
    """
    Calculates the maximum frequency per position for each category and maps it to the corresponding residue.
    """
    category_position_data = {category: {} for category in categories}

    for category in categories:
        # Derive the correct column name
        column_name = (
            "allele_frequency" if category == "allele_frequency" else f"{category}_allele_frequency"
        )
        if column_name not in df.columns:
            print(f"Warning: Column '{column_name}' not found in DataFrame. Skipping category '{category}'.")
            continue

        # Initialize category data as a defaultdict of dictionaries
        category_data = defaultdict(lambda: {"frequency": 0.0, "last_character": None})
        
        for _, row in df.iterrows():
            if pd.isna(row[column_name]) or row[column_name] <= 0:
                continue  # Skip invalid or zero frequencies

            try:
                # Extract position and mutation data
                position = extract_position(row["protein_variant"])
                frequency = row[column_name]
                last_character = row["protein_variant"][-1]

                # Ensure the structure at `position` is valid
                if not isinstance(category_data[position], dict):
                    raise TypeError(f"Expected dict for position {position}, got {type(category_data[position]).__name__}")

                # Update max frequency and corresponding last_character
                if frequency > category_data[position]["frequency"]:
                    category_data[position]["frequency"] = frequency
                    category_data[position]["last_character"] = last_character

            except ValueError as e:
                print(f"Skipping invalid row due to invalid variant format: {e}")
            except TypeError as e:
                print(f"TypeError encountered during row processing: {e}")
        
        # Store result for this category
        category_position_data[category] = dict(category_data)

    return category_position_data

def chainA_frequency(top_category_data: Dict[str, Dict[int, Dict[str, Any]]]) -> Dict[int, float]:
    """
    Creates a dictionary mapping each position to the maximum frequency across all categories.

    Args:
        top_category_data: Nested dictionary with frequency data for each category and position.
            Structure:
            {
                "category1": {
                    100: {"frequency": 0.05, "last_character": "A"},
                    150: {"frequency": 0.02, "last_character": "G"},
                },
                "category2": {
                    100: {"frequency": 0.1, "last_character": "T"},
                    150: {"frequency": 0.03, "last_character": "C"},
                },
            }

    Returns:
        A dictionary where the keys are positions and the values are the maximum frequencies across all categories.
        Structure:
        {
            100: 0.1,  # Maximum of frequencies from all categories for position 100
            150: 0.03, # Maximum of frequencies from all categories for position 150
        }
    """
    position_max_frequencies = {}

    # Loop through each category and its associated positions
    for category, positions in top_category_data.items():
        for position, data in positions.items():
            frequency = data.get("frequency", 0.0)  # Default to 0.0 if frequency is missing

            # Update the maximum frequency for each position
            if position not in position_max_frequencies:
                position_max_frequencies[position] = frequency
            else:
                position_max_frequencies[position] = max(position_max_frequencies[position], frequency)

    return position_max_frequencies

def mutate_residues(pose: Any, mutated_residues: Dict[int, str]) -> None:
    print("Performing residue mutations...")

    for position, target_residue in mutated_residues.items():
        try:
            if not isinstance(target_residue, str) or len(target_residue) != 1:
                print(f"Invalid target residue at position {position}: {target_residue}")
                continue

            rosetta_res_index = pose.pdb_info().pdb2pose('A', position)
            if rosetta_res_index == 0:
                print(f"Residue {position} not found in chain A. Skipping...")
                continue

            mutate_residue(pose, rosetta_res_index, target_residue)
            print(f"Mutated position {position} to {target_residue}")
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

def prepare_chains(pdb_file_path: str, top_category_data: Dict[str, Dict[int, Dict[str, Any]]], categories: List[str]) -> List[Tuple[Any, int]]:
    print("Initializing PyRosetta...")
    init(extra_options="-ex1 -ex2aro -use_input_sc -ignore_unrecognized_res")

    original_pose = pose_from_pdb(pdb_file_path)
    chain_labels = ascii_uppercase[21:25]  # V, W, X, Y

    poses_with_chains = []
    for chain_index, category in enumerate(categories[:len(chain_labels)]):
        if category not in top_category_data:
            print(f"No data for category '{category}', skipping...")
            continue

        category_data = top_category_data[category]
        mutated_residues = {
            pos: data.get("last_character")
            for pos, data in category_data.items()
            if data.get("last_character") is not None
        }

        if not mutated_residues:
            print(f"No mutation data for category '{category}', skipping...")
            continue

        print(f"Processing category '{category}' as Chain '{chain_labels[chain_index]}'...")
        chain_pose = original_pose.clone()

        mutate_residues(chain_pose, mutated_residues)
        repack_rotamers(chain_pose)
        minimize_clashes(chain_pose)

        poses_with_chains.append((chain_pose, chain_index))

    return poses_with_chains

def merge_and_save_with_bfactors(
    pdb_file_path: str,
    poses_with_chains: List[Tuple[Any, int]],
    output_dir: str,
    uniprot_id: str,
    top_category_data: Dict[str, Dict[int, Dict[str, Any]]],
    categories: List[str]
) -> None:
    """
    Merges mutated chains into the original PDB file, assigns B-factors based on maximum frequencies,
    updates Chain A's B-factors, and saves the final structure with metadata intact.

    Args:
        pdb_file_path (str): Path to the original PDB file containing chain A structure.
        poses_with_chains (List): Mutated Pose objects to be added to chain A.
        output_dir (str): Directory where the resulting structure is saved.
        uniprot_id (str): Protein identifier for naming output files.
        top_category_data (Dict): Dictionary with mutation data by category and positions.
        categories (List): List of categories (e.g., 'nfe', 'fin').
    """
    class NoHydrogensSelect(Select):
        """Exclude hydrogen atoms during PDB saving."""
        def accept_atom(self, atom):
            return not atom.element == "H"

    # Load the original structure
    parser = PDBParser(QUIET=True)
    pdbio = PDBIO()
    
    print("Loading original structure (Chain A only)...")
    original_structure = parser.get_structure("original", pdb_file_path)
    merged_model = original_structure[0]  # First and only model

    # Update Chain A's B-factors based on combined frequency data
    print("Updating B-factors for Chain A...")
    chain_a_frequencies = chainA_frequency(top_category_data)

    for chain in merged_model:
        if chain.id == "A":  # Ensure we're working with Chain A
            for residue in chain:
                residue_number = residue.id[1]  # Get the residue number
                if residue_number in chain_a_frequencies:
                    # Get the maximum frequency for this residue
                    frequency = chain_a_frequencies[residue_number]
                    # Map frequency to B-factor
                    bfactor = map_frequency_to_bin(frequency)
                    for atom in residue:
                        atom.bfactor = bfactor  # Update B-factor
                else:
                    # Assign a default B-factor if no frequency is found
                    for atom in residue:
                        atom.bfactor = 42.0  # Default fallback B-factor

    print("Merging additional chains into the structure...")
    # Dynamically assign chain letters V to Z (up to 5 categories)
    chain_labels = ascii_uppercase[21:25]  # V, W, X, Y

    for (pose, category_index), chain_label in zip(poses_with_chains, chain_labels):
        category = categories[category_index]

        print(f"Merging category '{category}', assigned as Chain '{chain_label}'...")

        with tempfile.NamedTemporaryFile(
            mode="w+", suffix=".pdb", prefix=f"Chain_{chain_label}_", delete=False
        ) as temp_file:
            pose.dump_pdb(temp_file.name)
            temp_file.seek(0)  # Reset the file pointer
            repacked_structure = parser.get_structure(chain_label, temp_file.name)

        # Process model
        for model in repacked_structure:
            for chain in model:
                chain.id = chain_label
                print(f"Processing chain {chain_label} for category '{category}'...")

                # Assign B-factors to residues
                for residue in chain:
                    residue_number = residue.id[1]
                    if residue_number in top_category_data[category]:
                        freq = top_category_data[category][residue_number]["frequency"]
                        bfactor = map_frequency_to_bin(freq)
                        for atom in residue:
                            atom.bfactor = bfactor
                    else:
                        for atom in residue:
                            atom.bfactor = 42.0  # Default fallback B-factor for unmapped residues

                # Merge this chain into the structure
                merged_model.add(chain)

    print("Preserving metadata and saving final output file...")
    # Read and process metadata from the original PDB file
    with open(pdb_file_path, 'r') as pdb_file:
        lines = pdb_file.readlines()

    # Separate metadata from atomic data
    metadata_lines = [line for line in lines if not line.startswith(("ATOM", "HETATM", "ANISOU", "MODEL", "END"))]

    # Update HEADER with today's date
    today_date = datetime.now().strftime("%d-%b-%y").upper()
    metadata_lines = [
        f"HEADER    GNOMAD CATEGORY MAPPING                          {today_date}\n"
        if line.startswith("HEADER") else line
        for line in metadata_lines
    ]

    # Add REMARK 2 metadata
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
        "REMARK   2 FOR VISUALIZATION OF MAX ALLELE FREQUENCY FROM SELECTED ANCESTRY OR SEX.\n",
        "REMARK   2\n",
        "REMARK   2 CHAINS V, W, AND X, ARE NOT REAL PROTEINS, AND MAY CONTAIN UNREALISTIC STRUCTURAL ARTIFACTS.\n",
        "REMARK   2 THEY ARE USED TO REPRESENT MOST COMMON VARIANTS FOR GIVEN ANCESTRY OR SEX AT EACH POSITION.\n",
    ]

    remark_2_lines += [
        f"REMARK   2 CHAIN {chain_labels[i]} CORRESPONDS TO CATEGORY: {categories[i].upper()}\n"
        for i in range(len(categories))
    ]

    # Insert REMARK 2 before REMARK 1
    first_remark_1_index = next(
        (i for i, line in enumerate(metadata_lines) if line.startswith("REMARK   1")),
        len(metadata_lines)  # Default to end if REMARK 1 not found
    )
    metadata_lines = metadata_lines[:first_remark_1_index] + remark_2_lines + metadata_lines[first_remark_1_index:]

    # Save final structure with metadata
    combined_output_path = os.path.join(output_dir, f"{uniprot_id}_gnomad_category-F1-model_v4.pdb")

    with open(combined_output_path, 'w') as output_file:
        output_file.writelines(metadata_lines)
        pdbio.set_structure(original_structure)
        pdbio.save(output_file, select=NoHydrogensSelect())  # Write atomic data

    print(f"Final PDB structure saved to: {combined_output_path}")
    
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

def process_file(input_file: str, pdb_file_path: str, output_dir: str, categories: List[str]) -> None:
    """
    Processes mutation data, prepares mutated chains for each category, and merges them into the final structure.

    Args:
        input_file (str): Path to the TSV file with mutation data.
        pdb_file_path (str): Path to the input PDB structure.
        output_dir (str): Directory where the output file will be saved.
        categories (List[str]): List of mutation categories to process.
    """
    uniprot_id = os.path.basename(input_file).split("_")[0]

    print("Processing the TSV file...")
    df_filtered = filter_tsv_file(input_file)

    print("Calculating top allele frequencies for categories...")
    top_category_data = calculate_top_frequencies_per_category(df_filtered, categories)

    print("Preparing mutated chains...")
    poses_with_chains = prepare_chains(pdb_file_path, top_category_data, categories)

    print("Merging and saving the final structure...")
    merged_file_path = os.path.join(output_dir, f"{uniprot_id}_gnomad_category-F1-model_v4.pdb")
    merge_and_save_with_bfactors(
        pdb_file_path, poses_with_chains, output_dir, uniprot_id, top_category_data, categories
    )

    print("Extending REMARK 2 with sequences...")
    extend_remark_with_sequences(
        pdb_file_path=merged_file_path,
        output_file_path=merged_file_path,  # Overwrite the same file
    )

    print(f"Process complete. Final file saved at {merged_file_path}.")   
    
if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script.py <input_tsv_file> <input_pdb_file> <output_directory> <categories>")
        print("Example: python script.py input.tsv input.pdb output/ nfe,fin,amr")
        sys.exit(1)

    # Parse arguments
    input_file = sys.argv[1]  # Path to the input TSV file
    pdb_file_path = sys.argv[2]  # Path to the input PDB file
    output_dir = sys.argv[3]  # Directory to save output files
    categories_arg = sys.argv[4]  # Comma-separated categories

    # Parse categories into a list
    categories = [category.strip() for category in categories_arg.split(",")]

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Call the main processing function with parsed arguments
    process_file(input_file, pdb_file_path, output_dir, categories)