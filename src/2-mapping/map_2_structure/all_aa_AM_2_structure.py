import re
import sys
import csv
from collections import defaultdict
from Bio.PDB import PDBParser, PDBIO

def extract_position(variant: str) -> int:
    """
    Extract the integer position from a protein variant string (e.g., 'M1A' -> 1, 'L3124Y' -> 3124).
    """
    match = re.search(r'\d+', variant)
    if match:
        return int(match.group())
    else:
        raise ValueError(f"Invalid variant format: {variant}")

def calculate_position_averages(tsv_file_path: str) -> dict:
    """
    Reads a TSV file, calculates the average am_pathogenicity values for rows sharing the same protein variant position,
    and returns a dictionary mapping residue number (position) to average values.
    """
    position_values = defaultdict(list)  # Dictionary to store am_pathogenicity values by position
    
    # Read the TSV and group am_pathogenicity values by position
    with open(tsv_file_path, newline='') as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')
        
        # Check if the necessary columns exist
        if 'protein_variant' not in reader.fieldnames or 'am_pathogenicity' not in reader.fieldnames:
            raise ValueError(f"Required columns 'protein_variant' or 'am_pathogenicity' not found in the file.")
        
        for row in reader:
            try:
                position = extract_position(row['protein_variant'])  # Extract the position (residue number)
                am_pathogenicity = float(row['am_pathogenicity'])  # Convert am_pathogenicity to float
                position_values[position].append(am_pathogenicity)
            except ValueError as e:
                print(f"Skipping invalid row: {e}")
    
    # Calculate the average for each position
    position_avg = {position: round(sum(values) / len(values), 2)  # Round to 2 decimal places
                    for position, values in position_values.items()}
    
    return position_avg


def update_pdb_bfactor(pdb_file_path: str, output_file_path: str, position_avg: dict) -> None:
    """
    Updates the B-factor values in the PDB file based on residue number and the provided pathogenicity scores.
    Ensures that HETATM residues are not modified and that CONECT records are preserved.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("PDB", pdb_file_path)

    updated_lines = []  # Store modified PDB lines

    # Read the PDB file and modify B-factors only for ATOM records
    with open(pdb_file_path, "r") as pdb_file:
        for line in pdb_file:
            record_type = line[:6].strip()  # Extract record type (ATOM, HETATM, CONECT, etc.)

            if record_type == "ATOM":  # Modify B-factors only for ATOM records
                try:
                    residue_number = int(line[22:26].strip())  # Extract residue number
                    
                    if residue_number in position_avg:
                        new_bfactor = f"{position_avg[residue_number]:6.2f}"  # Format B-factor
                        modified_line = line[:60] + new_bfactor.rjust(6) + line[66:]  # Replace B-factor
                        updated_lines.append(modified_line)
                        continue  # Skip to next line after modification

                except ValueError:
                    pass  # If residue number parsing fails, keep the line unchanged

            # Preserve all other records as they are (CONECT, HETATM, etc.)
            updated_lines.append(line)

    # Write updated PDB file
    with open(output_file_path, "w") as output_pdb:
        output_pdb.writelines(updated_lines)

    print(f"B-factors updated while preserving CONECT records in {output_file_path}")


if __name__ == "__main__":
    # Expect: script.py <tsv_file> <input_pdb> <output_pdb>
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <variants.tsv> <input.pdb> <output.pdb>")
        sys.exit(1)

    tsv_file_path = sys.argv[1]
    pdb_file_path = sys.argv[2]
    output_file_path = sys.argv[3]

    # Step 1: Calculate the average pathogenicity scores per position from the TSV file
    position_avg = calculate_position_averages(tsv_file_path)
    
    # Step 2: Update the B-factors in the PDB file using the calculated averages
    update_pdb_bfactor(pdb_file_path, output_file_path, position_avg)