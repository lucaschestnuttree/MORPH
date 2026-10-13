import os
import sys
import glob
import uuid 
import subprocess
import pandas as pd

##  Load & Transform data
def load_tsv(file_path: str) -> pd.DataFrame:
    """Load a TSV file into a pandas DataFrame."""
    try:
        df = pd.read_csv(file_path, sep='\t')
        print(f"Loaded data from {file_path}")

        # 🚨 Early exit if file has header only
        if df.empty:
            print(f"[INFO] Input file has header only (no rows): {file_path}")
            sys.exit(0)

        return df

    except FileNotFoundError:
        raise FileNotFoundError(f"The file {file_path} does not exist.")
    except pd.errors.ParserError:
        raise ValueError(f"Error parsing the file {file_path}.")

def transform_gnomad_data(df: pd.DataFrame) -> pd.DataFrame:
    cleaned = df['alleles'].astype(str).str.replace(r'[\[\]"]', '', regex=True)

    # Split first, then filter based on the result (more robust than count(',')==1)
    parts = cleaned.str.split(',', n=1, expand=True).reindex(columns=[0, 1])
    parts.columns = ['REF', 'ALT']

    # Keep rows where we truly have exactly 2 alleles (ALT exists) and ALT doesn't still contain commas
    ok = parts['ALT'].notna() & (parts['ALT'].str.contains(',', na=False) == False)

    df = df.loc[ok].copy()
    df[['REF', 'ALT']] = parts.loc[ok, ['REF', 'ALT']]

    # Split locus
    df[['#CHROM', 'POS']] = df['locus'].str.split(':', n=1, expand=True)
    df['POS'] = df['POS'].astype(int)

    df.drop(columns=['alleles', 'locus'], inplace=True)

    print(f"After allele parsing: {len(df)} rows remain.")
    return df

def filter_single_letter_alleles(df: pd.DataFrame) -> pd.DataFrame:
    """Keep only rows where both REF and ALT alleles are single characters."""
    df = df[(df['REF'].str.len() == 1) & (df['ALT'].str.len() == 1)]
    print(f"Filtered to {len(df)} rows with single-letter alleles.")
    return df

def outer_merge_gnomXmiss(df: pd.DataFrame, input_file: str, hg38_dir: str) -> pd.DataFrame:
    """
    Perform an outer merge between the current DataFrame and the hg38 missense file.
    The merge is done on the '#CHROM', 'POS', 'REF', and 'ALT' columns.
    """
    # Extract uniprot_id from the input file name
    uniprot_id = os.path.basename(input_file).split('_')[0]
    
    # Construct the file path for the hg38 data file
    hg38_file = os.path.join(hg38_dir, f"{uniprot_id}_hg38.tsv")
    
    # Check if the file exists
    if not os.path.exists(hg38_file):
        raise FileNotFoundError(f"File {hg38_file} not found.")
    
    # Load the hg38 file
    hg38_df = pd.read_csv(hg38_file, sep='\t')
    print(f"Loaded hg38 data from {hg38_file}")

    # Perform outer merge on '#CHROM', 'POS', 'REF', and 'ALT'
    merged_df = pd.merge(df, hg38_df, on=['#CHROM', 'POS', 'REF', 'ALT'], how='outer', suffixes=('_gnomad', '_hg38'))

    print(f"Performed outer merge on '#CHROM', 'POS', 'REF', 'ALT'.")
    
     # Fill missing values in the 'uniprot_id' column with the uniprot_id from the hg38 data
    if 'uniprot_id_hg38' in merged_df.columns:
        merged_df['uniprot_id'] = merged_df['uniprot_id_gnomad'].combine_first(merged_df['uniprot_id_hg38'])
        merged_df.drop(columns=['uniprot_id_gnomad', 'uniprot_id_hg38'], inplace=True)
    else:
        merged_df['uniprot_id'] = merged_df['uniprot_id'].fillna(uniprot_id)

    print(f"Filled missing 'uniprot_id' values with data from hg38 or with '{uniprot_id}'.")
    
    return merged_df

##Impute
# VEP Related Functions
def run_VEP(vcf_file: str, output_vcf: str, vep_sif: str, cache_dir: str, output_dir: str, clinvar_vcf: str):
    """Runs VEP (Variant Effect Predictor) on the given VCF file with ClinVar annotation."""
    # Get SLURM job and array task ID to ensure unique HTML output for each job
    slurm_job_id = os.getenv('SLURM_JOB_ID', 'default_job')
    slurm_array_task_id = os.getenv('SLURM_ARRAY_TASK_ID', 'default_task')
    
    # HTML output path
    html_output = os.path.join(output_dir, f"output_{slurm_job_id}_{slurm_array_task_id}.vcf_summary.html")

    # Ensure clinvar_vcf directory is properly bound
    clinvar_dir = os.path.dirname(os.path.abspath(clinvar_vcf))
    
    vep_command = [
        "apptainer", "exec", 
        "--bind", f"{cache_dir}:/opt/vep/cache",
        "--bind", f"{output_dir}:{output_dir}",
        "--bind", f"{clinvar_dir}:{clinvar_dir}",  # Bind ClinVar directory
        vep_sif,
        "vep",
        "--dir_cache", "/opt/vep/cache",
        "--cache",
        "--offline",
        "--format", "vcf",
        "--vcf",
        "--force_overwrite",
        "--input_file", vcf_file,
        "--output_file", output_vcf,
        "--symbol",
        "--protein",
        "--max_af",
        "--check_existing",  # Get existing variant info
        "--custom", f"{clinvar_vcf},ClinVar,vcf,exact,0,CLNSIG,CLNDN,CLNDISDB,CLNVI",  # Only needed ClinVar fields
        "--fork", "8"
    ]

    try:
        subprocess.run(vep_command, check=True)
        print(f"VEP successfully annotated variants with ClinVar data from {clinvar_vcf}")
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"VEP failed: {str(e)}")

    # Remove all .html files from the output directory
    html_files = glob.glob(os.path.join(output_dir, "*.html"))
    for html_file in html_files:
        if os.path.exists(html_file):
            os.remove(html_file)
            print(f"Removed HTML output file: {html_file}")

def extract_VEP_info(vep_output_vcf: str) -> pd.DataFrame:
    """Extract protein_variant, enst_id, max_af, and selected ClinVar information from the VEP VCF output."""
    protein_variants = []
    enst_ids = []
    max_afs = []
    existing_variations = []
    clinvar_significances = []
    clinvar_conditions = []
    clinvar_disease_ids = []
    clinvar_variant_ids = []

    with open(vep_output_vcf, 'r') as file:
        csq_fields_map = {}
        
        for line in file:
            if line.startswith('#'):
                if "##INFO=<ID=CSQ" in line:
                    csq_format_str = line.split("Format: ")[1].strip('">').split("|")
                    for i, field in enumerate(csq_format_str):
                        csq_fields_map[field] = i
                continue

            columns = line.strip().split('\t')
            info_field = columns[7]

            protein_variant = enst_id = max_af = existing_var = None
            clin_sig = clin_condition = clin_disease_id = clin_var_id = None

            if "CSQ=" in info_field:
                csq_parts = [part for part in info_field.split(';') if part.startswith('CSQ=')]
                if csq_parts:
                    csq_values = csq_parts[0].split('=')[1].split(',')[0].split('|')
                    
                    try:
                        protein_pos_idx = csq_fields_map.get('Protein_position', -1)
                        amino_acids_idx = csq_fields_map.get('Amino_acids', -1)
                        feature_idx = csq_fields_map.get('Feature', -1)
                        max_af_idx = csq_fields_map.get('MAX_AF', -1)
                        existing_var_idx = csq_fields_map.get('Existing_variation', -1)
                        
                        clin_sig_idx = csq_fields_map.get('ClinVar_CLNSIG', -1)
                        clin_condition_idx = csq_fields_map.get('ClinVar_CLNDN', -1)
                        clin_disease_id_idx = csq_fields_map.get('ClinVar_CLNDISDB', -1)
                        clin_var_id_idx = csq_fields_map.get('ClinVar_CLNVI', -1)
                        
                        if protein_pos_idx >= 0 and amino_acids_idx >= 0:
                            protein_pos = csq_values[protein_pos_idx] if len(csq_values) > protein_pos_idx else None
                            amino_acids = csq_values[amino_acids_idx] if len(csq_values) > amino_acids_idx else None
                            
                            if protein_pos and amino_acids and '/' in amino_acids:
                                aa_ref, aa_alt = amino_acids.split('/')
                                protein_variant = f"{aa_ref}{protein_pos}{aa_alt}"
                        
                        enst_id = csq_values[feature_idx] if feature_idx >= 0 and len(csq_values) > feature_idx else None
                        
                        max_af_str = csq_values[max_af_idx] if max_af_idx >= 0 and len(csq_values) > max_af_idx else None
                        max_af = float(max_af_str) if max_af_str and max_af_str != '' else None
                        
                        existing_var = csq_values[existing_var_idx] if existing_var_idx >= 0 and len(csq_values) > existing_var_idx else None
                        
                        clin_sig = csq_values[clin_sig_idx] if clin_sig_idx >= 0 and len(csq_values) > clin_sig_idx else None
                        clin_condition = csq_values[clin_condition_idx] if clin_condition_idx >= 0 and len(csq_values) > clin_condition_idx else None
                        clin_disease_id = csq_values[clin_disease_id_idx] if clin_disease_id_idx >= 0 and len(csq_values) > clin_disease_id_idx else None
                        clin_var_id = csq_values[clin_var_id_idx] if clin_var_id_idx >= 0 and len(csq_values) > clin_var_id_idx else None
                        
                    except (IndexError, ValueError) as e:
                        print(f"Error parsing CSQ field: {e}")

            protein_variants.append(protein_variant)
            enst_ids.append(enst_id)
            max_afs.append(max_af)
            existing_variations.append(existing_var)
            clinvar_significances.append(clin_sig)
            clinvar_conditions.append(clin_condition)
            clinvar_disease_ids.append(clin_disease_id)
            clinvar_variant_ids.append(clin_var_id)

    return pd.DataFrame({
        'protein_variant_vep': protein_variants,
        'enst': enst_ids,
        'max_af': max_afs,
        'existing_variation': existing_variations,
        'clinvar_significance': clinvar_significances,
        'clinvar_condition': clinvar_conditions,
        'clinvar_disease_id': clinvar_disease_ids,
        'clinvar_variant_id': clinvar_variant_ids
    })

def add_VEP_info(df: pd.DataFrame, vep_sif: str, cache_dir: str, output_dir: str, clinvar_vcf: str) -> pd.DataFrame:
    """Adds VEP annotations including ClinVar data by running VEP on a temporary VCF."""
    unique_id = uuid.uuid4()
    temp_vcf_name = os.path.join(output_dir, f"temp_{unique_id}.vcf")
    temp_output_vcf_name = os.path.join(output_dir, f"output_{unique_id}.vcf")
    
    with open(temp_vcf_name, 'w') as temp_vcf:
        temp_vcf.write("##fileformat=VCFv4.2\n")
        temp_vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for _, row in df.iterrows():
            temp_vcf.write(f"{row['#CHROM']}\t{row['POS']}\t.\t{row['REF']}\t{row['ALT']}\t.\t.\t.\n")
    
    run_VEP(temp_vcf_name, temp_output_vcf_name, vep_sif, cache_dir, output_dir, clinvar_vcf)
    vep_info_df = extract_VEP_info(temp_output_vcf_name)
    
    df['protein_variant'] = df['protein_variant'].combine_first(vep_info_df['protein_variant_vep'])
    df['enst'] = vep_info_df['enst']
    df['max_af'] = vep_info_df['max_af']
    df['existing_variation'] = vep_info_df['existing_variation']
    df['clinvar_significance'] = vep_info_df['clinvar_significance']
    df['clinvar_condition'] = vep_info_df['clinvar_condition']
    df['clinvar_disease_id'] = vep_info_df['clinvar_disease_id']
    df['clinvar_variant_id'] = vep_info_df['clinvar_variant_id']

    if os.path.exists(temp_vcf_name):
        os.remove(temp_vcf_name)
    if os.path.exists(temp_output_vcf_name):
        os.remove(temp_output_vcf_name)

    return df

# AM imputation
def rename_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Select and reorder the necessary columns."""
    
    if 'transcript_id' in df.columns:
        df.rename(columns={'transcript_id': 'am_transcript_id_gencode_v32'}, inplace=True)
    
    if 'enst' in df.columns:
        df.rename(columns={'enst': 'enst_gencode_v46'}, inplace=True)
    return df

def impute_am_transcript_id(df: pd.DataFrame) -> pd.DataFrame:
    """
    Imputes missing values in 'am_transcript_id_gencode_v32' using the 'enst_gencode_v46'
    values from other rows where 'am_transcript_id_gencode_v32' is not missing.
    After imputation, it removes rows where 'am_transcript_id_gencode_v32' is still missing.
    """
    # Create a mapping of 'enst_gencode_v46' to 'am_transcript_id_gencode_v32' from rows where 'am_transcript_id_gencode_v32' is not empty
    enst_to_am_transcript_map = df[df['am_transcript_id_gencode_v32'].notna() & df['enst_gencode_v46'].notna()].set_index('enst_gencode_v46')['am_transcript_id_gencode_v32'].to_dict()
    
    # Apply this mapping to fill missing 'am_transcript_id_gencode_v32' values where 'enst_gencode_v46' exists
    df['am_transcript_id_gencode_v32'] = df.apply(
        lambda row: enst_to_am_transcript_map.get(row['enst_gencode_v46'], row['am_transcript_id_gencode_v32']) 
        if pd.isna(row['am_transcript_id_gencode_v32']) else row['am_transcript_id_gencode_v32'],
        axis=1
    )
    
    # Remove rows where 'am_transcript_id_gencode_v32' is still missing
    df = df[df['am_transcript_id_gencode_v32'].notna()]
    
    return df

def impute_am_pathogenicity_class_splice(df: pd.DataFrame, input_file: str, aa_sub_dir: str) -> pd.DataFrame:
    """
    Imputes missing values in 'am_pathogenicity' and 'am_class' for missing values from an external file based on 'protein_variant' matches.
    Updates 'uniprot_id' to 'splice' for rows where 'protein_variant' does not match the format requirements in the amino acid substitution data.
    
    Parameters:
    df (pd.DataFrame): Main DataFrame where imputation is to be done.
    input_file (str): Path to the input file from which 'uniprot_id' will be extracted.
    aa_sub_dir (str): Directory containing the amino acid substitution (aa_sub) files for imputation.

    Returns:
    pd.DataFrame: Updated DataFrame with imputed values and potentially modified 'uniprot_id'.
    """
    
    # Extract uniprot_id from the input file name
    uniprot_id = os.path.basename(input_file).split('_')[0]
    
    # Construct the file path for the amino acid substitution (aa_sub) file
    aa_sub_file = os.path.join(aa_sub_dir, f"{uniprot_id}_aa_substitutions.tsv")
    
    # Check if the file exists
    if not os.path.exists(aa_sub_file):
        raise FileNotFoundError(f"File {aa_sub_file} not found.")

    # Load the amino acid substitution (aa_sub) file into a DataFrame
    aa_sub_df = pd.read_csv(aa_sub_file, sep='\t')
    print(f"Loaded amino acid substitution data from {aa_sub_file}")

    # Imputation logic: for rows with missing 'am_pathogenicity' and 'am_class', match using 'protein_variant'
    for index, row in df.iterrows():
        if pd.isna(row['am_pathogenicity']) and pd.isna(row['am_class']):
            # Find the matching row in the aa_sub_df based on 'protein_variant'
            match = aa_sub_df[aa_sub_df['protein_variant'] == row['protein_variant']]
            if not match.empty:
                # Impute values if a match is found
                df.at[index, 'am_pathogenicity'] = match['am_pathogenicity'].values[0]
                df.at[index, 'am_class'] = match['am_class'].values[0]

    # Additional feature: Check if the first letter and number of protein_variant match
    for index, row in df.iterrows():
        if pd.isna(row['am_pathogenicity']) and pd.isna(row['am_class']):
            protein_variant = row['protein_variant']
            if isinstance(protein_variant, str):
                # Extract first letter and number from protein_variant
                match = False
                for aa_sub_variant in aa_sub_df['protein_variant']:
                    if isinstance(aa_sub_variant, str):
                        if (protein_variant[0] == aa_sub_variant[0] and
                                protein_variant[1:-1] == aa_sub_variant[1:-1]):
                            match = True
                            break
                # Update uniprot_id if no match found
                if not match:
                    df.at[index, 'uniprot_id'] = 'splice'

    return df

def update_am_class_based_on_protein_variant(df: pd.DataFrame) -> pd.DataFrame:
    """
    Update the 'am_pathogenicity' column based on conditions related to the 'protein_variant' column.

    Conditions:
    1. If 'protein_variant' contains a '*', set 'am_pathogenicity' to "new_stop".
    2. If the first and last characters of 'protein_variant' are the same, or if it contains '*' at both ends,
       set 'am_pathogenicity' to "synonymous".

    Parameters:
    df (pd.DataFrame): DataFrame in which to perform the operation.

    Returns:
    pd.DataFrame: Updated DataFrame with 'am_pathogenicity' values updated based on conditions.
    """

    # Condition 1: Contains '*'
    df.loc[df['protein_variant'].str.contains(r'\*', na=False), 'am_class'] = "new_stop"

    # Condition 2: First and last characters are the same or format matches "*[number]*"
    mask_synonymous = df['protein_variant'].str.match(r'^\*.*\*$', na=False) | \
                      df['protein_variant'].str.match(r'^([A-Z])[0-9]*\1$', na=False)
    df.loc[mask_synonymous, 'am_class'] = "synonymous"

    return df

def filter_invalid_entries(df: pd.DataFrame) -> pd.DataFrame:
    """
    Filters out rows with:
    1. 'uniprot_id' values being 'splice'.
    2. 'am_class' values being 'new_stop' or 'synonymous'.
    3. 'protein_variant' matching the pattern '[A-Z]1[A-Z]', or empty, NaN, or None.
    """
    df = df[df['uniprot_id'] != 'splice']
    
    if 'am_class' in df.columns:
        df = df[(df['am_class'] != 'new_stop') & (df['am_class'] != 'synonymous')]
    
    df = df[~df['protein_variant'].isna() & 
            ~(df['protein_variant'].isin([None, ''])) & 
            ~df['protein_variant'].str.contains(r'[A-Z]1[A-Z]', na=False)]
    
    return df

## Tidy
def extract_and_sort_residue_number(df: pd.DataFrame) -> pd.DataFrame:
    """
    Extracts the numeric part of 'protein_variant', writes it to a new column 'residue_number',
    and sorts the DataFrame hierarchically by 'residue_number' and then the last character of 'protein_variant'.

    Parameters:
    df (pd.DataFrame): Input DataFrame with a 'protein_variant' column.

    Returns:
    pd.DataFrame: DataFrame with a new 'residue_number' column and sorted hierarchically.
    """
    # Extract numeric part from 'protein_variant' and convert to integer
    df['residue_number'] = df['protein_variant'].str.extract(r'(\d+)').astype('Int64')
    
    # Extract the last character of the protein_variant for secondary sorting
    df['last_character'] = df['protein_variant'].str[-1]
    
    # Sort the DataFrame by 'residue_number' and then by 'last_character'
    df = df.sort_values(by=['residue_number', 'last_character'], ascending=[True, True]).reset_index(drop=True)
    
    # Drop the helper column 'last_character'
    df = df.drop(columns=['last_character'])
    
    return df

def reorder_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Reorders the columns of a DataFrame according to a predefined order,
    excluding 'intron' and 'max_af' from the final output.

    Parameters:
    df (pd.DataFrame): The input DataFrame to be reordered.

    Returns:
    pd.DataFrame: The DataFrame with columns reordered and specific columns excluded.
    """
    # Define columns that should be excluded from the final output
    EXCLUDED_COLUMNS = ['genome','intron', 'max_af','nfe_allele_count', 'nfe_allele_number', 'fin_allele_count', 'fin_allele_number',
        'mid_allele_count', 'mid_allele_number', 'remaining_allele_count', 'remaining_allele_number',
        'ami_allele_count', 'ami_allele_number', 'asj_allele_count', 'asj_allele_number', 'afr_allele_count', 'afr_allele_number',
        'eas_allele_count', 'eas_allele_number', 'sas_allele_count', 'sas_allele_number', 'amr_allele_count', 'amr_allele_number',
        'XX_allele_count', 'XX_allele_number', 'XY_allele_count', 'XY_allele_number']
    
    # Create a copy of the dataframe without the excluded columns
    output_df = df.copy()
    for col in EXCLUDED_COLUMNS:
        if col in output_df.columns:
            output_df = output_df.drop(columns=[col])
    
    # Define the order of columns for the final output
    COLUMN_ORDER = [
        'uniprot_id', '#CHROM', 'POS', 'REF', 'ALT', 
        'protein_variant', 'residue_number', 'am_class', 'am_pathogenicity',
        # Only requested ClinVar fields
        'existing_variation',
        'clinvar_significance', 'clinvar_condition', 
        'clinvar_disease_id', 'clinvar_variant_id',
        'am_transcript_id_gencode_v32', 'enst_gencode_v46',
        # Continue with other columns
        'allele_count', 'allele_number', 'allele_frequency',
        'homozygote_count', 'nfe_allele_frequency', 'fin_allele_frequency', 'mid_allele_frequency', 
        'remaining_allele_frequency', 'ami_allele_frequency', 'asj_allele_frequency', 
        'afr_allele_frequency', 'eas_allele_frequency', 'sas_allele_frequency', 
        'amr_allele_frequency', 'XX_allele_frequency', 'XY_allele_frequency'
    ]
    
    # Reorder columns based on COLUMN_ORDER, keeping extra columns at the end
    reordered_columns = [col for col in COLUMN_ORDER if col in output_df.columns]
    extra_columns = [col for col in output_df.columns if col not in reordered_columns]
    
    return output_df[reordered_columns + extra_columns]

# Workflow
def process_file(input_file: str, output_dir: str, hg38_dir: str, aa_sub_dir: str, 
                 vep_sif: str, cache_dir: str, clinvar_vcf: str, filter_single_bp: bool = False) -> None:
    """Load, transform, filter, add 'uniprot_id', run VEP with ClinVar, and save GnomAD TSV data."""
    os.makedirs(output_dir, exist_ok=True)
    # Load & Transform data
    df = load_tsv(input_file)
    df = transform_gnomad_data(df)
    if filter_single_bp:
        df = filter_single_letter_alleles(df)
    df = outer_merge_gnomXmiss(df, input_file, hg38_dir)

    # Imputation
    # Prepare 'uniprot_id'
    base_name = os.path.basename(input_file).split('_')[0]
    df['uniprot_id'] = base_name
    
    # Add VEP annotations including ClinVar data
    df = add_VEP_info(df, vep_sif, cache_dir, output_dir, clinvar_vcf)
    df = rename_columns(df)
    
    # Impute alphamissense data 
    df = impute_am_transcript_id(df)
    df = impute_am_pathogenicity_class_splice(df, input_file, aa_sub_dir)
    df = update_am_class_based_on_protein_variant(df)
    df = filter_invalid_entries(df)
    
    # Tidy
    df = extract_and_sort_residue_number(df)
    df = reorder_columns(df)  # This will now exclude 'intron' and 'max_af' columns
    
    # Save Data
    output_file_name = f"{base_name}_amissenseXgnomAD.tsv" # Prepare output path
    output_path = os.path.join(output_dir, output_file_name)
    df.to_csv(output_path, sep='\t', index=False) # Save the merged data
    print(f"Saved annotated data with VEP and ClinVar to {output_path}")

if __name__ == "__main__":
    # Get command-line arguments using sys.argv
    input_file = sys.argv[1]
    output_dir = sys.argv[2]
    vep_sif = sys.argv[3]
    cache_dir = sys.argv[4]
    hg38_dir = sys.argv[5]
    aa_sub_dir = sys.argv[6]
    clinvar_vcf = sys.argv[7]  # Parameter for ClinVar VCF file
    filter_single_bp = bool(int(sys.argv[8])) if len(sys.argv) > 8 else False

    # Call the process_file function with the passed arguments
    process_file(input_file, output_dir, hg38_dir, aa_sub_dir, vep_sif, cache_dir, clinvar_vcf, filter_single_bp)