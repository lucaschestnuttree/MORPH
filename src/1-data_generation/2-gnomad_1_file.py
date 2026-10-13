import os
import sys
import glob
import time
import hail as hl
import pandas as pd
from typing import Dict, List, Tuple

def initialize_hail(log_dir: str, log_level: str = "WARN", quiet: bool = False) -> None:
    """Initialize the Hail environment."""
    start_time = time.time()
    
    # Set directories
    os.makedirs(log_dir, exist_ok=True)  # Ensure log directory exists
    os.environ["SPARK_LOCAL_DIRS"] = log_dir
    os.environ["HAIL_LOG_DIR"] = log_dir
    
    # Adjust Spark resource settings
    os.environ["PYSPARK_SUBMIT_ARGS"] = (
        '--driver-memory 4G --executor-memory 8G --executor-cores 4 '
        '--num-executors 2 --driver-java-options "-Dhail.log.dir={log_dir}" pyspark-shell'
    )
    
    # Initialize Hail with logging
    hl.init(log=os.path.join(log_dir, 'hail.log'), tmp_dir=log_dir, quiet=quiet)
    hl.spark_context().setLogLevel(log_level)
    
    end_time = time.time()
    print(f"Hail initialization took {end_time - start_time:.2f} seconds.")

def get_intervals_from_dataframe(df: pd.DataFrame) -> List[hl.Interval]:
    """Create a list of genomic intervals from a DataFrame with adjusted min and max positions."""
    chromosome_ranges = df.groupby('#CHROM')['POS'].agg(['min', 'max']).to_dict('index')
    
    intervals = [
        hl.Interval(hl.Locus(chrom, row['min']-100), hl.Locus(chrom, row['max']+100))
        for chrom, row in chromosome_ranges.items()
    ]
    
    return intervals

def read_and_filter_hail_table(ht_path: str, intervals: List[hl.Interval], n_partitions: int = 200) -> hl.Table:
    """Read and filter a Hail table using specified intervals, then repartition."""
    
    # First, read the table using the intervals only (no partitioning yet)
    if intervals:
        ht = hl.read_table(ht_path, _intervals=intervals)
    else:
        ht = hl.read_table(ht_path)
    
    # Then, repartition the table if n_partitions is specified
    if n_partitions:
        ht = ht.repartition(n_partitions)

    return ht

def create_frequency_table(ht: hl.Table, output_path: str, meta_keys_all: list):
    """Create and save the frequency table with subpopulation frequencies using Hail's efficient export."""
    start_time = time.time()

    # Collect the frequency metadata
    freq_meta_joint = ht.joint_globals.freq_meta.collect()[0]
    indices = {str(k): freq_meta_joint.index(k) for k in meta_keys_all}
    index = indices[str({'group': 'adj'})]

    # Annotate the Table with frequency information for the 'adj' group
    ht = ht.annotate(
        allele_count=ht.joint.freq[index].AC,
        allele_number=ht.joint.freq[index].AN,
        allele_frequency=ht.joint.freq[index].AF,
        homozygote_count=ht.joint.freq[index].homozygote_count
    )

    # Add annotations for each subpopulation frequency based on meta_keys_all
    for key in meta_keys_all[1:]:
        key_str = str(key)
        idx = indices[key_str]
        subpop_name = key.get('gen_anc', key.get('sex', ''))
        
        # Annotate the Hail Table with subpopulation frequencies
        ht = ht.annotate(**{
            f"{subpop_name}_allele_frequency": ht.joint.freq[idx].AF,
        })

    # Add annotations for each subpopulation frequency based on meta_keys_all
    for key in meta_keys_all[1:]:
        key_str = str(key)
        idx = indices[key_str]
        subpop_name = key.get('gen_anc', key.get('sex', ''))
        
        # Annotate the Hail Table with subpopulation frequencies
        ht = ht.annotate(**{
            f"{subpop_name}_allele_count": ht.joint.freq[idx].AC,
            f"{subpop_name}_allele_number": ht.joint.freq[idx].AN
        })

    # Drop unnecessary columns
    ht = ht.drop('region_flags', 'exomes', 'genomes', 'joint', 'freq_comparison_stats')

    # Export the Hail table directly to CSV
    ht.export(output_path)

    end_time = time.time()
    print(f"Creating and exporting frequency table took {end_time - start_time:.2f} seconds.")

def remove_crc_files(directory: str, base_name: str) -> None:
    """
    Remove .crc files matching the pattern <base_name>_gnomAD.tsv.crc in the given directory.

    :param directory: The directory where the .crc files are located.
    :param base_name: The dynamic base name to match (e.g., Q5VST9). Only files following the 
                      pattern <base_name>_gnomAD.tsv.crc will be removed.
    """
    # Construct the file pattern based on the base name
    pattern = os.path.join(directory, f".{base_name}_gnomAD.tsv.crc")
    
    # Find files matching the pattern
    crc_files = glob.glob(pattern)
    
    # Remove each matching .crc file
    for crc_file in crc_files:
        try:
            os.remove(crc_file)
            print(f"Removed: {crc_file}")
        except OSError as e:
            print(f"Error removing {crc_file}: {e}")

def process_file(tsv_path: str, ht_path: str, output_dir: str, log_dir: str, meta_keys_all: list,
                 read_header: bool = False, n_lines: int = 10, n_partitions: int = 200):
    """Process a single file with optional filtering by chromosome intervals."""
    total_start_time = time.time()
    
    # Extract the base name from the TSV path (e.g., P01834 from P01834_missense.tsv)
    base_name = os.path.basename(tsv_path).split('_')[0]
    
    # Create the output file name in the format P01834_gnomAD.tsv
    output_file_name = f"{base_name}_gnomAD.tsv"
    
    # Combine with output directory
    output_path = os.path.join(output_dir, output_file_name)
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Initialize Hail
    initialize_hail(log_dir)

    # Read the TSV file into a DataFrame
    if read_header:
        df = pd.read_csv(tsv_path, sep='\t', nrows=n_lines)
    else:
        df = pd.read_csv(tsv_path, sep='\t')
    
    # Get intervals from the chromosome ranges in the TSV file
    intervals = get_intervals_from_dataframe(df)

    # Read and filter the Hail table based on intervals, then repartition
    ht = read_and_filter_hail_table(ht_path, intervals, n_partitions)

    # Create and save the frequency table
    create_frequency_table(ht, output_path, meta_keys_all)

    remove_crc_files(output_dir,base_name)

    total_end_time = time.time()
    print(f"Total processing time: {total_end_time - total_start_time:.2f} seconds.")

if __name__ == "__main__":
    tsv_path = sys.argv[1]  # Get TSV path from command-line argument
    ht_path = sys.argv[2]
    output_dir = sys.argv[3]
    log_dir = sys.argv[4]

    # Define meta keys for filtering frequency data in the joint dataset
    meta_keys_all = [
        {'group': 'adj'},
        {'gen_anc': 'nfe', 'group': 'adj'},
        {'gen_anc': 'fin', 'group': 'adj'},
        {'gen_anc': 'mid', 'group': 'adj'},
        {'gen_anc': 'remaining', 'group': 'adj'},
        {'gen_anc': 'ami', 'group': 'adj'},
        {'gen_anc': 'asj', 'group': 'adj'},
        {'gen_anc': 'afr', 'group': 'adj'},
        {'gen_anc': 'eas', 'group': 'adj'},
        {'gen_anc': 'sas', 'group': 'adj'},
        {'gen_anc': 'amr', 'group': 'adj'},
        {'group': 'adj', 'sex': 'XX'},
        {'group': 'adj', 'sex': 'XY'}
    ]

    # Optional parameters
    read_header = False  # Set to True to read the first n lines, False to read the entire file
    n_lines = 100        # Number of lines to read if use_n_lines is True

    # Process the file with the specified options
    process_file(tsv_path, ht_path, output_dir, log_dir, meta_keys_all, read_header, n_lines)
