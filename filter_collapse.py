#!/usr/bin/env python3
"""
This script processes AIRR format BCRseq TSV data by filtering out non-functional rows and collapsing identical sequences. 
It determines whether the data is heavy or light chain and applies chain-specific filtering as needed. 
The output is a filtered TSV file and a log file documenting the processing steps.

Use `filter_collapse.py --help` for instructions to run the script.
"""

import argparse
import os
import pandas as pd
from datetime import datetime

def parse_arguments():
    """Parse command-line arguments for the input TSV file and minimum count threshold."""
    parser = argparse.ArgumentParser(
        description=("Takes AIRR format BCRseq data from either a heavy or light chain "
                     "amplicon, filters out rows presumed to be non-functional, collapses "
                     "identical sequences and counts their frequency. For heavy chain data, "
                     "rows lacking cdr3_aa or containing a stop codon in cdr3_aa are removed.")
    )
    parser.add_argument('tsv_path', type=str, help='Path to the input AIRR formatted TSV file')
    mincount_default = 2
    parser.add_argument('--threshold', type=int, default=mincount_default,
                        help=f'Minimum count threshold for nucleotide sequences (default: {mincount_default})')
    return parser.parse_args()

def setup_logging(tsv_path):
    """Create and return a log file handle and its path for the given TSV file."""
    base_path = os.path.dirname(os.path.abspath(tsv_path))
    # Create a new directory 'filter_collapse' within the input file's directory if it doesn't exist
    log_dir = os.path.join(base_path, "filter_collapse")
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    log_filename = datetime.now().strftime("log_filtering_%Y-%m-%d_%H-%M-%S.txt")
    log_path = os.path.join(log_dir, log_filename)
    log_handle = open(log_path, 'w')
    return log_handle, log_path

def determine_chain_type(df, log_file):
    """Determine whether the data is heavy or light chain by tallying gene call columns and logging the results."""
    heavy_tally = 0
    light_tally = 0
    gene_cols = ['v_call', 'd_call', 'j_call', 'c_call']
    
    for col in gene_cols:
        calls = df[col].dropna().apply(lambda x: x.split(',')[0])
        tally_heavy = calls.apply(lambda s: 1 if len(s) >= 3 and s[2] == 'H' else 0).sum()
        tally_light = calls.apply(lambda s: 1 if len(s) >= 3 and s[2] in ['K', 'L'] else 0).sum()
        log_file.write(f"{col} tally - Heavy: {tally_heavy}, Light: {tally_light}\n")
        heavy_tally += tally_heavy
        light_tally += tally_light

    log_file.write(f"Total tally - Heavy: {heavy_tally}, Light: {light_tally}\n")
    return heavy_tally > light_tally

def filter(df, log_file, is_heavy_chain, threshold):
    """Apply filtering steps to the DataFrame and log each processing step."""
    df = df[df['sequence'].notna()].copy()
    log_file.write(f"After removing rows with NaN 'sequence': {len(df)} rows left\n")
    df = df[df['stop_codon'] != 'T'].copy()
    log_file.write(f"After removing rows where 'stop_codon' is 'T': {len(df)} rows left\n")
    if is_heavy_chain:
        log_file.write("Data determined to be heavy chain. Applying heavy chain filters.\n")
        df = df[df['cdr3_aa'].notna()].copy()
        log_file.write(f"After removing rows with NaN 'cdr3_aa': {len(df)} rows left\n")
        df = df[df['cdr3_aa'].apply(lambda x: '*' not in x)].copy()
        log_file.write(f"After removing rows with asterisk in 'cdr3_aa': {len(df)} rows left\n")
    else:
        log_file.write("Data determined to be light chain. Skipping heavy chain specific filters.\n")
    seq_counts = df['sequence'].value_counts()
    seq_counts_df = seq_counts.reset_index()
    seq_counts_df.columns = ['sequence', 'nt_seq_count']
    df = df.merge(seq_counts_df, on='sequence', how='left')
    df.drop_duplicates(subset=['sequence'], inplace=True)
    log_file.write(f"Row count after removing duplicate sequence values: {len(df)}\n")
    df = df[df['nt_seq_count'] >= threshold]
    log_file.write(f"Row count after filtering sequences with count < {threshold}: {len(df)}\n")
    return df

def main():
    """Process the TSV file, filter its contents, and output a new filtered TSV file with logging."""
    args = parse_arguments()
    tsv_path = args.tsv_path
    threshold = args.threshold

    log_file, log_path = setup_logging(tsv_path)
    log_file.write(f"Processing TSV: {tsv_path}\n")
    log_file.write(f"Minimum sequence count threshold: {threshold}\n")

    try:
        df = pd.read_csv(tsv_path, sep='\t', dtype=str)
    except Exception as e:
        log_file.write(f"Error reading TSV file: {e}\n")
        log_file.close()
        raise

    initial_count = len(df)
    log_file.write(f"Initial number of rows: {initial_count}\n")

    is_heavy_chain = determine_chain_type(df, log_file)
    df = filter(df, log_file, is_heavy_chain, threshold)

    base, ext = os.path.splitext(tsv_path)
    output_path = base + "_filtered.tsv"
    try:
        df.to_csv(output_path, sep='\t', index=False)
        log_file.write(f"Filtered TSV written to: {output_path}\n")
    except Exception as e:
        log_file.write(f"Error writing output TSV: {e}\n")

    log_file.write("Processing complete.\n")
    log_file.close()

if __name__ == "__main__":
    main()
