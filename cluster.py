'''
This script takes an IgBLAST output TSV file as input.
Primarily this script clusters CDRH3 sequences into estimated clonal lineages.
The output from this script will be a TSV file similar to the input file, but 
now including clustering data.

USAGE:
python cluster.py --help
'''

import argparse
import os
import pandas as pd
from datetime import datetime

def process_tsv(tsv_path, executable_path, threshold):
    # Setup log file
    base_path = os.path.dirname(tsv_path)
    log_filename = datetime.now().strftime("log_clustering_%Y-%m-%d_%H-%M-%S.txt")
    log_path = os.path.join(base_path, log_filename)
    
    with open(log_path, 'w') as log_file:
        log_file.write(f"Processing started for {tsv_path} using executable {executable_path}\n")

        # Load TSV file
        df = pd.read_csv(tsv_path, sep='\t')
        initial_row_count = len(df)
        log_file.write(f"Initial row count: {initial_row_count}\n")

        # Filter out rows where 'sequence' is NaN
        filtered_df = df[df['sequence'].notna()].copy()
        log_file.write(f"After removing rows with NaN 'sequence': {len(filtered_df)} reads left\n")

        # Filter out rows where 'cdr3_aa' is NaN
        filtered_df = filtered_df[filtered_df['cdr3_aa'].notna()].copy()
        log_file.write(f"After removing rows with NaN 'cdr3_aa': {len(filtered_df)} reads left\n")

        # Filter out rows where 'stop_codon' is 'T'
        filtered_df = filtered_df[filtered_df['stop_codon'] != 'T'].copy()
        log_file.write(f"After removing rows where 'stop_codon' is 'T': {len(filtered_df)} reads left\n")

        # Filter out rows where 'cdr3_aa' contains an asterisk
        filtered_df = filtered_df[filtered_df['cdr3_aa'].apply(lambda x: '*' not in x)].copy()
        log_file.write(f"After removing rows with asterisk in 'cdr3_aa': {len(filtered_df)} reads left\n")

        # Add nucleotide sequence counts and remove duplicates
        sequence_counts = filtered_df['sequence'].value_counts()
        sequence_counts_df = sequence_counts.reset_index()
        sequence_counts_df.columns = ['sequence', 'nt_seq_count']
        filtered_df = filtered_df.merge(sequence_counts_df, on='sequence', how='left')
        filtered_df.drop_duplicates(subset=['sequence'], inplace=True)
        after_dup_removal_row_count = len(filtered_df)
        log_file.write(f"Row count after removing duplicates: {after_dup_removal_row_count}\n")

        # Filter out rows where nt_seq_count is below the user-specified threshold
        filtered_df = filtered_df[filtered_df['nt_seq_count'] >= threshold]
        after_mincount_row_count = len(filtered_df)
        log_file.write(f"Row count after removing seqs with count < {threshold}: {after_mincount_row_count}\n")

        # Extract unique CDR3 sequences
        unique_cdr3 = filtered_df['cdr3_aa'].drop_duplicates().dropna()
        aa_seqs_file = os.path.splitext(tsv_path)[0] + '_unique_AAs.txt'
        unique_cdr3.to_csv(aa_seqs_file, index=False, header=False)

        # Run the clustering executable
        clustered_file = os.path.splitext(tsv_path)[0] + '_clustered.tsv'
        os.system(f"{executable_path} -f {aa_seqs_file} > {clustered_file}")

        # Read and merge the cluster results
        cluster_results = pd.read_csv(clustered_file, sep='\t', header=None, names=['cdr3_aa', 'ClusterID'])
        merged_df = pd.merge(filtered_df, cluster_results, on='cdr3_aa', how='left')

        # Sort dataframe
        sorted_df = merged_df.sort_values(by=['ClusterID', 'nt_seq_count', 'cdr1_aa', 'cdr2_aa', 'cdr3_aa'], ascending=[True, False, True, True, True])

        # Save the final dataframe
        sorted_df.to_csv(clustered_file, sep='\t', index=False)
        log_file.write(f"Final sorted dataframe saved to {clustered_file}\n")
        log_file.write("Processing completed successfully.\n")

def main():
    parser = argparse.ArgumentParser(description='Process a TSV file containing cDNA sequences.')
    parser.add_argument('annotation_tsv', type=str, help='Path to the input TSV file from BCR cDNA read annotation by IgBLAST.')
    parser.add_argument('jeffclust_executable', type=str, help='Path to the executable file for clustering called CDR3_Clonotyping.')
    parser.add_argument('threshold', type=int, help='Minimum count of full cDNA nucleotide sequences required to keep a row')

    args = parser.parse_args()

    process_tsv(args.annotation_tsv, args.jeffclust_executable, args.threshold)

if __name__ == "__main__":
    main()
