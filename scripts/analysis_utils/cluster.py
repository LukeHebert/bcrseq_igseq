'''
This script takes largely AIRR-formatted BCRseq heavy chain data that has been 
filtered by `filter.py`.
Primarily this script clusters CDRH3 sequences into estimated clonal lineages.
The output from this script will be a TSV file similar to the input file, but 
now including cluster (estimated clonal lineage) assignments.

USAGE:
python cluster.py --help
'''

import argparse
import os
import pandas as pd
from datetime import datetime

def parse_arguments():
    """Parse and return command-line arguments."""
    parser = argparse.ArgumentParser(description='Process a TSV file containing cDNA sequences.')
    parser.add_argument('annotation_tsv', type=str,
                        help='Path to the input TSV file from BCR cDNA read annotation by IgBLAST.')
    parser.add_argument('clustering_executable', type=str,
                        help='Path to the executable file for clustering. Clustering executable must take a TXT file with single unlabeled column of CDR3 AA seqs & create a TSV file of two unlabeled columns, the input CDR3 AA seqs & their new cluster IDs.')
    return parser.parse_args()

def setup_logging(tsv_path):
    """Create and return a log file handle for the TSV processing along with its path."""
    base_path = os.path.dirname(tsv_path)
    log_filename = datetime.now().strftime("log_clustering_%Y-%m-%d_%H-%M-%S.txt")
    log_path = os.path.join(base_path, log_filename)
    return open(log_path, 'w'), log_path

def cluster(df, tsv_path, executable_path):
    """Cluster the CDR3 sequences and return the sorted merged DataFrame."""
    # Extract unique CDR3 sequences and write to file
    unique_cdr3 = df['cdr3_aa'].drop_duplicates().dropna()
    cdr3s_file = os.path.splitext(tsv_path)[0] + '_cdr3s_aa.txt'
    unique_cdr3.to_csv(cdr3s_file, index=False, header=False)

    # Run the clustering executable
    clusterIDs_file = os.path.splitext(tsv_path)[0] + '_clusterIDs.tsv'
    os.system(f'"{executable_path}" -f "{cdr3s_file}" > "{clusterIDs_file}"')

    # Merge the cluster assignments with the original DataFrame
    cluster_results = pd.read_csv(clusterIDs_file, sep='\t', header=None, names=['cdr3_aa', 'ClusterID'])
    merged_df = pd.merge(df, cluster_results, on='cdr3_aa', how='left')

    # Sort the DataFrame
    clustered_df = merged_df.sort_values(
        by=['ClusterID', 'nt_seq_count', 'cdr1_aa', 'cdr2_aa', 'cdr3_aa'],
        ascending=[True, False, True, True, True]
    )
    return clustered_df

def main():
    """Main function that orchestrates the clustering process."""
    args = parse_arguments()
    log_file, log_path = setup_logging(args.annotation_tsv)
    try:
        log_file.write(f"Processing started for {args.annotation_tsv} using executable {args.clustering_executable}\n")

        # Load TSV file directly in main for simplicity
        df = pd.read_csv(args.annotation_tsv, sep='\t')
        log_file.write(f"Initial row count: {len(df)}\n")

        # Cluster the data (includes writing unique CDR3s, running clustering, merging results, and sorting)
        clustered_df = cluster(df, args.annotation_tsv, args.clustering_executable)

        # Save the final sorted DataFrame
        clustered_file = os.path.splitext(args.annotation_tsv)[0] + '_clustered.tsv'
        clustered_df.to_csv(clustered_file, sep='\t', index=False)

        log_file.write(f"Final sorted dataframe saved to {clustered_file}\n")
        log_file.write("Processing completed successfully.\n")
    finally:
        log_file.close()

if __name__ == "__main__":
    main()