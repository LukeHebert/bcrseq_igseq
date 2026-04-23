#!/usr/bin/env python3
"""
This script clusters B cell receptor sequences into clonal lineages using
an approach inspired by Gupta et al. 2017.
It groups sequences by identical V gene, J gene, and CDRH3 amino acid (cdr3_aa) length.
Within each group, it computes pairwise Hamming distances between CDRH3 amino acid sequences,
and then uses single–linkage hierarchical clustering to assign sequences to clonal lineages.
Importantly, the script can automatically determine the Hamming distance threshold 
via a "distance-to-nearest" analysis.

Gupta NT, Adams KD, Briggs AW, Timberlake SC, Vigneault F, Kleinstein SH. 
Hierarchical Clustering Can Identify B Cell Clones with High Confidence in Ig 
Repertoire Sequencing Data. J Immunol. 2017 Mar 15;198(6):2489-2499. 
doi: 10.4049/jimmunol.1601850. Epub 2017 Feb 8. PMID: 28179494; PMCID: 
PMC5340603.

Use `python gupta_cluster.py --help` for instructions
"""


import argparse
import os
import pandas as pd
import numpy as np
from datetime import datetime
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt

def parse_arguments():
    """Parse and return command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Cluster BCR sequences into clonal lineages based on V/J gene, '
                    'CDRH3 amino acid length, and Hamming distance.')
    parser.add_argument('annotation_tsv', type=str,
                        help='Path to the input TSV file from IgBLAST annotation.')
    parser.add_argument('--threshold', type=int, default=None,
                        help='Hamming distance threshold for clustering (number of mismatches).')
    parser.add_argument('--auto_threshold', action='store_true',
                        help='Automatically determine the Hamming distance threshold '
                             'using the distance-to-nearest method.')
    return parser.parse_args()

def create_output_dirs(tsv_path):
    """
    Create a 'clustering/' subdirectory within the input TSV's directory and
    return its path.
    """
    base_path = os.path.dirname(tsv_path)
    clustering_dir = os.path.join(base_path, 'clustering')
    os.makedirs(clustering_dir, exist_ok=True)
    return clustering_dir

def setup_logging(tsv_path, clustering_dir):
    """Create and return a log file handle and its path inside clustering_dir."""
    log_filename = datetime.now().strftime("log_clustering_%Y-%m-%d_%H-%M-%S.txt")
    log_path = os.path.join(clustering_dir, log_filename)
    return open(log_path, 'w'), log_path

def extract_gene(gene_call):
    """
    Extract the gene name from the 'v_call' or 'j_call' field.
    Assumes the gene call is a comma-delimited string and returns the first value.
    """
    if pd.isnull(gene_call):
        return None
    return gene_call.split(',')[0]

def hamming_distance(seq1, seq2):
    """
    Compute the Hamming distance between two sequences of equal length.
    """
    return sum(ch1 != ch2 for ch1, ch2 in zip(seq1, seq2))

def cluster_group(seqs, threshold):
    """
    Given a list of CDRH3 amino acid sequences (all the same length),
    compute the pairwise Hamming distance matrix, perform single-linkage clustering,
    and return an array of cluster labels.
    
    The error previously encountered (a 1D array) is fixed by converting each
    sequence into a list of characters to form a 2D array.
    """
    if len(seqs) == 0:
        return np.array([], dtype=int)
    if len(seqs) == 1:
        return np.array([1], dtype=int)
    
    # Convert each sequence string to a list of characters to obtain a 2D array.
    arr = np.array([list(seq) for seq in seqs])
    distances = pdist(arr, lambda u, v: hamming_distance(''.join(u), ''.join(v)))
    Z = linkage(distances, method='single')
    clusters = fcluster(Z, t=threshold, criterion='distance')
    return clusters

def compute_auto_threshold(df, clustering_dir, log_file):
    """
    Automatically determine the Hamming distance threshold based on the 
    "distance-to-nearest" distribution.
    
    For each group (defined by same V gene, J gene, and CDRH3 amino acid length)
    with at least two sequences, compute the minimum (nearest neighbor) Hamming distance.
    Pool these distances, estimate a kernel density (KDE) using a Gaussian kernel,
    then compute the first and second derivatives of the KDE.
    The threshold is selected as the first grid point after the main peak where the first derivative
    changes from negative to positive (indicating a local minimum) and the second derivative is positive.
    A PNG visualization (dpi=600) is saved to the 'clustering/' subdirectory.
    """
    distances = []
    groups = df.groupby(['v_gene', 'j_gene', 'cdr3_aa_length'])
    for name, group in groups:
        if len(group) < 2:
            continue
        seqs = group['cdr3_aa'].tolist()
        n = len(seqs)
        for i in range(n):
            dists = [hamming_distance(seqs[i], seqs[j]) for j in range(n) if j != i]
            if dists:
                distances.append(min(dists))
    distances = np.array(distances)
    if len(distances) == 0:
        log_file.write("No groups with >=2 sequences found; using default threshold of 2.\n")
        return 2

    # Estimate the density using Gaussian KDE.
    kde = gaussian_kde(distances)
    grid = np.linspace(distances.min(), distances.max(), 1000)
    kde_vals = kde(grid)
    
    # Compute first and second derivatives of the KDE.
    d1 = np.gradient(kde_vals, grid)
    d2 = np.gradient(d1, grid)
    
    # Identify the index of the main peak (the maximum of the KDE).
    main_peak_idx = np.argmax(kde_vals)
    
    # Find the first index after the main peak where the first derivative changes from negative to positive 
    # and the second derivative is positive.
    candidate_idx = None
    for i in range(main_peak_idx + 1, len(grid)):
        if d1[i-1] < 0 and d1[i] >= 0 and d2[i] > 0:
            candidate_idx = i
            break
    if candidate_idx is None:
        # If no candidate is found, default to the median.
        auto_thresh = int(np.round(np.median(distances)))
        log_file.write(f"No clear valley found; using median threshold {auto_thresh}.\n")
    else:
        auto_thresh = grid[candidate_idx]
        auto_thresh = int(np.round(auto_thresh))
        log_file.write(f"Automatically determined threshold based on derivatives: {auto_thresh} mismatches.\n")
    
    # Plot KDE along with its first derivative zero crossing and the selected threshold.
    plt.figure(figsize=(8, 6))
    plt.plot(grid, kde_vals, label='KDE of distances')
    plt.axvline(x=grid[main_peak_idx], color='green', linestyle=':', label='Main peak')
    plt.axvline(x=auto_thresh, color='red', linestyle='--', label=f"Threshold = {auto_thresh}")
    plt.xlabel("Nearest neighbor Hamming distance")
    plt.ylabel("Density")
    plt.title("Distance-to-Nearest Distribution with Threshold")
    plt.legend()
    png_filename = os.path.join(clustering_dir, "distance_to_nearest.png")
    plt.savefig(png_filename, dpi=600, bbox_inches='tight')
    plt.close()
    log_file.write(f"KDE plot saved to {png_filename}\n")
    return auto_thresh

def perform_clustering(df, threshold, log_file):
    """
    Group the DataFrame by V gene, J gene, and CDRH3 amino acid length,
    perform hierarchical clustering within each group using the specified threshold,
    and return the DataFrame with an added 'ClusterID' column.
    
    It is assumed that the DataFrame already contains 'v_gene', 'j_gene',
    and 'cdr3_aa_length' columns.
    """
    groups = df.groupby(['v_gene', 'j_gene', 'cdr3_aa_length'])
    log_file.write(f"Number of groups to cluster: {len(groups)}\n")
    
    df['ClusterID'] = np.nan
    global_cluster_id = 1
    
    for name, group in groups:
        indices = group.index
        seqs = group['cdr3_aa'].tolist()
        clusters = cluster_group(seqs, threshold)
        clusters_global = clusters + global_cluster_id - 1
        df.loc[indices, 'ClusterID'] = clusters_global
        global_cluster_id += clusters.max()
    df['ClusterID'] = df['ClusterID'].astype(int)
    return df

def main():
    """Main function that orchestrates the clustering process."""
    args = parse_arguments()
    clustering_dir = create_output_dirs(args.annotation_tsv)
    log_file, log_path = setup_logging(args.annotation_tsv, clustering_dir)
    try:
        log_file.write(f"Processing started for {args.annotation_tsv}\n")
        
        # Load the input TSV file.
        df = pd.read_csv(args.annotation_tsv, sep='\t')
        log_file.write(f"Initial row count: {len(df)}\n")
        
        # Preprocess: add v_gene, j_gene, and cdr3_aa_length columns.
        df['v_gene'] = df['v_call'].apply(extract_gene)
        df['j_gene'] = df['j_call'].apply(extract_gene)
        df['cdr3_aa_length'] = df['cdr3_aa'].apply(lambda x: len(x) if pd.notnull(x) else 0)
        
        # Determine threshold.
        if args.auto_threshold:
            threshold = compute_auto_threshold(df, clustering_dir, log_file)
        elif args.threshold is not None:
            threshold = args.threshold
            log_file.write(f"Using user-specified threshold: {threshold}\n")
        else:
            threshold = 2
            log_file.write("No threshold specified; defaulting to 2 mismatches.\n")
        
        # Perform clustering.
        clustered_df = perform_clustering(df, threshold, log_file)
        clustered_df = clustered_df.sort_values(by=['ClusterID'], ascending=True)
        
        # Save the final clustered TSV file in the same directory as the input TSV.
        output_file = os.path.splitext(args.annotation_tsv)[0] + '_clustered.tsv'
        clustered_df.to_csv(output_file, sep='\t', index=False)
        log_file.write(f"Final sorted dataframe saved to {output_file}\n")
        log_file.write("Processing completed successfully.\n")
    finally:
        log_file.close()

if __name__ == "__main__":
    main()