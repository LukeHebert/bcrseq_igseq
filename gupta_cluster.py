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
from scipy.signal import find_peaks
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
    """
    if len(seqs) == 0:
        return np.array([], dtype=int)
    if len(seqs) == 1:
        return np.array([1], dtype=int)
    
    # Convert each sequence (string) to a list of its characters.
    # This results in a 2D numpy array with shape (number_of_sequences, sequence_length)
    arr = np.array([list(seq) for seq in seqs])
    
    # Use pdist with a custom lambda function to compute the pairwise Hamming distance.
    # The lambda function receives two 1D arrays (each representing a sequence as a list of characters)
    # and computes the Hamming distance after joining them back to strings.
    distances = pdist(arr, lambda u, v: hamming_distance(''.join(u), ''.join(v)))
    
    # Perform single-linkage hierarchical clustering on the distance matrix.
    Z = linkage(distances, method='single')
    
    # Form flat clusters with the given threshold (max number of mismatches allowed).
    clusters = fcluster(Z, t=threshold, criterion='distance')
    return clusters

def compute_auto_threshold(df, clustering_dir, log_file):
    """
    Automatically determine the Hamming distance threshold based on the 
    "distance-to-nearest" distribution.
    
    For each group (defined by same V gene, J gene, and CDRH3 amino acid length)
    with at least two sequences, compute the minimum (nearest neighbor) Hamming distance.
    Pool these distances, estimate a kernel density (KDE), and choose the threshold
    as the valley between the first two peaks. Also save a PNG visualization
    of the KDE with detected peaks and the selected threshold.
    """
    distances = []
    groups = df.groupby(['v_gene', 'j_gene', 'cdr3_aa_length'])
    for name, group in groups:
        if len(group) < 2:
            continue
        seqs = group['cdr3_aa'].tolist()
        n = len(seqs)
        for i in range(n):
            # Compute Hamming distances to all other sequences in the group.
            dists = [hamming_distance(seqs[i], seqs[j]) for j in range(n) if j != i]
            if dists:
                distances.append(min(dists))
    distances = np.array(distances)
    if len(distances) == 0:
        log_file.write("No groups with >=2 sequences found; using default threshold of 2.\n")
        return 2

    # Compute a kernel density estimate (KDE) over the distances.
    kde = gaussian_kde(distances)
    grid = np.linspace(distances.min(), distances.max(), 1000)
    kde_vals = kde(grid)
    
    # Identify peaks in the KDE.
    peaks, _ = find_peaks(kde_vals)
    if len(peaks) < 2:
        auto_thresh = int(np.round(np.median(distances)))
        log_file.write(f"Less than two peaks found; using median threshold {auto_thresh}.\n")
    else:
        first_peak = peaks[0]
        second_peak = peaks[1]
        # The valley is the grid point with the minimum density between the first two peaks.
        # These peaks should belong to fellow B cell lineage members & non-members
        valley_index = np.argmin(kde_vals[first_peak:second_peak]) + first_peak
        auto_thresh = grid[valley_index]
        auto_thresh = int(np.round(auto_thresh))
        log_file.write(f"Automatically determined threshold based on distance-to-nearest: {auto_thresh} mismatches.\n")
    
    # Plot the KDE and mark the detected peaks and threshold.
    plt.figure(figsize=(8, 6))
    plt.plot(grid, kde_vals, label='KDE of distances')
    plt.xlabel("Nearest neighbor Hamming distance")
    plt.ylabel("Density")
    plt.title("Distance-to-Nearest Distribution")
    if len(peaks) > 0:
        plt.plot(grid[peaks], kde_vals[peaks], "x", label="Peaks")
    if len(peaks) >= 2:
        plt.axvline(x=auto_thresh, color='red', linestyle='--', label=f"Threshold = {auto_thresh}")
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
    """
    groups = df.groupby(['v_gene', 'j_gene', 'cdr3_aa_length'])
    log_file.write(f"Number of groups to cluster: {len(groups)}\n")
    
    df['ClusterID'] = np.nan
    global_cluster_id = 1
    
    for name, group in groups:
        indices = group.index
        seqs = group['cdr3_aa'].tolist()
        clusters = cluster_group(seqs, threshold)
        # Assign global cluster IDs by offsetting with global_cluster_id.
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
        
        # Determine threshold: either automatically or use user-specified value.
        if args.auto_threshold:
            threshold = compute_auto_threshold(df, clustering_dir, log_file)
        elif args.threshold is not None:
            threshold = args.threshold
            log_file.write(f"Using user-specified threshold: {threshold}\n")
        else:
            threshold = 2
            log_file.write("No threshold specified; defaulting to 2 mismatches.\n")
        
        # Perform clustering on the DataFrame.
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