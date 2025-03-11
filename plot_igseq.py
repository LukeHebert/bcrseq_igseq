#!/usr/bin/env python3
"""
This script processes an antibody mass spectroscopy TSV file merged with BCR-seq
transcript data and produces bar plots showing the relative abundance of B cell
lineages (by ClusterID for heavy chain) and gene usage (v_call and j_call) for 
heavy, kappa, and light chain subsets.

A companion TSV is generated for each plot, containing the actual X tick labels
(one column) and the corresponding Y axis values (another column).
 
Use `python plot_igseq.py --help` for instructions.
"""

import os
import sys
import argparse
import logging
from datetime import datetime

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def create_output_dir(input_file):
    """Create a new 'plot_igseq' directory in the parent directory of the input file."""
    parent_dir = os.path.dirname(os.path.abspath(input_file))
    output_dir = os.path.join(parent_dir, "plot_igseq")
    os.makedirs(output_dir, exist_ok=True)
    return output_dir

def setup_logging(output_dir, cmd):
    """Setup logging with a timestamped log file in the output directory."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_filename = os.path.join(output_dir, f"log_{timestamp}.txt")
    logging.basicConfig(
        filename=log_filename,
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    logging.info("Command executed: " + " ".join(cmd))
    logging.info("Output files will be saved in: " + output_dir)
    return log_filename

def load_data(tsv_file):
    """Load the TSV file into a pandas DataFrame."""
    df = pd.read_csv(tsv_file, sep="\t")
    logging.info(f"Loaded {len(df)} rows from {tsv_file}")
    return df

def split_chains(df):
    """
    Split the data into heavy, kappa, lambda, and combined light chain subsets
    by searching the v_call column for IGHV, IGKV, and IGLV substrings.
    """
    heavy_df = df[df['v_call'].str.contains("IGHV", na=False)]
    kappa_df = df[df['v_call'].str.contains("IGKV", na=False)]
    lambda_df = df[df['v_call'].str.contains("IGLV", na=False)]
    light_df = pd.concat([kappa_df, lambda_df])
    logging.info(
        f"Heavy chain rows: {len(heavy_df)}; "
        f"Kappa chain rows: {len(kappa_df)}; "
        f"Lambda chain rows: {len(lambda_df)}; "
        f"Combined light chain rows: {len(light_df)}"
    )
    return heavy_df, kappa_df, lambda_df, light_df

def plot_heavy_cluster(heavy_df, output_dir):
    """Generate the heavy chain ClusterID bar plot and companion files."""
    # Filter heavy chain rows with cdr3_coverage_percentage > 0
    heavy_filtered = heavy_df[heavy_df['cdr3_coverage_percentage'] > 0]
    if heavy_filtered.empty:
        logging.info("No heavy chain rows with cdr3_coverage_percentage > 0; skipping heavy cluster plot.")
        return

    # Group by ClusterID and sum Precursor Abundance
    cluster_group = heavy_filtered.groupby("ClusterID")["Precursor Abundance"].sum().reset_index()
    total_cluster_abundance = cluster_group["Precursor Abundance"].sum()
    # Calculate percent abundance for each ClusterID
    cluster_group["percent_abundance"] = (
        cluster_group["Precursor Abundance"] / total_cluster_abundance
    ) * 100

    # Order the data by percent_abundance descending
    cluster_group = cluster_group.sort_values(by="percent_abundance", ascending=False)

    # Create bar plot using seaborn with ordered x-axis and vertical labels
    plt.figure(figsize=(10, 6))
    sns.barplot(
        x="ClusterID",
        y="percent_abundance",
        data=cluster_group,
        order=cluster_group["ClusterID"].tolist()
    )
    plt.ylabel("Percent Abundance")
    plt.title("Relative Abundance of B cell lineages (Heavy chain by ClusterID)")
    plt.xticks(rotation=90)
    heavy_cluster_plot_file = os.path.join(output_dir, "heavy_cluster_plot.png")
    plt.tight_layout()
    plt.savefig(heavy_cluster_plot_file, dpi=600)
    plt.close()
    logging.info("Saved heavy chain cluster plot to " + heavy_cluster_plot_file)

    # Save the aggregated data (x-axis & y-axis) as a companion TSV
    heavy_cluster_plot_data = cluster_group[["ClusterID", "percent_abundance"]]
    heavy_cluster_plot_data_file = os.path.join(output_dir, "heavy_cluster_plot_data.tsv")
    heavy_cluster_plot_data.to_csv(heavy_cluster_plot_data_file, sep="\t", index=False)
    logging.info("Saved heavy chain cluster plot data to " + heavy_cluster_plot_data_file)

    # Save companion TSV file with the filtered heavy chain data used for grouping
    heavy_cluster_tsv = os.path.join(output_dir, "heavy_cluster_data.tsv")
    heavy_filtered.to_csv(heavy_cluster_tsv, sep="\t", index=False)
    logging.info("Saved heavy chain filtered data TSV to " + heavy_cluster_tsv)

    # Write companion TXT file with stats
    total_heavy_abundance = heavy_df["Precursor Abundance"].sum()
    abundance_ratio = (
        (total_cluster_abundance / total_heavy_abundance) * 100
        if total_heavy_abundance > 0
        else 0
    )

    n_total = len(heavy_df)
    n_filtered = len(heavy_filtered)
    peptide_ratio = (n_filtered / n_total * 100) if n_total > 0 else 0

    heavy_stats_file = os.path.join(output_dir, "heavy_cluster_stats.txt")
    with open(heavy_stats_file, "w") as f:
        f.write("Heavy Chain ClusterID Plot Statistics\n")
        f.write("---------------------------------------\n")
        f.write(
            f"Total Precursor Abundance in filtered heavy chain (cdr3_coverage_percentage > 0): {total_cluster_abundance}\n"
        )
        f.write(
            f"Total Precursor Abundance in entire heavy chain dataset: {total_heavy_abundance}\n"
        )
        f.write(
            f"Percentage of filtered total vs. entire heavy chain: {abundance_ratio:.2f}%\n\n"
        )
        f.write(f"Number of peptides with cdr3_coverage_percentage > 0: {n_filtered}\n")
        f.write(f"Total number of peptides in heavy chain dataset: {n_total}\n")
        f.write(f"Percentage of peptides used: {peptide_ratio:.2f}%\n")
    logging.info("Saved heavy chain cluster stats to " + heavy_stats_file)

def plot_gene_usage(chain_name, chain_df, gene_col, output_dir):
    """Generate a gene usage bar plot and companion stats file for a given chain and gene column."""
    # Make a copy to avoid modifying the original DF
    chain_df = chain_df.copy()

    # Process the gene column: split on comma and take the first element
    chain_df[gene_col] = chain_df[gene_col].apply(
        lambda x: x.split(',')[0].strip() if isinstance(x, str) else x
    )

    # Filter rows where the gene column is not empty or null
    gene_df = chain_df[chain_df[gene_col].notnull() & (chain_df[gene_col] != "")]
    if gene_df.empty:
        logging.info(f"No rows with non-empty {gene_col} for {chain_name} chain; skipping {gene_col} plot.")
        return

    # Group by the gene column and sum Precursor Abundance
    gene_group = gene_df.groupby(gene_col)["Precursor Abundance"].sum().reset_index()
    total_gene_abundance = gene_group["Precursor Abundance"].sum()
    gene_group["percent_abundance"] = (
        gene_group["Precursor Abundance"] / total_gene_abundance
    ) * 100

    # Order the data by percent_abundance descending
    gene_group = gene_group.sort_values(by="percent_abundance", ascending=False)

    # Create bar plot using seaborn with ordered x-axis and vertical labels
    plt.figure(figsize=(10, 6))
    sns.barplot(
        x=gene_col,
        y="percent_abundance",
        data=gene_group,
        order=gene_group[gene_col].tolist()
    )
    plt.ylabel("Percent Abundance")
    plt.title(f"Relative Abundance by {gene_col} for {chain_name} chain")
    plt.xticks(rotation=90)
    plot_filename = os.path.join(output_dir, f"{chain_name}_{gene_col}_plot.png")
    plt.tight_layout()
    plt.savefig(plot_filename, dpi=600)
    plt.close()
    logging.info(f"Saved {chain_name} {gene_col} plot to " + plot_filename)

    # Save the aggregated data (x-axis & y-axis) as a companion TSV
    plot_data = gene_group[[gene_col, "percent_abundance"]]
    plot_data_file = os.path.join(output_dir, f"{chain_name}_{gene_col}_plot_data.tsv")
    plot_data.to_csv(plot_data_file, sep="\t", index=False)
    logging.info(f"Saved {chain_name} {gene_col} plot data to " + plot_data_file)

    # Write companion TXT file with statistics
    total_peptides = len(chain_df)
    peptides_used = len(gene_df)
    peptide_percentage = (peptides_used / total_peptides * 100) if total_peptides > 0 else 0
    chain_total_abundance = chain_df["Precursor Abundance"].sum()
    abundance_percentage = (
        total_gene_abundance / chain_total_abundance * 100
        if chain_total_abundance > 0
        else 0
    )

    stats_filename = os.path.join(output_dir, f"{chain_name}_{gene_col}_stats.txt")
    with open(stats_filename, "w") as f:
        f.write(f"{chain_name.capitalize()} Chain {gene_col.upper()} Plot Statistics\n")
        f.write("----------------------------------------\n")
        f.write(f"Total peptides for {chain_name} chain: {total_peptides}\n")
        f.write(f"Peptides used for {gene_col} plot: {peptides_used} ({peptide_percentage:.2f}%)\n")
        f.write(f"Total Precursor Abundance for plotted subset: {total_gene_abundance}\n")
        f.write(f"Total Precursor Abundance for {chain_name} chain: {chain_total_abundance}\n")
        f.write(f"Percentage of plotted subset's Precursor Abundance: {abundance_percentage:.2f}%\n")
    logging.info(f"Saved {chain_name} {gene_col} stats to " + stats_filename)

def main():
    """Main function to parse arguments and run the plotting procedures."""
    parser = argparse.ArgumentParser(
        description="Plot relative quantities of B cell lineages and gene usage from antibody mass spec data merged with BCR-seq data."
    )
    parser.add_argument("--input_file", required=True,
                        help="Path to input TSV file containing antibody mass spec and BCR-seq data.")
    args = parser.parse_args()

    # Create output directory and set up logging
    output_dir = create_output_dir(args.input_file)
    setup_logging(output_dir, sys.argv)

    # Load data
    df = load_data(args.input_file)

    # Split into heavy, kappa, lambda, and light chain subsets (using v_call)
    heavy_df, kappa_df, lambda_df, light_df = split_chains(df)

    # Generate heavy chain ClusterID plot (if heavy_df is non-empty)
    if not heavy_df.empty:
        plot_heavy_cluster(heavy_df, output_dir)
    else:
        logging.info("No heavy chain data available; skipping heavy cluster plot.")

    # For gene usage plots, we handle heavy, kappa, and light chain subsets separately.
    # Now using 'v_call' and 'j_call' instead of 'v_gene' and 'j_gene'.
    chain_dict = {
        "heavy": heavy_df,
        "kappa": kappa_df,
        "light": light_df
    }
    for chain_name, chain_df in chain_dict.items():
        if chain_df.empty:
            logging.info(f"No {chain_name} chain data available; skipping gene usage plots for {chain_name}.")
            continue
        # v_call plot
        plot_gene_usage(chain_name, chain_df, "v_call", output_dir)
        # j_call plot
        plot_gene_usage(chain_name, chain_df, "j_call", output_dir)

if __name__ == "__main__":
    main()
