#!/usr/bin/env python3
"""
Visualize BCR peptide mapping overlap percentages by CDR region.

This script reads one or more TSV files containing merged IgBLAST and Proteome Discoverer annotations,
parses the `region_overlap_percentages` field for CDR regions, and produces two plots:
1. A clustered box-and-whisker plot showing the distribution of overlap percentages by CDR region
   (with colored outlier markers matching each box).
2. A clustered bar plot showing the percentage of total peptide mappings per CDR region for each input file.

Specify custom titles for each input file via --titles, and provide an output file prefix via -o.
The script will save `<prefix>_box.png` for the boxplot and `<prefix>_bar.png` for the barplot.
"""

import argparse
from pathlib import Path

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.collections as mcoll
import numpy as np


def parse_args():
    """Parse command-line arguments for input TSV files, titles, and output prefix."""
    parser = argparse.ArgumentParser(
        description='Visualize BCR peptide mapping overlap percentages by CDR region.'
    )
    parser.add_argument(
        'input_files',
        nargs='+',
        help='One or more input TSV files containing a region_overlap_percentages column.'
    )
    parser.add_argument(
        '--titles',
        nargs='+',
        help='Optional titles for each input file, in the same order as input_files.'
    )
    parser.add_argument(
        '-o', '--output',
        help='Output file prefix (e.g. "plot" → generates plot_box.png and plot_bar.png).',
        required=True
    )
    return parser.parse_args()


def load_data(file_paths, titles=None):
    """Load and concatenate TSVs, tagging each row with its source title."""
    frames = []
    if titles and len(titles) != len(file_paths):
        raise ValueError("Number of titles must match number of input files")
    for idx, fp in enumerate(file_paths):
        df = pd.read_csv(fp, sep='\t')
        title = titles[idx] if titles else Path(fp).stem
        df['source'] = title
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def parse_region_overlap(df):
    """Transform region_overlap_percentages into a tidy DataFrame of region vs. percentage."""
    tmp = df[['source', 'region_overlap_percentages']].copy()
    tmp = tmp.assign(kv=tmp['region_overlap_percentages'].str.split(','))
    tmp = tmp.explode('kv').dropna(subset=['kv'])
    kv_split = tmp['kv'].str.split(':', expand=True)
    tmp['region'] = kv_split[0]
    tmp['percentage'] = kv_split[1].str.rstrip('%').astype(float)
    cdrs = ['cdr1_aa', 'cdr2_aa', 'cdr3_aa']
    return tmp[tmp['region'].isin(cdrs)][['source', 'region', 'percentage']]


def plot_boxplot(raw_df, tidy_df, output_prefix):
    """Create and save a box-and-whisker plot with colored outliers."""
    sns.set(style='whitegrid')
    plt.figure(figsize=(8, 6))
    sources = list(raw_df['source'].unique())
    palette = sns.color_palette('colorblind', n_colors=len(sources))

    # Draw boxplot with small flier markers
    flier_props = dict(marker='o', markersize=4, linestyle='none')
    ax = sns.boxplot(
        data=tidy_df,
        x='region',
        y='percentage',
        hue='source',
        dodge=True,
        palette=palette,
        flierprops=flier_props
    )

    # Recolor flier markers to match their box
    fliers = [c for c in ax.collections if isinstance(c, mcoll.PathCollection)]
    for idx, coll in enumerate(fliers):
        color = palette[idx % len(sources)]
        coll.set_edgecolor(color)
        coll.set_facecolor(color)
        coll.set_linewidth(0.5)

    # Styling
    plt.ylabel('Overlap Percentage')
    plt.xticks([0, 1, 2], ['CDR1', 'CDR2', 'CDR3'])
    plt.title('Distribution of Peptide Mapping Overlap by CDR Region')
    plt.legend(title='Source Files', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()

    box_file = f"{output_prefix}_box.png"
    plt.savefig(box_file, dpi=600)
    print(f"Boxplot saved to {box_file}")
    plt.close()


def plot_barplot(raw_df, tidy_df, output_prefix):
    """Create and save a clustered bar plot of mapping frequency percentages."""
    sns.set(style='whitegrid')

    # Calculate mapping frequencies per source-region
    total_counts = raw_df.groupby('source').size()
    region_counts = tidy_df.groupby(['source', 'region']).size()
    mapping_pct = (
        (region_counts / total_counts * 100)
        .round()
        .astype(int)
        .reset_index(name='mapping_pct')
    )

    plt.figure(figsize=(8, 6))
    palette = sns.color_palette('colorblind', n_colors=len(raw_df['source'].unique()))
    ax = sns.barplot(
        data=mapping_pct,
        x='region',
        y='mapping_pct',
        hue='source',
        palette=palette
    )
    plt.ylabel('Mapping Frequency (%)')
    plt.xticks([0, 1, 2], ['CDR1', 'CDR2', 'CDR3'])
    plt.title('Percentage of Total Peptide Mappings per CDR Region')
    plt.legend(title='Source Files', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()

    bar_file = f"{output_prefix}_bar.png"
    plt.savefig(bar_file, dpi=600)
    print(f"Barplot saved to {bar_file}")
    plt.close()


def main():
    """Execute data loading, parsing, and plotting workflow."""
    args = parse_args()
    raw_df = load_data(args.input_files, args.titles)
    tidy_df = parse_region_overlap(raw_df)
    prefix = Path(args.output).stem
    plot_boxplot(raw_df, tidy_df, prefix)
    plot_barplot(raw_df, tidy_df, prefix)


if __name__ == '__main__':
    main()
