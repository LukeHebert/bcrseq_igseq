#!/usr/bin/env python3
"""
Visualize BCR peptide mapping overlap percentages by CDR region and clonal family counts.

This script reads one or more TSV files containing merged IgBLAST and Proteome Discoverer annotations,
parses the `region_overlap_percentages` field for CDR regions, and produces three plots:
1. `<prefix>_box.png`    – box-and-whisker of overlap % by CDR region (with colored outliers).
2. `<prefix>_bar.png`    – bar-plot of % of total peptides mapping to each CDR region.
3. `<prefix>_clonal.png` – bar-plot of number of unique Clonal Families (ClusterID) with cdr3 > 0.

Titles (via `--titles`) determine the legend labels, and `-o PREFIX` sets the filename stem.
Each source always gets the same color in all plots.
"""

import argparse
from pathlib import Path

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.collections as mcoll
import numpy as np


def parse_args():
    """Parse command-line arguments for input TSV files, optional titles, and output prefix."""
    parser = argparse.ArgumentParser(
        description='Visualize BCR peptide mapping overlap and clonal families.'
    )
    parser.add_argument(
        'input_files', nargs='+',
        help='One or more TSV files with region_overlap_percentages and ClusterID.'
    )
    parser.add_argument(
        '--titles', nargs='+',
        help='Optional labels for each input file, in the same order.'
    )
    parser.add_argument(
        '-o', '--output', required=True,
        help='Output file prefix (e.g. "plot" → plot_box.png, plot_bar.png, plot_clonal.png).'
    )
    return parser.parse_args()


def load_data(file_paths, titles=None):
    """Load and concatenate TSVs, tagging each row with its source title."""
    frames = []
    if titles and len(titles) != len(file_paths):
        raise ValueError("Number of titles must match number of input files")
    for idx, fp in enumerate(file_paths):
        df = pd.read_csv(fp, sep='\t')
        source = titles[idx] if titles else Path(fp).stem
        df['source'] = source
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def parse_region_overlap(df):
    """Explode and parse region_overlap_percentages into tidy rows."""
    tmp = df[['source', 'region_overlap_percentages']].copy()
    tmp['kv'] = tmp['region_overlap_percentages'].str.split(',')
    tmp = tmp.explode('kv').dropna(subset=['kv'])
    # FIXED SPLIT:
    kv_split = tmp['kv'].str.split(':', n=1, expand=True)
    tmp['region']     = kv_split[0]
    tmp['percentage'] = kv_split[1].str.rstrip('%').astype(float)
    return tmp[tmp['region'].isin(['cdr1_aa','cdr2_aa','cdr3_aa'])][['source','region','percentage']]



def plot_boxplot(raw_df, tidy_df, prefix, sources, palette_dict):
    """Box-and-whisker of overlap % by CDR, with colored outliers."""
    sns.set(style='whitegrid')
    plt.figure(figsize=(8,6))

    flier_props = dict(marker='o', markersize=4, linestyle='none')
    ax = sns.boxplot(
        data=tidy_df, x='region', y='percentage',
        hue='source', hue_order=sources,
        palette=palette_dict, flierprops=flier_props, dodge=True
    )

    # recolor outliers
    fliers = [c for c in ax.collections if isinstance(c, mcoll.PathCollection)]
    for idx, coll in enumerate(fliers):
        src = sources[idx % len(sources)]
        c = palette_dict[src]
        coll.set_edgecolor(c)
        coll.set_facecolor(c)
        coll.set_linewidth(0.5)

    ax.set_ylabel('Overlap Percentage')
    ax.set_xticklabels(['CDR1','CDR2','CDR3'])
    ax.set_title('Distribution of Peptide Mapping Overlap by CDR Region')
    ax.legend(title='Source Files', bbox_to_anchor=(1.05,1), loc='upper left')
    plt.tight_layout()

    out = f"{prefix}_box.png"
    plt.savefig(out, dpi=600)
    print(f"Boxplot saved to {out}")
    plt.close()


def plot_barplot(raw_df, tidy_df, prefix, sources, palette_dict):
    """Bar-plot of % of total peptides mapping to each CDR region."""
    sns.set(style='whitegrid')
    total = raw_df.groupby('source').size()
    counts = tidy_df.groupby(['source','region']).size()
    pct = (counts/total * 100).round().astype(int).reset_index(name='mapping_pct')

    plt.figure(figsize=(8,6))
    ax = sns.barplot(
        data=pct, x='region', y='mapping_pct',
        hue='source', hue_order=sources,
        palette=palette_dict
    )
    ax.set_ylabel('Mapping Frequency (%)')
    ax.set_xticklabels(['CDR1','CDR2','CDR3'])
    ax.set_title('Percentage of Total Peptide Mappings per CDR Region')
    ax.legend(title='Source Files', bbox_to_anchor=(1.05,1), loc='upper left')
    plt.tight_layout()

    out = f"{prefix}_bar.png"
    plt.savefig(out, dpi=600)
    print(f"Barplot saved to {out}")
    plt.close()


def plot_clonal_barplot(raw_df, tidy_df, prefix, sources, palette_dict):
    """Bar-plot of # unique Clonal Families with cdr3 coverage > 0 per source."""
    sns.set(style='whitegrid')

    # find indices where cdr3_aa > 0
    idx = tidy_df.loc[(tidy_df.region=='cdr3_aa') & (tidy_df.percentage>0)].index.unique()
    sub = raw_df.loc[idx]
    fam_counts = sub.groupby('source')['ClusterID'].nunique()
    # ensure all sources appear
    fam_counts = fam_counts.reindex(sources, fill_value=0).reset_index(name='n_families')

    plt.figure(figsize=(8,6))
    ax = sns.barplot(
        data=fam_counts, x='source', y='n_families',
        order=sources, palette=palette_dict
    )
    ax.set_ylabel('Number of Clonal Families (cdr3 > 0)')
    ax.set_xlabel('')
    ax.set_title('Clonal Families with cdr3 Coverage > 0')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()

    out = f"{prefix}_clonal.png"
    plt.savefig(out, dpi=600)
    print(f"Clonal families barplot saved to {out}")
    plt.close()


def main():
    args = parse_args()
    raw_df  = load_data(args.input_files, args.titles)
    tidy_df = parse_region_overlap(raw_df)
    prefix  = Path(args.output).stem

    # derive a consistent source order & palette
    sources     = raw_df['source'].drop_duplicates().tolist()
    palette_list = sns.color_palette('colorblind', n_colors=len(sources))
    palette_dict = dict(zip(sources, palette_list))

    plot_boxplot(raw_df, tidy_df, prefix, sources, palette_dict)
    plot_barplot(raw_df, tidy_df, prefix, sources, palette_dict)
    plot_clonal_barplot(raw_df, tidy_df, prefix, sources, palette_dict)


if __name__ == '__main__':
    main()
