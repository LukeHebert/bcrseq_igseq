#!/usr/bin/env python3
"""
Process one or more antibody-mass-spec TSVs merged with BCR-seq to produce
either clustered bar plots (by file) or clustered boxplots (by user-defined groups)
for:
  - Heavy-chain lineage by rank (clonal families)
  - Gene usage (v_call, j_call) across heavy, kappa, light chains

Usage:
  # Bar plots by file:
  python plot_igseq.py --input_files a.tsv b.tsv … [--split_names] [--max_cf N] [--max_v N] [--max_j M]

  # Boxplots by group:
  python plot_igseq.py \
    --groups Group1:a1.tsv,a2.tsv Group2:b1.tsv … \
    [--max_cf N] [--max_v N] [--max_j M]
"""
import os
import sys
import argparse
import logging
from datetime import datetime

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def create_output_dir(first_input):
    """Create a single 'plot_igseq' directory next to the first input file."""
    parent = os.path.dirname(os.path.abspath(first_input))
    out = os.path.join(parent, "plot_igseq")
    os.makedirs(out, exist_ok=True)
    return out

def setup_logging(output_dir, cmd):
    """Timestamped log file in output_dir recording the command and outputs."""
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    logf = os.path.join(output_dir, f"log_{ts}.txt")
    logging.basicConfig(
        filename=logf,
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s"
    )
    logging.info("Command: " + " ".join(cmd))
    logging.info("Output dir: " + output_dir)
    return logf

def split_chains(df):
    """
    Split on v_call into heavy (IGHV), kappa (IGKV), lambda (IGLV),
    plus combined light = kappa+lambda.
    """
    hv = df[df['v_call'].str.contains("IGHV", na=False)]
    kp = df[df['v_call'].str.contains("IGKV", na=False)]
    lv = df[df['v_call'].str.contains("IGLV", na=False)]
    lt = pd.concat([kp, lv], ignore_index=True)
    logging.info(
        f"Chains: heavy={len(hv)}, kappa={len(kp)}, "
        f"lambda={len(lv)}, light={len(lt)}"
    )
    return hv, kp, lv, lt

def gather_heavy_lineages(label, df):
    """Compute percent_abundance by ClusterID→ranked for one label (file or group member)."""
    sub = df[df['cdr3_coverage_percentage'] > 0]
    grp = (sub.groupby("ClusterID")["Precursor Abundance"]
              .sum().reset_index(name="abundance"))
    tot = grp["abundance"].sum()
    grp["percent_abundance"] = grp["abundance"] / tot * 100
    grp = grp.sort_values("percent_abundance", ascending=False).reset_index(drop=True)
    grp["rank"] = grp.index + 1
    return grp[["rank", "percent_abundance"]].assign(label=label)

def gather_gene_usage(label, df, gene_col):
    """Compute percent_abundance by gene_call (first allele) for one label."""
    d = df.copy()
    d[gene_col] = d[gene_col].apply(
        lambda x: x.split(",")[0].strip() if isinstance(x, str) else x
    )
    d = d[d[gene_col].notna() & (d[gene_col] != "")]
    grp = (d.groupby(gene_col)["Precursor Abundance"]
            .sum().reset_index(name="abundance"))
    tot = grp["abundance"].sum()
    grp["percent_abundance"] = grp["abundance"] / tot * 100
    return grp[[gene_col, "percent_abundance"]] \
              .rename(columns={gene_col:"gene_call"}) \
              .assign(label=label)

def parse_groups(group_args):
    """
    Parse --groups of the form Name:file1.tsv,file2.tsv ...
    Returns dict name→[filepaths]. Allows even single-file groups.
    """
    groups = {}
    for g in group_args:
        if ":" not in g:
            raise ValueError(f"Invalid group spec: {g}")
        name, files = g.split(":",1)
        flist = files.split(",")
        groups[name] = flist
    if len(groups) < 2:
        raise ValueError("Must specify ≥2 groups for boxplots")
    return groups

def plot_bar(df, x, y, hue, order, xlabel, ylabel, title, out_png, out_tsv, split_names):
    """Clustered bar plot by file."""
    if split_names:
        df[hue] = df[hue].apply(
            lambda f: os.path.splitext(os.path.basename(f))[0].split("_")[0]
        )
    plt.figure(figsize=(10,6))
    sns.barplot(data=df, x=x, y=y, hue=hue, palette="colorblind", order=order)
    plt.title(title); plt.xlabel(xlabel); plt.ylabel(ylabel)
    plt.xticks(rotation=90)
    plt.legend(title="File", bbox_to_anchor=(1.02,1), loc="upper left")
    plt.tight_layout()
    plt.savefig(out_png, dpi=600); plt.close()
    df[[x,hue,y]].to_csv(out_tsv, sep="\t", index=False)
    logging.info(f"Wrote {out_png} + {out_tsv}")

def plot_box(df, x, y, hue, order, xlabel, ylabel, title, out_png, out_tsv):
    """Clustered box‐and‐whisker plot by group, overlaid with points so single‐value groups show color."""
    # Colorblind palette; blue‐ish & orange‐ish for first two groups, etc.
    pal = sns.color_palette("colorblind", n_colors=df[hue].nunique())

    plt.figure(figsize=(10,6))
    ax = sns.boxplot(
        data=df, x=x, y=y, hue=hue,
        palette=pal, order=order
    )
    # overlay the actual points
    sns.stripplot(
        data=df, x=x, y=y, hue=hue,
        palette=pal, dodge=True,
        size=6, edgecolor="gray", linewidth=0.5,
        ax=ax, legend=False
    )

    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xticks(rotation=90)
    # only one legend (from boxplot)
    plt.legend(title="Group", bbox_to_anchor=(1.02,1), loc="upper left")
    plt.tight_layout()
    plt.savefig(out_png, dpi=600)
    plt.close()

    # companion TSV
    df[[x, hue, y]].to_csv(out_tsv, sep="\t", index=False)
    logging.info(f"Wrote {out_png} + {out_tsv}")


def main():
    p = argparse.ArgumentParser(
        description="Clustered barplots or boxplots of heavy lineage & gene usage"
    )
    p.add_argument("--input_files", nargs="*", help="TSV files (bar mode)")
    p.add_argument(
        "--groups", nargs="*", metavar="Name:f1,f2,…",
        help="Groups for boxplot mode: Name:file1.tsv,… (≥2 groups)"
    )
    p.add_argument(
        "--split_names", action="store_true",
        help="In bar mode, abbreviate file labels by splitting on '_'"
    )
    p.add_argument("--max_cf", type=int, default=None,
                   help="Max number of clonal families (lineage ranks) to plot")
    p.add_argument("--max_v", type=int, default=None, help="Max v_call genes to plot")
    p.add_argument("--max_j", type=int, default=None, help="Max j_call genes to plot")
    args = p.parse_args()

    # Choose mode
    if args.groups:
        groups = parse_groups(args.groups)
        all_files = [f for fl in groups.values() for f in fl]
        bar_mode = False
    else:
        if not args.input_files:
            p.error("Provide either --input_files (bar mode) or --groups (box mode)")
        all_files = args.input_files
        groups = None
        bar_mode = True

    outdir = create_output_dir(all_files[0])
    setup_logging(outdir, sys.argv)

    # Collect
    lineage_dfs = []
    gene_dfs = { gc:{ch:[] for ch in ('heavy','kappa','light')} for gc in ('v_call','j_call') }

    # In box mode we’ll gather per-file then re-label as group below
    for f in all_files:
        df = pd.read_csv(f, sep="\t")
        hv, kp, lv, lt = split_chains(df)

        if bar_mode:
            lineage_dfs.append(gather_heavy_lineages(f, hv))
            for chain, sub in (('heavy',hv),('kappa',kp),('light',lt)):
                for gc in ('v_call','j_call'):
                    gene_dfs[gc][chain].append(gather_gene_usage(f, sub, gc))
        else:
            # still gather per-file for grouping
            lineage_dfs.append(gather_heavy_lineages(f, hv).assign(file=f))
            for chain,sub in (('heavy',hv),('kappa',kp),('light',lt)):
                for gc in ('v_call','j_call'):
                    gene_dfs[gc][chain].append(
                        gather_gene_usage(f, sub, gc).assign(file=f)
                    )

    # ==== BAR MODE ====
    if bar_mode:
        # --- Heavy lineage bar ---
        all_lin = pd.concat(lineage_dfs, ignore_index=True)
        # enforce full (file,rank) grid
        files = all_lin['label'].unique()
        ranks = sorted(all_lin['rank'].unique())
        if args.max_cf:
            ranks = ranks[:args.max_cf]
            all_lin = all_lin[all_lin['rank'].isin(ranks)]
        idx = pd.MultiIndex.from_product([files, ranks], names=['label','rank'])
        all_lin = (all_lin.set_index(['label','rank'])
                         .reindex(idx, fill_value=0)
                         .reset_index())
        plot_bar(
            all_lin, 'rank','percent_abundance','label', ranks,
            'Lineage rank','Percent abundance',
            'Heavy‐chain lineage by rank',
            os.path.join(outdir,'heavy_lineage_bar.png'),
            os.path.join(outdir,'heavy_lineage_bar_data.tsv'),
            args.split_names
        )

        # --- Gene usage bars ---
        for gc, cap in (('v_call',args.max_v), ('j_call',args.max_j)):
            for chain in ('heavy','kappa','light'):
                parts = gene_dfs[gc][chain]
                if not parts:
                    logging.info(f"Skip bar {chain}/{gc}, no data")
                    continue
                df_long = pd.concat(parts, ignore_index=True)
                avg = (df_long.groupby('gene_call')['percent_abundance']
                              .mean()
                              .reset_index()
                              .sort_values('percent_abundance',ascending=False))
                genes = avg['gene_call'].tolist()[:cap] if cap else avg['gene_call'].tolist()
                idx2 = pd.MultiIndex.from_product(
                    [df_long['label'].unique(), genes],
                    names=['label','gene_call']
                )
                df_long = (df_long.set_index(['label','gene_call'])
                                   .reindex(idx2, fill_value=0)
                                   .reset_index())
                plot_bar(
                    df_long,'gene_call','percent_abundance','label', genes,
                    gc,'Percent abundance',
                    f'{chain.capitalize()}‐chain {gc} usage',
                    os.path.join(outdir,f'{chain}_{gc}_bar.png'),
                    os.path.join(outdir,f'{chain}_{gc}_bar_data.tsv'),
                    args.split_names
                )

    # ==== BOX MODE ====
    else:
        # --- Heavy lineage box ---
        recs = []
        for grp_name, files in groups.items():
            for f in files:
                tmp = lineage_dfs.pop(0)  # in same order
                # tmp has columns rank,percent_abundance,label
                dfb = tmp.rename(columns={'percent_abundance':'val','label':'file'})
                dfb['group'] = grp_name
                recs.append(dfb[['rank','val','group']])
        df_box = pd.concat(recs, ignore_index=True)
        if args.max_cf:
            df_box = df_box[df_box['rank'] <= args.max_cf]
        order = sorted(df_box['rank'].unique())
        plot_box(
            df_box,'rank','val','group', order,
            'Lineage rank','Percent abundance',
            'Heavy‐chain lineage by rank (box)',
            os.path.join(outdir,'heavy_lineage_box.png'),
            os.path.join(outdir,'heavy_lineage_box_data.tsv')
        )

        # --- Gene usage boxes ---
        for gc, cap in (('v_call',args.max_v), ('j_call',args.max_j)):
            for chain in ('heavy','kappa','light'):
                recs = []
                for grp_name, files in groups.items():
                    for _ in files:
                        dfc = gene_dfs[gc][chain].pop(0)  # same order
                        dfb = dfc.rename(columns={'percent_abundance':'val','label':'file'})
                        dfb['group'] = grp_name
                        recs.append(dfb[['gene_call','val','group']])
                if not recs:
                    logging.info(f"Skip box {chain}/{gc}, no data")
                    continue
                dfg = pd.concat(recs, ignore_index=True)
                avg = (dfg.groupby('gene_call')['val']
                          .mean().reset_index()
                          .sort_values('val', ascending=False))
                genes = avg['gene_call'].tolist()[:cap] if cap else avg['gene_call'].tolist()
                dfg = dfg[dfg['gene_call'].isin(genes)]
                plot_box(
                    dfg,'gene_call','val','group', genes,
                    gc,'Percent abundance',
                    f'{chain.capitalize()}‐chain {gc} usage (box)',
                    os.path.join(outdir,f'{chain}_{gc}_box.png'),
                    os.path.join(outdir,f'{chain}_{gc}_box_data.tsv')
                )

if __name__ == "__main__":
    main()
