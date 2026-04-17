#!/usr/bin/env python3
"""
Plot lineage abundance and per-lineage BCR sequence/peptide coverage summaries.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import pandas as pd


REGIONS = ["fwr1_aa", "cdr1_aa", "fwr2_aa", "cdr2_aa", "fwr3_aa", "cdr3_aa", "fwr4_aa"]
AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")
MIN_LOGO_FIGURE_WIDTH = 18
LOGO_WIDTH_PER_RESIDUE = 0.22
LOGO_FIGURE_HEIGHT = 4.5
LOGO_COVERAGE_HEIGHT_RATIO = 1.2
LOGO_SEQUENCE_HEIGHT_RATIO = 1.7


def timestamp() -> str:
    return datetime.now().strftime("%Y-%m-%d_%H-%M-%S")


def make_out_dir(input_path: Path, out_dir: str | None) -> Path:
    if out_dir:
        path = Path(out_dir).expanduser()
    else:
        path = input_path.parent / f"plot_lineage_repertoire_{timestamp()}"
    path.mkdir(parents=True, exist_ok=True)
    return path


def read_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    df.columns = df.columns.str.strip('"')
    return df


def to_number(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def build_full_sequence(row: pd.Series) -> tuple[str, list[dict[str, Any]]]:
    seqs: list[str] = []
    positions: list[dict[str, Any]] = []
    start = 0
    for region in REGIONS:
        seq = str(row.get(region, "") if pd.notna(row.get(region, "")) else "").upper()
        seqs.append(seq)
        end = start + len(seq)
        positions.append({"region": region, "start": start, "end": end})
        start = end
    return "".join(seqs), positions


def plot_lineage_abundance(lineage: pd.DataFrame, out_path: Path, threshold: float) -> None:
    plot_df = lineage[lineage["relative_abundance"] >= threshold].copy()
    plot_df = plot_df.sort_values("relative_abundance", ascending=False)
    plt.figure(figsize=(max(8, len(plot_df) * 0.35), 5))
    plt.bar(plot_df["ClusterID"].astype(str), plot_df["relative_abundance"], color="#4C78A8")
    plt.ylabel("Relative abundance")
    plt.xlabel("ClusterID")
    plt.title(f"Lineage abundance (relative abundance >= {threshold:.3g})")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()


def aa_frequency_matrix(seqs: list[str]) -> pd.DataFrame:
    max_len = max((len(s) for s in seqs), default=0)
    rows = []
    for pos in range(max_len):
        counts = Counter(s[pos] for s in seqs if pos < len(s) and s[pos] in AA_ORDER)
        total = sum(counts.values()) or 1
        row = {aa: counts.get(aa, 0) / total for aa in AA_ORDER}
        row["position"] = pos + 1
        rows.append(row)
    return pd.DataFrame(rows)


def require_logomaker():
    try:
        import logomaker  # type: ignore
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "plot_lineage_repertoire.py now requires the 'logomaker' package for "
            "letter-based sequence logos. Install/update the conda environment "
            "from environment.yml or run: conda install -c conda-forge logomaker"
        ) from exc
    return logomaker


def plot_lineage_detail(
    cluster_id: str,
    bcr_rows: pd.DataFrame,
    peptide_rows: pd.DataFrame,
    out_path: Path,
) -> dict[str, Any]:
    seqs_and_pos = [build_full_sequence(row) for _, row in bcr_rows.iterrows()]
    seqs = [sp[0] for sp in seqs_and_pos if sp[0]]
    if not seqs:
        return {"cluster_id": cluster_id, "warning": "No BCR sequences available."}

    logomaker = require_logomaker()
    lengths = sorted({len(s) for s in seqs})
    ref_seq, ref_positions = seqs_and_pos[0]
    coverage = [0.0] * len(ref_seq)

    for _, row in peptide_rows.iterrows():
        if str(row.get("match_found", "")).lower() not in {"true", "1"}:
            continue
        try:
            start = int(float(row["match_start"]))
            end = int(float(row["match_end"]))
        except Exception:
            continue
        abundance = float(row.get("peptide_abundance", 0) or 0)
        for idx in range(max(0, start), min(len(coverage), end)):
            coverage[idx] += abundance

    freq = aa_frequency_matrix(seqs)
    logo_matrix = freq.set_index("position")[AA_ORDER]
    coverage_df = pd.DataFrame(
        {
            "position": list(range(1, len(coverage) + 1)),
            "reference_aa": list(ref_seq),
            "peptide_abundance_coverage": coverage,
        }
    )
    for comp in ref_positions:
        coverage_df.loc[
            (coverage_df["position"] >= comp["start"] + 1) & (coverage_df["position"] <= comp["end"]),
            "region",
        ] = comp["region"]

    matrix_path = out_path.with_name(out_path.stem + "_logo_matrix.tsv")
    coverage_path = out_path.with_name(out_path.stem + "_coverage.tsv")
    sequences_path = out_path.with_name(out_path.stem + "_sequences.tsv")
    freq.to_csv(matrix_path, sep="\t", index=False)
    coverage_df.to_csv(coverage_path, sep="\t", index=False)
    pd.DataFrame(
        {
            "ClusterID": cluster_id,
            "sequence_index": list(range(1, len(seqs) + 1)),
            "full_bcr_seq": seqs,
            "sequence_length": [len(s) for s in seqs],
        }
    ).to_csv(sequences_path, sep="\t", index=False)

    figure_width = max(MIN_LOGO_FIGURE_WIDTH, len(ref_seq) * LOGO_WIDTH_PER_RESIDUE)
    fig, axes = plt.subplots(
        2,
        1,
        figsize=(figure_width, LOGO_FIGURE_HEIGHT),
        sharex=True,
        gridspec_kw={"height_ratios": [LOGO_COVERAGE_HEIGHT_RATIO, LOGO_SEQUENCE_HEIGHT_RATIO]},
    )
    axes[0].bar(range(1, len(coverage) + 1), coverage, width=1.0, color="#F58518")
    axes[0].set_ylabel("Peptide abundance")
    axes[0].set_title(f"Cluster {cluster_id}: peptide coverage and BCR sequence logo")

    logomaker.Logo(logo_matrix, ax=axes[1])
    axes[1].set_ylabel("AA frequency")
    axes[1].set_xlabel("Residue position")
    tick_step = max(1, len(ref_seq) // 20)
    tick_positions = list(range(0, len(ref_seq), tick_step))
    axes[1].set_xticks(tick_positions)
    axes[1].set_xticklabels([str(i + 1) for i in tick_positions], rotation=0)

    ymax = max(coverage) * 1.05 if coverage and max(coverage) > 0 else 1
    for comp in ref_positions:
        if comp["end"] <= comp["start"]:
            continue
        mid = (comp["start"] + comp["end"]) / 2 + 1
        axes[0].axvspan(comp["start"] + 1, comp["end"], alpha=0.08, color="#54A24B")
        axes[0].text(mid, ymax, comp["region"].replace("_aa", "").upper(), ha="center", va="top", fontsize=8)
    axes[0].set_ylim(0, ymax * 1.15)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close(fig)
    return {
        "cluster_id": cluster_id,
        "bcr_sequence_count": int(len(seqs)),
        "sequence_lengths": lengths,
        "logo_matrix_tsv": str(matrix_path),
        "coverage_tsv": str(coverage_path),
        "sequences_tsv": str(sequences_path),
        "warning": "Sequence lengths differ within lineage." if len(lengths) > 1 else "",
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot quantified lineage abundance and per-lineage coverage summaries."
    )
    parser.add_argument("mapped_peptides", help="Mapped peptide TSV from quantify_map_peptides.py.")
    parser.add_argument("bcrseq_tsv", help="BCRseq annotation TSV.")
    parser.add_argument("--min-cdr3-overlap-aa", type=int, default=0)
    parser.add_argument("--relative-abundance-threshold", type=float, default=0.01)
    parser.add_argument("--out-dir", help="Output directory. Default: timestamped directory next to mapped peptides.")
    args = parser.parse_args()

    mapped_path = Path(args.mapped_peptides).expanduser()
    bcr_path = Path(args.bcrseq_tsv).expanduser()
    out_dir = make_out_dir(mapped_path, args.out_dir)
    stem = mapped_path.stem

    mapped = read_table(mapped_path)
    bcr = read_table(bcr_path)
    mapped["peptide_abundance_num"] = to_number(mapped["peptide_abundance"]).fillna(0.0)
    mapped["cdr3_cover_aa_num"] = to_number(mapped["cdr3_cover_aa"]).fillna(0)
    mapped["match_start_num"] = to_number(mapped["match_start"])
    mapped["match_end_num"] = to_number(mapped["match_end"])

    retained = mapped[
        (mapped["match_found"].astype(str).str.lower().isin(["true", "1"]))
        & (mapped["cdr3_cover_aa_num"] >= args.min_cdr3_overlap_aa)
    ].copy()

    abundance_units = retained.drop_duplicates(["ClusterID", "peptide_sequence"]).copy()
    lineage = (
        abundance_units.groupby("ClusterID", dropna=False)["peptide_abundance_num"]
        .sum()
        .reset_index(name="lineage_abundance")
    )
    total = lineage["lineage_abundance"].sum()
    lineage["relative_abundance"] = lineage["lineage_abundance"] / total if total else 0.0
    lineage = lineage.sort_values("relative_abundance", ascending=False)

    lineage_tsv = out_dir / f"{stem}_lineage_abundance.tsv"
    lineage_png = out_dir / f"{stem}_lineage_abundance.png"
    lineage_plot_tsv = out_dir / f"{stem}_lineage_abundance_plotted.tsv"
    lineage.to_csv(lineage_tsv, sep="\t", index=False)
    lineage[lineage["relative_abundance"] >= args.relative_abundance_threshold].to_csv(
        lineage_plot_tsv, sep="\t", index=False
    )
    plot_lineage_abundance(lineage, lineage_png, args.relative_abundance_threshold)

    detail_dir = out_dir / "lineage_details"
    detail_dir.mkdir(exist_ok=True)
    selected = lineage[lineage["relative_abundance"] >= args.relative_abundance_threshold]["ClusterID"].astype(str).tolist()
    detail_logs = []
    for cluster_id in selected:
        bcr_rows = bcr[bcr["ClusterID"].astype(str) == cluster_id].copy()
        pep_rows = retained[retained["ClusterID"].astype(str) == cluster_id].copy()
        pep_rows = pep_rows.drop_duplicates(["ClusterID", "peptide_sequence", "match_start", "match_end"])
        safe_id = re.sub(r"[^A-Za-z0-9_.-]+", "_", cluster_id)
        detail_logs.append(plot_lineage_detail(cluster_id, bcr_rows, pep_rows, detail_dir / f"cluster_{safe_id}_coverage_logo.png"))

    log = {
        "command": " ".join(sys.argv),
        "mapped_peptides": str(mapped_path),
        "bcrseq_tsv": str(bcr_path),
        "min_cdr3_overlap_aa": args.min_cdr3_overlap_aa,
        "relative_abundance_threshold": args.relative_abundance_threshold,
        "input_mapped_rows": int(len(mapped)),
        "retained_rows": int(len(retained)),
        "retained_abundance": float(abundance_units["peptide_abundance_num"].sum()),
        "selected_lineages": selected,
        "lineage_detail_notes": detail_logs,
        "outputs": {
            "lineage_tsv": str(lineage_tsv),
            "lineage_plot_tsv": str(lineage_plot_tsv),
            "lineage_png": str(lineage_png),
            "lineage_details": str(detail_dir),
        },
    }
    log_path = out_dir / f"{stem}_plot_lineage_repertoire_log_{timestamp()}.json"
    log_path.write_text(json.dumps(log, indent=2))

    print(f"Wrote lineage abundance TSV: {lineage_tsv}")
    print(f"Wrote lineage abundance plot: {lineage_png}")
    print(f"Wrote lineage detail plots: {detail_dir}")
    print(f"Wrote log: {log_path}")


if __name__ == "__main__":
    main()
