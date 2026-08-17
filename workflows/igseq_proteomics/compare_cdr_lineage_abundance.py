#!/usr/bin/env python3
"""
Compare lineage abundance estimates derived from mapped peptide CDR coverage.

The input is the ``*_mapped_peptides.tsv`` output from
``quantify_map_peptides.py``.  Each quantitative peptide observation is counted
once per lineage, even when it maps to more than one BCR-seq member.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import warnings
from datetime import datetime
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import ConstantInputWarning, spearmanr


REGIONS = ["fwr1_aa", "cdr1_aa", "fwr2_aa", "cdr2_aa", "fwr3_aa", "cdr3_aa", "fwr4_aa"]
CDR_REGIONS = ["cdr1_aa", "cdr2_aa", "cdr3_aa"]
SELECTIONS = ["cdr1", "cdr2", "cdr1_cdr2", "cdr3", "all_peptides"]
SELECTION_LABELS = {
    "cdr1": "CDR1-only",
    "cdr2": "CDR2-only",
    "cdr1_cdr2": "Combined CDR1/CDR2",
    "cdr3": "CDR3",
    "all_peptides": "All peptides",
}
COMPARISONS = [
    ("cdr1", "cdr3"),
    ("cdr2", "cdr3"),
    ("cdr1_cdr2", "cdr3"),
    ("cdr1", "all_peptides"),
    ("cdr2", "all_peptides"),
    ("cdr1_cdr2", "all_peptides"),
    ("cdr3", "all_peptides"),
]
REQUIRED_COLUMNS = ["peptide_sequence", "ClusterID", "match_found", "match_start", "match_end"] + REGIONS


def timestamp() -> str:
    return datetime.now().strftime("%Y-%m-%d_%H-%M-%S")


def make_out_dir(input_path: Path, out_dir: str | None) -> Path:
    path = Path(out_dir).expanduser() if out_dir else input_path.parent / f"compare_cdr_lineage_abundance_{timestamp()}"
    path.mkdir(parents=True, exist_ok=True)
    return path


def read_table(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t", dtype=str)
    table.columns = table.columns.str.strip('"')
    return table


def require_columns(table: pd.DataFrame) -> str:
    missing = [column for column in REQUIRED_COLUMNS if column not in table.columns]
    if missing:
        raise ValueError("Mapped peptide TSV is missing required column(s): " + ", ".join(missing))
    if "peptide_abundance_total" in table.columns:
        return "peptide_abundance_total"
    if "peptide_abundance" in table.columns:
        return "peptide_abundance"
    raise ValueError("Mapped peptide TSV needs 'peptide_abundance_total' or 'peptide_abundance'.")


def as_match_found(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip().str.lower().isin({"true", "t", "1", "yes", "y"})


def region_coordinates(row: pd.Series) -> dict[str, tuple[int, int]]:
    """Return zero-based, half-open positions for all framework and CDR regions."""
    coordinates: dict[str, tuple[int, int]] = {}
    start = 0
    for region in REGIONS:
        sequence = "" if pd.isna(row[region]) else str(row[region]).strip()
        end = start + len(sequence)
        coordinates[region] = (start, end)
        start = end
    return coordinates


def overlap_length(start: int, end: int, region_start: int, region_end: int) -> int:
    return max(0, min(end, region_end) - max(start, region_start))


def classify_mapped_rows(table: pd.DataFrame, min_overlap: int) -> pd.DataFrame:
    """Add CDR overlap counts and qualification flags to successful mapping rows."""
    classified = table.loc[as_match_found(table["match_found"])].copy()
    classified["match_start"] = pd.to_numeric(classified["match_start"], errors="coerce")
    classified["match_end"] = pd.to_numeric(classified["match_end"], errors="coerce")
    if classified[["match_start", "match_end"]].isna().any().any():
        raise ValueError("Successful mapped rows must have numeric match_start and match_end values.")

    for index, row in classified.iterrows():
        coordinates = region_coordinates(row)
        start, end = int(row["match_start"]), int(row["match_end"])
        for region in CDR_REGIONS:
            region_start, region_end = coordinates[region]
            classified.loc[index, f"{region}_overlap_aa"] = overlap_length(start, end, region_start, region_end)

    for region in CDR_REGIONS:
        classified[f"{region}_qualifies"] = classified[f"{region}_overlap_aa"] >= min_overlap
    return classified


def collapse_quantitative_evidence(classified: pd.DataFrame, abundance_column: str) -> pd.DataFrame:
    """Keep one peptide abundance contribution per lineage and peptide sequence."""
    evidence = classified.copy()
    evidence["abundance"] = pd.to_numeric(evidence[abundance_column], errors="coerce").fillna(0.0)
    grouped = evidence.groupby(["ClusterID", "peptide_sequence"], dropna=False, sort=False)
    collapsed = grouped["abundance"].first().reset_index()
    for region in CDR_REGIONS:
        collapsed[f"{region}_qualifies"] = grouped[f"{region}_qualifies"].any().to_numpy()
    collapsed["cdr1"] = collapsed["cdr1_aa_qualifies"] & ~collapsed["cdr2_aa_qualifies"]
    collapsed["cdr2"] = collapsed["cdr2_aa_qualifies"] & ~collapsed["cdr1_aa_qualifies"]
    collapsed["cdr1_cdr2"] = collapsed["cdr1_aa_qualifies"] | collapsed["cdr2_aa_qualifies"]
    collapsed["cdr3"] = collapsed["cdr3_aa_qualifies"]
    collapsed["all_peptides"] = True
    return collapsed


def lineage_abundance_table(evidence: pd.DataFrame) -> pd.DataFrame:
    lineages = pd.DataFrame({"ClusterID": sorted(evidence["ClusterID"].dropna().astype(str).unique())})
    evidence = evidence.copy()
    evidence["ClusterID"] = evidence["ClusterID"].astype(str)
    for selection in SELECTIONS:
        selected = evidence.loc[evidence[selection]].groupby("ClusterID", as_index=False)["abundance"].sum()
        selected = selected.rename(columns={"abundance": f"{selection}_abundance"})
        lineages = lineages.merge(selected, on="ClusterID", how="left")
        abundance_column = f"{selection}_abundance"
        lineages[abundance_column] = lineages[abundance_column].fillna(0.0)
        total = lineages[abundance_column].sum()
        lineages[f"{selection}_relative_abundance"] = lineages[abundance_column] / total if total > 0 else 0.0
        lineages[f"{selection}_detected"] = lineages[abundance_column] > 0
    return lineages


def comparison_slug(left: str, right: str) -> str:
    return f"{left}_vs_{right}"


def top_n_table(lineages: pd.DataFrame, selection: str, top_n: int) -> pd.DataFrame:
    columns = ["ClusterID", f"{selection}_relative_abundance"]
    ranked = lineages.loc[lineages[f"{selection}_detected"], columns].copy()
    ranked = ranked.sort_values([f"{selection}_relative_abundance", "ClusterID"], ascending=[False, True], kind="stable")
    ranked = ranked.head(top_n).reset_index(drop=True)
    ranked["rank"] = ranked.index + 1
    return ranked


def comparison_tables(
    lineages: pd.DataFrame, left: str, right: str, top_n: int
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    left_relative = f"{left}_relative_abundance"
    right_relative = f"{right}_relative_abundance"
    left_detected = f"{left}_detected"
    right_detected = f"{right}_detected"
    pair = lineages[["ClusterID", left_relative, right_relative, left_detected, right_detected]].copy()
    pair["detection_status"] = "neither"
    pair.loc[pair[left_detected] & ~pair[right_detected], "detection_status"] = "left_only"
    pair.loc[~pair[left_detected] & pair[right_detected], "detection_status"] = "right_only"
    pair.loc[pair[left_detected] & pair[right_detected], "detection_status"] = "both"

    shared = pair.loc[pair["detection_status"] == "both"].copy()
    rho: float | None = None
    p_value: float | None = None
    if len(shared) >= 2:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ConstantInputWarning)
            result = spearmanr(shared[left_relative], shared[right_relative])
        if math.isfinite(float(result.statistic)):
            rho = float(result.statistic)
        if math.isfinite(float(result.pvalue)):
            p_value = float(result.pvalue)

    left_top = top_n_table(lineages, left, top_n).rename(columns={"rank": "left_rank", left_relative: "left_relative_abundance"})
    right_top = top_n_table(lineages, right, top_n).rename(columns={"rank": "right_rank", right_relative: "right_relative_abundance"})
    membership = left_top.merge(right_top, on="ClusterID", how="outer")
    membership["in_left_top_n"] = membership["left_rank"].notna()
    membership["in_right_top_n"] = membership["right_rank"].notna()
    membership["in_intersection"] = membership["in_left_top_n"] & membership["in_right_top_n"]
    membership = membership.sort_values(["in_intersection", "ClusterID"], ascending=[False, True], kind="stable")
    overlap_count = int(membership["in_intersection"].sum())
    union_count = int(len(membership))
    metrics = {
        "shared_detection_count": int(len(shared)),
        "spearman_rho": rho,
        "spearman_p_value": p_value,
        "top_n_overlap_count": overlap_count,
        "top_n_jaccard_similarity": (overlap_count / union_count) if union_count else 0.0,
    }
    return pair, membership, metrics


def plot_comparison(pair: pd.DataFrame, left: str, right: str, output_path: Path) -> None:
    shared = pair.loc[pair["detection_status"] == "both"]
    left_relative = f"{left}_relative_abundance"
    right_relative = f"{right}_relative_abundance"
    plt.figure(figsize=(6, 6))
    plt.scatter(shared[left_relative], shared[right_relative], color="#4C78A8")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(f"{SELECTION_LABELS[left]} relative abundance")
    plt.ylabel(f"{SELECTION_LABELS[right]} relative abundance")
    plt.title(f"{SELECTION_LABELS[left]} vs {SELECTION_LABELS[right]}")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def write_tsv(table: pd.DataFrame, path: Path) -> None:
    table.to_csv(path, sep="\t", index=False, na_rep="NA")


def run(mapped_path: Path, out_dir: Path, top_n: int, min_overlap: int) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    table = read_table(mapped_path)
    abundance_column = require_columns(table)
    classified = classify_mapped_rows(table, min_overlap)
    evidence = collapse_quantitative_evidence(classified, abundance_column)
    lineages = lineage_abundance_table(evidence)
    stem = mapped_path.stem

    outputs: dict[str, Path] = {}
    abundance_path = out_dir / f"{stem}_cdr_lineage_abundance.tsv"
    write_tsv(lineages, abundance_path)
    outputs["lineage_abundance"] = abundance_path

    summary_rows: list[dict[str, Any]] = []
    for left, right in COMPARISONS:
        slug = comparison_slug(left, right)
        pair, membership, metrics = comparison_tables(lineages, left, right, top_n)
        pair_path = out_dir / f"{stem}_{slug}_lineages.tsv"
        top_path = out_dir / f"{stem}_{slug}_top_{top_n}_membership.tsv"
        write_tsv(pair, pair_path)
        write_tsv(membership, top_path)
        outputs[f"{slug}_lineages"] = pair_path
        outputs[f"{slug}_top_membership"] = top_path
        scatter_path: Path | None = None
        if metrics["shared_detection_count"] >= 2:
            scatter_path = out_dir / f"{stem}_{slug}_scatter.png"
            plot_comparison(pair, left, right, scatter_path)
            outputs[f"{slug}_scatter"] = scatter_path
        summary_rows.append({
            "comparison": slug,
            "left_selection": left,
            "right_selection": right,
            **metrics,
            "scatterplot_path": str(scatter_path) if scatter_path else "",
        })

    summary_path = out_dir / f"{stem}_cdr_lineage_comparison_summary.tsv"
    write_tsv(pd.DataFrame(summary_rows), summary_path)
    outputs["comparison_summary"] = summary_path
    log_path = out_dir / f"{stem}_compare_cdr_lineage_abundance_log_{timestamp()}.json"
    outputs["run_log"] = log_path
    log = {
        "command": " ".join(sys.argv),
        "mapped_peptides_tsv": str(mapped_path),
        "abundance_column": abundance_column,
        "top_n": top_n,
        "min_cdr_overlap_aa": min_overlap,
        "input_rows": int(len(table)),
        "successful_mapping_rows": int(len(classified)),
        "unique_lineage_peptides": int(len(evidence)),
        "lineage_count": int(len(lineages)),
        "selection_labels": SELECTION_LABELS,
        "comparisons": summary_rows,
        "outputs": {name: str(path) for name, path in outputs.items()},
    }
    log_path.write_text(json.dumps(log, indent=2))
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare CDR-derived lineage abundance estimates from mapped peptides.")
    parser.add_argument("mapped_peptides_tsv", help="Mapped peptide TSV from quantify_map_peptides.py.")
    parser.add_argument("--top-n", type=int, default=10, help="Number of top lineages to compare (default: 10).")
    parser.add_argument("--min-cdr-overlap-aa", type=int, default=3, help="Minimum peptide/CDR overlap in amino acids (default: 3).")
    parser.add_argument("--out-dir", help="Output directory (default: timestamped directory next to the input).")
    args = parser.parse_args()
    if args.top_n < 1:
        parser.error("--top-n must be at least 1.")
    if args.min_cdr_overlap_aa < 1:
        parser.error("--min-cdr-overlap-aa must be at least 1.")

    mapped_path = Path(args.mapped_peptides_tsv).expanduser()
    outputs = run(mapped_path, make_out_dir(mapped_path, args.out_dir), args.top_n, args.min_cdr_overlap_aa)
    for name, path in outputs.items():
        print(f"Wrote {name}: {path}")


if __name__ == "__main__":
    main()
