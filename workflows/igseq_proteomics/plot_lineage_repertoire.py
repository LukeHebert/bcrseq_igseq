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
from Bio.Align import PairwiseAligner


REGIONS = ["fwr1_aa", "cdr1_aa", "fwr2_aa", "cdr2_aa", "fwr3_aa", "cdr3_aa", "fwr4_aa"]
AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")
MIN_LOGO_FIGURE_WIDTH = 18
LOGO_WIDTH_PER_RESIDUE = 0.22
LOGO_FIGURE_HEIGHT = 4.5
LOGO_COVERAGE_HEIGHT_RATIO = 1.2
LOGO_SEQUENCE_HEIGHT_RATIO = 1.7
METRIC_CHOICES = ["abundance_total", "abundance_average", "psm_total", "psm_average"]


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


def build_raw_region_map(row: pd.Series) -> dict[str, str]:
    region_map: dict[str, str] = {}
    for region in REGIONS:
        value = row.get(region, "")
        region_map[region] = str(value).upper() if pd.notna(value) else ""
    return region_map


def metric_config(metric: str) -> dict[str, str]:
    configs = {
        "abundance_total": {
            "primary_column": "peptide_abundance_total",
            "fallback_column": "peptide_abundance",
            "lineage_y_label": "Relative abundance",
            "coverage_y_label": "Peptide abundance",
            "detail_title_label": "peptide coverage",
        },
        "abundance_average": {
            "primary_column": "avg_abundance_per_injection",
            "fallback_column": "",
            "lineage_y_label": "Relative average abundance",
            "coverage_y_label": "Average peptide abundance",
            "detail_title_label": "average peptide coverage",
        },
        "psm_total": {
            "primary_column": "psm_count_total",
            "fallback_column": "psm_count",
            "lineage_y_label": "Relative total PSM count",
            "coverage_y_label": "Total PSM count",
            "detail_title_label": "total PSM coverage",
        },
        "psm_average": {
            "primary_column": "avg_psm_per_injection",
            "fallback_column": "",
            "lineage_y_label": "Relative average PSM count",
            "coverage_y_label": "Average PSM count",
            "detail_title_label": "average PSM coverage",
        },
    }
    return configs[metric]


def resolve_metric_series(df: pd.DataFrame, metric: str) -> pd.Series:
    config = metric_config(metric)
    primary = config["primary_column"]
    fallback = config["fallback_column"]
    if primary in df.columns:
        return to_number(df[primary]).fillna(0.0)
    if fallback and fallback in df.columns:
        return to_number(df[fallback]).fillna(0.0)
    raise ValueError(
        f"Mapped peptide table is missing the column needed for metric '{metric}'. "
        f"Expected '{primary}'"
        + (f" or '{fallback}'." if fallback else ".")
    )


def plot_lineage_abundance(
    lineage: pd.DataFrame,
    out_path: Path,
    threshold: float,
    metric: str,
) -> None:
    config = metric_config(metric)
    plot_df = lineage[lineage["relative_abundance"] >= threshold].copy()
    plot_df = plot_df.sort_values("relative_abundance", ascending=False)
    plt.figure(figsize=(max(8, len(plot_df) * 0.35), 5))
    plt.bar(plot_df["ClusterID"].astype(str), plot_df["relative_abundance"], color="#4C78A8")
    plt.ylabel(config["lineage_y_label"])
    plt.xlabel("ClusterID")
    plt.title(f"Lineage abundance by {metric} ({config['lineage_y_label'].lower()} >= {threshold:.3g})")
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
        row["aligned_position"] = pos + 1
        rows.append(row)
    return pd.DataFrame(rows)


def require_logomaker():
    try:
        import logomaker  # type: ignore
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "plot_lineage_repertoire.py requires the 'logomaker' package for "
            "letter-based sequence logos. Install project dependencies from "
            "requirements.txt or install logomaker into your active Python environment."
        ) from exc
    return logomaker


def make_pairwise_aligner() -> PairwiseAligner:
    aligner = PairwiseAligner(mode="global")
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -2.0
    aligner.extend_gap_score = -0.5
    aligner.target_end_gap_score = 0.0
    aligner.query_end_gap_score = 0.0
    return aligner


def aligned_strings_from_alignment(target: str, query: str, alignment: Any) -> tuple[str, str]:
    coords = alignment.coordinates
    target_parts: list[str] = []
    query_parts: list[str] = []
    target_prev = int(coords[0, 0])
    query_prev = int(coords[1, 0])

    for idx in range(1, coords.shape[1]):
        target_curr = int(coords[0, idx])
        query_curr = int(coords[1, idx])
        target_step = target_curr - target_prev
        query_step = query_curr - query_prev

        if target_step > 0 and query_step > 0:
            target_parts.append(target[target_prev:target_curr])
            query_parts.append(query[query_prev:query_curr])
        elif target_step > 0 and query_step == 0:
            target_parts.append(target[target_prev:target_curr])
            query_parts.append("-" * target_step)
        elif target_step == 0 and query_step > 0:
            target_parts.append("-" * query_step)
            query_parts.append(query[query_prev:query_curr])

        target_prev = target_curr
        query_prev = query_curr

    return "".join(target_parts), "".join(query_parts)


def profile_consensus(profile_aligned: dict[str, str], profile_order: list[str]) -> str:
    if not profile_order:
        return ""
    consensus_chars: list[str] = []
    width = len(profile_aligned[profile_order[0]])
    for col in range(width):
        counts = Counter(
            profile_aligned[member_id][col]
            for member_id in profile_order
            if profile_aligned[member_id][col] != "-"
        )
        consensus_chars.append(counts.most_common(1)[0][0] if counts else "-")
    return "".join(consensus_chars)


def raw_to_aligned_map(aligned_seq: str) -> list[int]:
    mapping: list[int] = []
    for aligned_idx, char in enumerate(aligned_seq):
        if char != "-":
            mapping.append(aligned_idx)
    return mapping


def align_region_progressive(region_name: str, region_sequences: dict[str, str]) -> dict[str, Any]:
    member_ids = list(region_sequences.keys())
    nonempty = [(member_id, region_sequences[member_id]) for member_id in member_ids if region_sequences[member_id]]
    if not nonempty:
        return {
            "region": region_name,
            "aligned_sequences": {member_id: "" for member_id in member_ids},
            "raw_to_aligned": {member_id: [] for member_id in member_ids},
            "aligned_width": 0,
            "seed_member_id": "",
            "columns_added_before": 0,
            "columns_added_after": 0,
            "raw_lengths": {member_id: 0 for member_id in member_ids},
        }

    seed_member_id, seed_seq = max(nonempty, key=lambda item: (len(item[1]), item[0]))
    profile_order = [seed_member_id]
    profile_aligned = {seed_member_id: seed_seq}
    columns_added_before = 0
    columns_added_after = 0
    aligner = make_pairwise_aligner()

    for member_id in member_ids:
        if member_id == seed_member_id:
            continue
        seq = region_sequences[member_id]
        if not seq:
            continue

        consensus = profile_consensus(profile_aligned, profile_order)
        alignment = aligner.align(consensus, seq)[0]
        aligned_consensus, aligned_seq = aligned_strings_from_alignment(consensus, seq, alignment)

        expanded: dict[str, list[str]] = {existing_id: [] for existing_id in profile_order}
        expanded[member_id] = []
        old_col = 0
        old_width = len(consensus)

        for cons_char, seq_char in zip(aligned_consensus, aligned_seq):
            if cons_char == "-":
                for existing_id in profile_order:
                    expanded[existing_id].append("-")
                expanded[member_id].append(seq_char)
                if old_col == 0:
                    columns_added_before += 1
                elif old_col == old_width:
                    columns_added_after += 1
            else:
                for existing_id in profile_order:
                    expanded[existing_id].append(profile_aligned[existing_id][old_col])
                expanded[member_id].append(seq_char)
                old_col += 1

        profile_order.append(member_id)
        profile_aligned = {existing_id: "".join(chars) for existing_id, chars in expanded.items()}

    aligned_sequences: dict[str, str] = {}
    raw_maps: dict[str, list[int]] = {}
    seed_width = len(profile_aligned[seed_member_id])
    for member_id in member_ids:
        aligned_seq = profile_aligned.get(member_id, "-" * seed_width)
        if not region_sequences[member_id]:
            aligned_seq = "-" * seed_width
        aligned_sequences[member_id] = aligned_seq
        raw_maps[member_id] = raw_to_aligned_map(aligned_seq)

    return {
        "region": region_name,
        "aligned_sequences": aligned_sequences,
        "raw_to_aligned": raw_maps,
        "aligned_width": seed_width,
        "seed_member_id": seed_member_id,
        "columns_added_before": columns_added_before,
        "columns_added_after": columns_added_after,
        "raw_lengths": {member_id: len(region_sequences[member_id]) for member_id in member_ids},
    }


def build_lineage_alignment(bcr_rows: pd.DataFrame) -> dict[str, Any]:
    member_ids = bcr_rows["bcrseq_id"].astype(str).tolist()
    raw_regions = {member_id: build_raw_region_map(row) for member_id, (_, row) in zip(member_ids, bcr_rows.iterrows())}
    region_alignments = []
    aligned_spans: list[dict[str, Any]] = []
    running_offset = 0

    for region in REGIONS:
        region_sequences = {member_id: raw_regions[member_id][region] for member_id in member_ids}
        region_alignment = align_region_progressive(region, region_sequences)
        region_alignments.append(region_alignment)
        width = int(region_alignment["aligned_width"])
        aligned_spans.append(
            {
                "region": region,
                "start": running_offset,
                "end": running_offset + width,
                "width": width,
            }
        )
        running_offset += width

    members: dict[str, dict[str, Any]] = {}
    for member_id in member_ids:
        aligned_parts: list[str] = []
        raw_parts: list[str] = []
        full_raw_to_aligned: list[int] = []
        full_region_records: list[dict[str, Any]] = []
        running_offset = 0

        for region_alignment in region_alignments:
            region = str(region_alignment["region"])
            raw_seq = raw_regions[member_id][region]
            aligned_seq = region_alignment["aligned_sequences"][member_id]
            local_map = region_alignment["raw_to_aligned"][member_id]
            full_raw_to_aligned.extend([running_offset + idx for idx in local_map])
            full_region_records.append(
                {
                    "region": region,
                    "raw_sequence": raw_seq,
                    "aligned_sequence": aligned_seq,
                    "aligned_width": int(region_alignment["aligned_width"]),
                }
            )
            aligned_parts.append(aligned_seq)
            raw_parts.append(raw_seq)
            running_offset += int(region_alignment["aligned_width"])

        members[member_id] = {
            "bcrseq_id": member_id,
            "raw_full_seq": "".join(raw_parts),
            "aligned_full_seq": "".join(aligned_parts),
            "raw_to_aligned": full_raw_to_aligned,
            "region_records": full_region_records,
        }

    return {
        "members": members,
        "region_alignments": region_alignments,
        "aligned_region_spans": aligned_spans,
        "aligned_width": sum(span["width"] for span in aligned_spans),
        "member_order": member_ids,
    }


def project_group_to_aligned_columns(peptide_group: pd.DataFrame, members: dict[str, dict[str, Any]]) -> tuple[set[int], int]:
    covered_columns: set[int] = set()
    rows = peptide_group.drop_duplicates(["bcrseq_id", "match_start", "match_end"]).copy()
    suffix_skipped = 0

    for _, row in rows.iterrows():
        member_id = str(row.get("bcrseq_id", ""))
        if member_id not in members:
            continue
        member = members[member_id]
        raw_to_aligned = member["raw_to_aligned"]
        try:
            start = int(float(row["match_start"]))
            end = int(float(row["match_end"]))
        except Exception:
            continue

        clipped_start = max(0, start)
        clipped_end = min(len(raw_to_aligned), end)
        if clipped_start >= clipped_end:
            suffix_skipped += 1
            continue
        covered_columns.update(raw_to_aligned[clipped_start:clipped_end])

    return covered_columns, suffix_skipped


def plot_lineage_detail(
    cluster_id: str,
    bcr_rows: pd.DataFrame,
    peptide_rows: pd.DataFrame,
    out_path: Path,
    metric: str,
) -> dict[str, Any]:
    if bcr_rows.empty:
        return {"cluster_id": cluster_id, "warning": "No BCR sequences available."}

    config = metric_config(metric)
    logomaker = require_logomaker()
    alignment = build_lineage_alignment(bcr_rows)
    members = alignment["members"]
    member_order = alignment["member_order"]
    aligned_width = int(alignment["aligned_width"])
    if aligned_width == 0:
        return {"cluster_id": cluster_id, "warning": "No aligned region content available."}

    aligned_sequences = [members[member_id]["aligned_full_seq"] for member_id in member_order]
    coverage = [0.0] * aligned_width
    suffix_skip_count = 0
    union_projection_count = 0

    peptide_rows = peptide_rows[peptide_rows["match_found"].astype(str).str.lower().isin(["true", "1"])].copy()
    for _, peptide_group in peptide_rows.groupby("peptide_sequence", dropna=False):
        covered_columns, group_suffix_skipped = project_group_to_aligned_columns(peptide_group, members)
        suffix_skip_count += group_suffix_skipped
        if len(peptide_group.drop_duplicates(["bcrseq_id", "match_start", "match_end"])) > 1:
            union_projection_count += 1
        if not covered_columns:
            continue
        abundance = float(peptide_group.iloc[0]["metric_value_num"])
        for aligned_idx in covered_columns:
            coverage[aligned_idx] += abundance

    freq = aa_frequency_matrix(aligned_sequences)
    logo_matrix = freq.set_index("aligned_position")[AA_ORDER]

    coverage_df = pd.DataFrame(
        {
            "aligned_position": list(range(1, aligned_width + 1)),
            "coverage_value": coverage,
        }
    )
    coverage_df["metric"] = metric
    coverage_df["region"] = ""
    coverage_df["column_has_residue"] = [
        any(seq[pos] in AA_ORDER for seq in aligned_sequences) for pos in range(aligned_width)
    ]

    for span in alignment["aligned_region_spans"]:
        if span["end"] <= span["start"]:
            continue
        coverage_df.loc[
            (coverage_df["aligned_position"] > span["start"])
            & (coverage_df["aligned_position"] <= span["end"]),
            "region",
        ] = span["region"]

    matrix_path = out_path.with_name(out_path.stem + "_logo_matrix.tsv")
    coverage_path = out_path.with_name(out_path.stem + "_coverage.tsv")
    sequences_path = out_path.with_name(out_path.stem + "_sequences.tsv")
    freq.to_csv(matrix_path, sep="\t", index=False)
    coverage_df.to_csv(coverage_path, sep="\t", index=False)

    sequence_rows = []
    for member_id in member_order:
        member = members[member_id]
        row: dict[str, Any] = {
            "ClusterID": cluster_id,
            "bcrseq_id": member_id,
            "raw_full_seq": member["raw_full_seq"],
            "aligned_full_seq": member["aligned_full_seq"],
            "raw_sequence_length": len(member["raw_full_seq"]),
            "aligned_sequence_length": len(member["aligned_full_seq"]),
        }
        for record in member["region_records"]:
            row[record["region"]] = record["raw_sequence"]
            row[f"aligned_{record['region']}"] = record["aligned_sequence"]
        sequence_rows.append(row)
    pd.DataFrame(sequence_rows).to_csv(sequences_path, sep="\t", index=False)

    figure_width = max(MIN_LOGO_FIGURE_WIDTH, aligned_width * LOGO_WIDTH_PER_RESIDUE)
    fig, axes = plt.subplots(
        2,
        1,
        figsize=(figure_width, LOGO_FIGURE_HEIGHT),
        sharex=True,
        gridspec_kw={"height_ratios": [LOGO_COVERAGE_HEIGHT_RATIO, LOGO_SEQUENCE_HEIGHT_RATIO]},
    )
    x_positions = list(range(1, aligned_width + 1))
    axes[0].bar(x_positions, coverage, width=1.0, color="#F58518")
    axes[0].set_ylabel(config["coverage_y_label"])
    axes[0].set_title(f"Cluster {cluster_id}: {config['detail_title_label']} and aligned BCR sequence logo")

    logomaker.Logo(logo_matrix, ax=axes[1])
    axes[1].set_ylabel("AA frequency")
    axes[1].set_xlabel("Aligned residue position")
    tick_step = max(1, aligned_width // 20)
    tick_positions = list(range(1, aligned_width + 1, tick_step))
    axes[1].set_xticks(tick_positions)
    axes[1].set_xticklabels([str(pos) for pos in tick_positions], rotation=0)

    ymax = max(coverage) * 1.05 if coverage and max(coverage) > 0 else 1
    for span in alignment["aligned_region_spans"]:
        if span["end"] <= span["start"]:
            continue
        axes[0].axvspan(span["start"] + 0.5, span["end"] + 0.5, alpha=0.08, color="#54A24B")
        mid = (span["start"] + span["end"] + 1) / 2
        axes[0].text(mid, ymax, span["region"].replace("_aa", "").upper(), ha="center", va="top", fontsize=8)
    axes[0].set_ylim(0, ymax * 1.15)
    axes[0].set_xlim(0.5, aligned_width + 0.5)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close(fig)

    region_alignment_notes = []
    for region_alignment in alignment["region_alignments"]:
        region_alignment_notes.append(
            {
                "region": region_alignment["region"],
                "aligned_width": int(region_alignment["aligned_width"]),
                "seed_member_id": region_alignment["seed_member_id"],
                "columns_added_before": int(region_alignment["columns_added_before"]),
                "columns_added_after": int(region_alignment["columns_added_after"]),
                "raw_lengths": region_alignment["raw_lengths"],
            }
        )

    return {
        "cluster_id": cluster_id,
        "metric": metric,
        "bcr_sequence_count": int(len(member_order)),
        "aligned_lineage_width": aligned_width,
        "suffix_skipped_peptide_rows": suffix_skip_count,
        "union_projection_count": union_projection_count,
        "region_alignment": region_alignment_notes,
        "logo_matrix_tsv": str(matrix_path),
        "coverage_tsv": str(coverage_path),
        "sequences_tsv": str(sequences_path),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot quantified lineage abundance and per-lineage coverage summaries."
    )
    parser.add_argument("mapped_peptides", help="Mapped peptide TSV from quantify_map_peptides.py.")
    parser.add_argument("bcrseq_tsv", help="BCRseq annotation TSV.")
    parser.add_argument("--min-cdr3-overlap-aa", type=int, default=0)
    parser.add_argument("--relative-abundance-threshold", type=float, default=0.01)
    parser.add_argument(
        "--metric",
        choices=METRIC_CHOICES,
        default="abundance_total",
        help="Peptide-level metric to use for lineage abundance and coverage plotting.",
    )
    parser.add_argument("--out-dir", help="Output directory. Default: timestamped directory next to mapped peptides.")
    args = parser.parse_args()

    mapped_path = Path(args.mapped_peptides).expanduser()
    bcr_path = Path(args.bcrseq_tsv).expanduser()
    out_dir = make_out_dir(mapped_path, args.out_dir)
    stem = mapped_path.stem

    mapped = read_table(mapped_path)
    bcr = read_table(bcr_path)
    mapped["metric_value_num"] = resolve_metric_series(mapped, args.metric)
    mapped["cdr3_cover_aa_num"] = to_number(mapped["cdr3_cover_aa"]).fillna(0)

    retained = mapped[
        (mapped["match_found"].astype(str).str.lower().isin(["true", "1"]))
        & (mapped["cdr3_cover_aa_num"] >= args.min_cdr3_overlap_aa)
    ].copy()

    abundance_units = retained.drop_duplicates(["ClusterID", "peptide_sequence"]).copy()
    lineage = (
        abundance_units.groupby("ClusterID", dropna=False)["metric_value_num"]
        .sum()
        .reset_index(name="lineage_metric_value")
    )
    lineage["metric"] = args.metric
    total = lineage["lineage_metric_value"].sum()
    lineage["relative_abundance"] = lineage["lineage_metric_value"] / total if total else 0.0
    lineage = lineage.sort_values("relative_abundance", ascending=False)

    lineage_tsv = out_dir / f"{stem}_lineage_abundance.tsv"
    lineage_png = out_dir / f"{stem}_lineage_abundance.png"
    lineage_plot_tsv = out_dir / f"{stem}_lineage_abundance_plotted.tsv"
    lineage.to_csv(lineage_tsv, sep="\t", index=False)
    lineage[lineage["relative_abundance"] >= args.relative_abundance_threshold].to_csv(
        lineage_plot_tsv, sep="\t", index=False
    )
    plot_lineage_abundance(lineage, lineage_png, args.relative_abundance_threshold, args.metric)

    detail_dir = out_dir / "lineage_details"
    detail_dir.mkdir(exist_ok=True)
    selected = lineage[lineage["relative_abundance"] >= args.relative_abundance_threshold]["ClusterID"].astype(str).tolist()
    detail_logs = []
    for cluster_id in selected:
        bcr_rows = bcr[bcr["ClusterID"].astype(str) == cluster_id].copy()
        pep_rows = retained[retained["ClusterID"].astype(str) == cluster_id].copy()
        safe_id = re.sub(r"[^A-Za-z0-9_.-]+", "_", cluster_id)
        detail_logs.append(
            plot_lineage_detail(
                cluster_id,
                bcr_rows,
                pep_rows,
                detail_dir / f"cluster_{safe_id}_coverage_logo.png",
                args.metric,
            )
        )

    log = {
        "command": " ".join(sys.argv),
        "mapped_peptides": str(mapped_path),
        "bcrseq_tsv": str(bcr_path),
        "metric": args.metric,
        "min_cdr3_overlap_aa": args.min_cdr3_overlap_aa,
        "relative_abundance_threshold": args.relative_abundance_threshold,
        "input_mapped_rows": int(len(mapped)),
        "retained_rows": int(len(retained)),
        "retained_metric_value": float(abundance_units["metric_value_num"].sum()),
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
