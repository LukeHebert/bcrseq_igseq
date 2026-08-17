#!/usr/bin/env python3
"""Simulate protease digestion around heavy-chain CDR3 sequences.

The input is an AIRR/IgBLAST-style TSV with ``fwr3_aa``, ``cdr3_aa``,
``fwr4_aa``, and ``nt_seq_count`` columns.  The tool applies each protease
rule to the concatenated FR3-CDR3-FR4 sequence and reports how often a full
CDR3 is retained in a resulting peptide.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
import re
from typing import Iterable

import matplotlib.pyplot as plt
import pandas as pd


REQUIRED_SEQUENCE_COLUMNS = ("fwr3_aa", "cdr3_aa", "fwr4_aa")
REQUIRED_INPUT_COLUMNS = REQUIRED_SEQUENCE_COLUMNS + ("nt_seq_count",)
RULE_COLUMNS = ("short_name", "P4", "P3", "P2", "P1", "P1'", "P2'")
AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWYBXZJUO")


@dataclass(frozen=True)
class ProteaseRule:
    """One accepted cleavage-site context for a protease."""

    positions: dict[str, str]


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""
    default_rules = Path(__file__).with_name("data") / "proteases.tsv"
    parser = argparse.ArgumentParser(
        description=(
            "Simulate digestion of AIRR/IgBLAST heavy-chain FR3-CDR3-FR4 "
            "sequences and summarize CDR3-spanning peptides."
        )
    )
    parser.add_argument("annotation_tsv", type=Path, help="AIRR/IgBLAST annotation TSV.")
    parser.add_argument(
        "--protease-rules",
        type=Path,
        default=default_rules,
        help=f"Protease specificity TSV (default: {default_rules}).",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help=(
            "Directory for results. Default: a timestamped protease_simulation "
            "folder beside the input TSV."
        ),
    )
    return parser.parse_args()


def default_output_directory(annotation_tsv: Path) -> Path:
    """Return a reproducible timestamped output location beside the input."""
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    return annotation_tsv.parent / "protease_simulation" / timestamp


def validate_amino_acid_sequence(value: object, column: str, row_number: int) -> str:
    """Return an upper-case amino-acid sequence or raise a helpful error."""
    sequence = str(value).strip().upper()
    # Some historical tables mark a truncated FR4 with one trailing underscore.
    if column == "fwr4_aa":
        sequence = sequence.rstrip("_")
    if not sequence:
        raise ValueError(f"Row {row_number}: {column} is empty.")
    invalid = sorted(set(sequence) - AMINO_ACIDS)
    if invalid:
        raise ValueError(
            f"Row {row_number}: {column} contains unsupported characters: {''.join(invalid)}."
        )
    return sequence


def read_annotation_table(path: Path) -> pd.DataFrame:
    """Load and validate the AIRR/IgBLAST sequence and count columns."""
    try:
        table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False).fillna("")
    except OSError as error:
        raise ValueError(f"Could not read annotation TSV {path}: {error}") from error

    missing = [column for column in REQUIRED_INPUT_COLUMNS if column not in table.columns]
    if missing:
        raise ValueError(
            "Annotation TSV is missing required AIRR/IgBLAST columns: "
            + ", ".join(missing)
            + "."
        )

    if table.empty:
        raise ValueError("Annotation TSV contains no records.")

    validated = table.copy()
    for row_index, row in validated.iterrows():
        row_number = row_index + 2
        for column in REQUIRED_SEQUENCE_COLUMNS:
            validated.at[row_index, column] = validate_amino_acid_sequence(
                row[column], column, row_number
            )

    counts = pd.to_numeric(validated["nt_seq_count"], errors="coerce")
    invalid_counts = counts.isna() | (counts <= 0) | (counts % 1 != 0)
    if invalid_counts.any():
        bad_rows = ", ".join(str(index + 2) for index in validated.index[invalid_counts][:5])
        raise ValueError(
            "nt_seq_count must be a positive integer for every row. "
            f"Invalid value(s) found at TSV row(s): {bad_rows}."
        )
    validated["nt_seq_count"] = counts.astype(int)
    return validated


def validate_requirement(requirement: str, column: str, row_number: int) -> str:
    """Validate the legacy semicolon-separated specificity syntax."""
    normalized = requirement.strip().upper()
    if not normalized:
        return ""

    values = normalized.split(";")
    negative = all(value.startswith("0") for value in values)
    positive = all(not value.startswith("0") for value in values)
    if not (negative or positive) or any(not value or value == "0" for value in values):
        raise ValueError(
            f"Protease-rule row {row_number}, {column}: use either AA;AA or 0AA;0AA."
        )
    residue_values = [value[1:] if negative else value for value in values]
    if any(set(value) - AMINO_ACIDS for value in residue_values):
        raise ValueError(
            f"Protease-rule row {row_number}, {column} contains unsupported amino-acid code(s)."
        )
    return normalized


def read_protease_rules(path: Path) -> dict[str, list[ProteaseRule]]:
    """Load a named protease-rule TSV into rules grouped by protease name."""
    try:
        table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False).fillna("")
    except OSError as error:
        raise ValueError(f"Could not read protease-rule TSV {path}: {error}") from error

    missing = [column for column in RULE_COLUMNS if column not in table.columns]
    if missing:
        raise ValueError("Protease-rule TSV is missing columns: " + ", ".join(missing) + ".")
    if table.empty:
        raise ValueError("Protease-rule TSV contains no rules.")

    proteases: dict[str, list[ProteaseRule]] = {}
    position_columns = RULE_COLUMNS[1:]
    for row_index, row in table.iterrows():
        row_number = row_index + 2
        name = row["short_name"].strip()
        if not name:
            raise ValueError(f"Protease-rule row {row_number} has an empty short_name.")
        positions = {
            column: validate_requirement(row[column], column, row_number)
            for column in position_columns
        }
        proteases.setdefault(name, []).append(ProteaseRule(positions))
    return proteases


def requirement_matches(requirement: str, observed: str) -> bool:
    """Check one legacy rule field against its observed sequence context."""
    if not requirement:
        return True
    options = requirement.split(";")
    if options[0].startswith("0"):
        return observed not in {option[1:] for option in options}
    return observed in set(options)


def context_at_boundary(sequence: str, boundary: int) -> dict[str, str]:
    """Return Schechter-Berger positions for a cut between two residues."""
    def left(offset: int) -> str:
        index = boundary - offset
        return sequence[index] if index >= 0 else ""

    def right(offset: int) -> str:
        index = boundary + offset - 1
        return sequence[index] if index < len(sequence) else ""

    return {
        "P4": left(4),
        "P3": left(3),
        "P2": left(2),
        "P1": left(1),
        "P1'": right(1),
        "P2'": right(2),
    }


def find_cleavage_boundaries(sequence: str, rules: Iterable[ProteaseRule]) -> list[int]:
    """Return unique boundaries where at least one rule permits cleavage."""
    boundaries: list[int] = []
    for boundary in range(len(sequence) + 1):
        context = context_at_boundary(sequence, boundary)
        if any(
            all(requirement_matches(rule.positions[position], context[position]) for position in rule.positions)
            for rule in rules
        ):
            boundaries.append(boundary)
    return boundaries


def peptides_from_boundaries(sequence: str, cleavage_boundaries: Iterable[int]) -> list[str]:
    """Split a sequence at cleavage boundaries without duplicating residues."""
    boundaries = sorted({0, len(sequence), *cleavage_boundaries})
    return [sequence[start:end] for start, end in zip(boundaries, boundaries[1:]) if start != end]


def nearest_distance(boundaries: list[int], target: int, upstream: bool) -> int | None:
    """Find the nearest cleavage boundary upstream or downstream of target."""
    candidates = [boundary for boundary in boundaries if boundary <= target]
    if not upstream:
        candidates = [boundary for boundary in boundaries if boundary >= target]
    if not candidates:
        return None
    boundary = max(candidates) if upstream else min(candidates)
    return target - boundary if upstream else boundary - target


def digest_record(row: pd.Series, protease: str, rules: list[ProteaseRule]) -> dict[str, object]:
    """Digest one FR3-CDR3-FR4 sequence and return CDR3-focused measurements."""
    fr3, cdr3, fr4 = (row[column] for column in REQUIRED_SEQUENCE_COLUMNS)
    sequence = fr3 + cdr3 + fr4
    cdr3_start = len(fr3)
    cdr3_end = cdr3_start + len(cdr3)
    boundaries = find_cleavage_boundaries(sequence, rules)
    peptides = peptides_from_boundaries(sequence, boundaries)
    cdr3_peptides = [peptide for peptide in peptides if cdr3 in peptide]
    cdr3_peptide = cdr3_peptides[0] if cdr3_peptides else None
    internal_cuts = [boundary for boundary in boundaries if cdr3_start < boundary < cdr3_end]

    result: dict[str, object] = {
        "input_row": row.name + 2,
        "protease": protease,
        "nt_seq_count": row["nt_seq_count"],
        "fwr3_aa": fr3,
        "cdr3_aa": cdr3,
        "fwr4_aa": fr4,
        "cleavage_boundaries": ";".join(str(boundary) for boundary in boundaries),
        "peptides": ";".join(peptides),
        "cuts_inside_cdr3": bool(internal_cuts),
        "cdr3_spanning_peptide": cdr3_peptide or "",
        "cdr3_spanning_peptide_length": len(cdr3_peptide) if cdr3_peptide else pd.NA,
        "upstream_cut_distance_aa": nearest_distance(boundaries, cdr3_start, upstream=True),
        "downstream_cut_distance_aa": nearest_distance(boundaries, cdr3_end, upstream=False),
    }
    if "bcrseq_id" in row.index:
        result["bcrseq_id"] = row["bcrseq_id"]
    return result


def weighted_mean(values: pd.Series, weights: pd.Series) -> float | None:
    """Calculate a weighted mean while excluding missing measurements."""
    present = values.notna()
    if not present.any():
        return None
    return float((values[present] * weights[present]).sum() / weights[present].sum())


def create_summary(details: pd.DataFrame) -> pd.DataFrame:
    """Create one weighted CDR3-digestion summary row per protease."""
    records: list[dict[str, object]] = []
    for protease, group in details.groupby("protease", sort=True):
        weights = group["nt_seq_count"]
        total = int(weights.sum())
        cdr3_spanning = group["cdr3_spanning_peptide_length"].notna()
        internal_cuts = group["cuts_inside_cdr3"]
        records.append(
            {
                "protease": protease,
                "unique_sequence_count": len(group),
                "weighted_sequence_count": total,
                "weighted_cdr3_spanning_count": int(weights[cdr3_spanning].sum()),
                "weighted_cdr3_spanning_percent": 100 * weights[cdr3_spanning].sum() / total,
                "weighted_internal_cdr3_cut_count": int(weights[internal_cuts].sum()),
                "weighted_internal_cdr3_cut_percent": 100 * weights[internal_cuts].sum() / total,
                "weighted_mean_cdr3_spanning_peptide_length": weighted_mean(
                    group["cdr3_spanning_peptide_length"], weights
                ),
                "weighted_mean_upstream_cut_distance_aa": weighted_mean(
                    group["upstream_cut_distance_aa"], weights
                ),
                "weighted_mean_downstream_cut_distance_aa": weighted_mean(
                    group["downstream_cut_distance_aa"], weights
                ),
            }
        )
    return pd.DataFrame(records)


def safe_filename(text: str) -> str:
    """Make a protease name safe for use in a filename."""
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text).strip("_")


def plot_weighted_histogram(
    details: pd.DataFrame, column: str, xlabel: str, title: str, output_path: Path
) -> None:
    """Write an integer-valued histogram weighted by transcript count."""
    plotted = details.loc[details[column].notna(), [column, "nt_seq_count"]]
    figure, axis = plt.subplots(figsize=(8, 5))
    if plotted.empty:
        axis.text(0.5, 0.5, "No qualifying peptides", ha="center", va="center", transform=axis.transAxes)
        axis.set_xlabel(xlabel)
        axis.set_ylabel("Weighted BCRseq transcript count")
        axis.set_title(title)
        figure.tight_layout()
        figure.savefig(output_path, dpi=300)
        plt.close(figure)
        return
    values = plotted[column].astype(int)
    maximum = int(values.max())
    axis.hist(
        values,
        bins=range(0, maximum + 2),
        weights=plotted["nt_seq_count"],
        align="left",
        color="#4c78a8",
        edgecolor="white",
    )
    axis.set_xlabel(xlabel)
    axis.set_ylabel("Weighted BCRseq transcript count")
    axis.set_title(title)
    axis.set_xticks(range(0, maximum + 1))
    figure.tight_layout()
    figure.savefig(output_path, dpi=300)
    plt.close(figure)


def write_plots(details: pd.DataFrame, output_directory: Path) -> None:
    """Write the three CDR3-focused plots for every protease."""
    plots = (
        ("upstream_cut_distance_aa", "Upstream cut distance from CDR3 (aa)", "upstream_cut_distance"),
        ("downstream_cut_distance_aa", "Downstream cut distance from CDR3 (aa)", "downstream_cut_distance"),
        (
            "cdr3_spanning_peptide_length",
            "CDR3-spanning peptide length (aa)",
            "cdr3_spanning_peptide_length",
        ),
    )
    for protease, group in details.groupby("protease", sort=True):
        for column, xlabel, suffix in plots:
            plot_weighted_histogram(
                group,
                column,
                xlabel,
                f"{protease}: {xlabel}",
                output_directory / f"{safe_filename(protease)}_{suffix}.png",
            )


def main() -> None:
    """Run the complete in-silico digestion analysis."""
    args = parse_arguments()
    annotations = read_annotation_table(args.annotation_tsv)
    proteases = read_protease_rules(args.protease_rules)
    output_directory = args.out_dir or default_output_directory(args.annotation_tsv)
    output_directory.mkdir(parents=True, exist_ok=False)

    details = pd.DataFrame(
        digest_record(row, name, rules)
        for name, rules in proteases.items()
        for _, row in annotations.iterrows()
    )
    details.to_csv(output_directory / "protease_digestion_details.tsv", sep="\t", index=False)
    create_summary(details).to_csv(
        output_directory / "protease_digestion_summary.tsv", sep="\t", index=False
    )
    write_plots(details, output_directory)
    print(f"Wrote protease-digestion results to: {output_directory}")


if __name__ == "__main__":
    try:
        main()
    except ValueError as error:
        raise SystemExit(f"Error: {error}") from error
