#!/usr/bin/env python3
"""
Filter Proteome Discoverer PSM exports into heavy-chain, single-lineage PSMs.

This script keeps quantification at the scan/PSM level while removing rows that
cannot be safely assigned to one heavy-chain BCR lineage.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from datetime import datetime
from pathlib import Path
from typing import Iterable

import pandas as pd


REQUIRED_COLUMNS = [
    "Annotated Sequence",
    "Protein Accessions",
    "Spectrum File",
    "First Scan",
    "Percolator PEP",
    "XCorr",
    "Rank",
    "Search Engine Rank",
    "Precursor Abundance",
]


def timestamp() -> str:
    return datetime.now().strftime("%Y-%m-%d_%H-%M-%S")


def make_out_dir(input_path: Path, out_dir: str | None) -> Path:
    if out_dir:
        path = Path(out_dir).expanduser()
    else:
        path = input_path.parent / f"filter_psms_{timestamp()}"
    path.mkdir(parents=True, exist_ok=True)
    return path


def read_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    df.columns = df.columns.str.strip('"')
    return df


def require_columns(df: pd.DataFrame, columns: Iterable[str]) -> None:
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise ValueError("Missing required PSM column(s): " + ", ".join(missing))


def to_number(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def extract_peptide(annotated_sequence: object) -> str:
    text = str(annotated_sequence)
    parts = text.split(".")
    peptide = parts[1] if len(parts) >= 3 else text
    return "".join(re.findall(r"[A-Za-z]", peptide)).upper()


def il_normalize(seq: object) -> str:
    return re.sub(r"[IL]", "X", str(seq).upper())


def split_accessions(accessions: object) -> list[str]:
    if pd.isna(accessions):
        return []
    return [a.strip().strip('"') for a in str(accessions).split(";") if a.strip().strip('"')]


def is_contaminant_row(row: pd.Series) -> bool:
    value = str(row.get("Contaminant", "")).strip().lower()
    if value in {"true", "t", "1", "yes", "y"}:
        return True
    return any(a.startswith("Cont_") for a in row["accession_list"])


def parse_heavy_cluster_id(accession: str) -> str | None:
    """Return ClusterID for heavy-chain BCRseq accessions, otherwise None."""
    if not accession.startswith("BCRseq_"):
        return None
    parts = accession.split("_")
    if len(parts) < 4:
        return None
    candidate = parts[2]
    if candidate.lower() in {"heavy", "igh", "kappa", "lambda", "igk", "igl"}:
        return None
    return candidate if re.fullmatch(r"\d+", candidate) else None


def is_bcrseq(accession: str) -> bool:
    return accession.startswith("BCRseq_")


def is_light_accession(accession: str) -> bool:
    if not accession.startswith("BCRseq_"):
        return False
    parts = accession.split("_")
    return len(parts) >= 3 and parts[2].lower() in {"kappa", "lambda", "igk", "igl"}


def accession_class(accessions: list[str], contaminant: bool) -> str:
    if contaminant:
        return "contaminant"
    if not accessions or not any(is_bcrseq(a) for a in accessions):
        return "non_bcrseq"

    all_bcr = all(is_bcrseq(a) for a in accessions)
    heavy_clusters = [parse_heavy_cluster_id(a) for a in accessions]
    heavy = [c for c in heavy_clusters if c is not None]
    lights = [a for a in accessions if is_light_accession(a)]

    if all_bcr and len(lights) == len(accessions):
        return "light_chain"
    if all_bcr and len(heavy) == len(accessions):
        return "heavy"
    return "mixed_heavy_nonheavy"


def deduplicate_scans(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    work = df.copy()
    work["_pep_num"] = to_number(work["Percolator PEP"])
    work["_xcorr_num"] = to_number(work["XCorr"])
    work["_rank_num"] = to_number(work["Rank"])
    work["_se_rank_num"] = to_number(work["Search Engine Rank"])
    work["_input_order"] = range(len(work))

    sorted_df = work.sort_values(
        by=[
            "Spectrum File",
            "First Scan",
            "_pep_num",
            "_xcorr_num",
            "_rank_num",
            "_se_rank_num",
            "_input_order",
        ],
        ascending=[True, True, True, False, True, True, True],
        na_position="last",
    )
    keep = sorted_df.drop_duplicates(["Spectrum File", "First Scan"], keep="first")
    dropped = work.loc[~work.index.isin(keep.index)].copy()
    helper_cols = [c for c in work.columns if c.startswith("_")]
    return keep.drop(columns=helper_cols), dropped.drop(columns=helper_cols)


def add_normalized_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["peptide_sequence"] = out["Annotated Sequence"].map(extract_peptide)
    out["peptide_sequence_il_norm"] = out["peptide_sequence"].map(il_normalize)
    out["accession_list"] = out["Protein Accessions"].map(split_accessions)
    out["accession_list_str"] = out["accession_list"].map(lambda xs: "; ".join(xs))
    out["n_accessions"] = out["accession_list"].map(len)
    out["psm_abundance"] = to_number(out["Precursor Abundance"]).fillna(0.0)
    out["scan_key"] = out["Spectrum File"].astype(str) + ":" + out["First Scan"].astype(str)
    out["is_contaminant"] = out.apply(is_contaminant_row, axis=1)
    out["accession_class"] = out.apply(lambda r: accession_class(r["accession_list"], r["is_contaminant"]), axis=1)

    def clusters(accessions: list[str]) -> list[str]:
        return sorted({c for c in (parse_heavy_cluster_id(a) for a in accessions) if c is not None}, key=int)

    out["cluster_id_list"] = out["accession_list"].map(clusters)
    out["cluster_id_list_str"] = out["cluster_id_list"].map(lambda xs: "; ".join(xs))
    out["n_cluster_ids"] = out["cluster_id_list"].map(len)
    out["cluster_id"] = out["cluster_id_list"].map(lambda xs: xs[0] if len(xs) == 1 else "")
    return out


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    out = df.copy()
    for col in ["accession_list", "cluster_id_list"]:
        if col in out.columns:
            out = out.drop(columns=[col])
    out.to_csv(path, sep="\t", index=False)


def summarize(df: pd.DataFrame) -> dict[str, float | int]:
    return {
        "rows": int(len(df)),
        "unique_scan_keys": int(df["scan_key"].nunique()) if "scan_key" in df.columns else 0,
        "psm_abundance_sum": float(df["psm_abundance"].sum()) if "psm_abundance" in df.columns else 0.0,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Filter Proteome Discoverer PSMs to strict heavy-chain, single-lineage rows."
    )
    parser.add_argument("psm_file", help="Proteome Discoverer PSM .txt/.tsv export.")
    parser.add_argument("--out-dir", help="Output directory. Default: timestamped directory next to input.")
    args = parser.parse_args()

    input_path = Path(args.psm_file).expanduser()
    out_dir = make_out_dir(input_path, args.out_dir)
    stem = input_path.stem

    df = read_table(input_path)
    require_columns(df, REQUIRED_COLUMNS)
    df = add_normalized_columns(df)

    deduped, scan_duplicates = deduplicate_scans(df)
    deduped = add_normalized_columns(deduped)
    scan_duplicates = add_normalized_columns(scan_duplicates) if not scan_duplicates.empty else scan_duplicates

    contaminants = deduped[deduped["accession_class"] == "contaminant"].copy()
    light_chain = deduped[deduped["accession_class"] == "light_chain"].copy()
    mixed = deduped[deduped["accession_class"] == "mixed_heavy_nonheavy"].copy()
    non_bcrseq = deduped[deduped["accession_class"] == "non_bcrseq"].copy()
    heavy = deduped[deduped["accession_class"] == "heavy"].copy()
    multi_lineage = heavy[heavy["n_cluster_ids"] != 1].copy()
    filtered = heavy[heavy["n_cluster_ids"] == 1].copy()

    outputs = {
        "filtered_heavy_single_lineage": filtered,
        "excluded_scan_duplicates": scan_duplicates,
        "excluded_contaminants": contaminants,
        "excluded_light_chain": light_chain,
        "excluded_mixed_heavy_nonheavy": mixed,
        "excluded_non_bcrseq": non_bcrseq,
        "excluded_multi_lineage_heavy": multi_lineage,
    }

    output_paths: dict[str, str] = {}
    for name, frame in outputs.items():
        path = out_dir / f"{stem}_{name}.tsv"
        write_tsv(frame, path)
        output_paths[name] = str(path)

    log = {
        "command": " ".join(sys.argv),
        "input_file": str(input_path),
        "output_dir": str(out_dir),
        "initial": summarize(df),
        "after_scan_deduplication": summarize(deduped),
        "outputs": {name: summarize(frame) for name, frame in outputs.items()},
        "output_paths": output_paths,
        "notes": [
            "Scan deduplication key is Spectrum File + First Scan.",
            "Strict heavy-only rows require every accession to parse as a heavy BCRseq accession.",
            "Rows with more than one heavy ClusterID are excluded from quantification.",
        ],
    }
    log_path = out_dir / f"{stem}_filter_psms_log_{timestamp()}.json"
    log_path.write_text(json.dumps(log, indent=2))

    print(f"Wrote filtered PSMs: {output_paths['filtered_heavy_single_lineage']}")
    print(f"Wrote log: {log_path}")


if __name__ == "__main__":
    main()
