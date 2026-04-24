#!/usr/bin/env python3
"""
Quantify filtered heavy-chain PSMs by peptide and map them to exact BCRseq accessions.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd


REGIONS = ["fwr1_aa", "cdr1_aa", "fwr2_aa", "cdr2_aa", "fwr3_aa", "cdr3_aa", "fwr4_aa"]
REQUIRED_PSM_COLUMNS = [
    "peptide_sequence",
    "peptide_sequence_il_norm",
    "cluster_id",
    "accession_list_str",
    "psm_abundance",
    "Spectrum File",
    "First Scan",
]


def timestamp() -> str:
    return datetime.now().strftime("%Y-%m-%d_%H-%M-%S")


def make_out_dir(input_path: Path, out_dir: str | None) -> Path:
    if out_dir:
        path = Path(out_dir).expanduser()
    else:
        path = input_path.parent / f"quantify_map_peptides_{timestamp()}"
    path.mkdir(parents=True, exist_ok=True)
    return path


def read_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    df.columns = df.columns.str.strip('"')
    return df


def to_number(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def split_accessions(accessions: object) -> list[str]:
    if pd.isna(accessions):
        return []
    return [a.strip().strip('"') for a in str(accessions).split(";") if a.strip().strip('"')]


def il_normalize(seq: object) -> str:
    return re.sub(r"[IL]", "X", str(seq).upper())


def read_suffixes(path: Path) -> list[str]:
    seqs = [line.strip().upper() for line in path.read_text().splitlines() if line.strip()]
    if not seqs:
        raise ValueError(f"No suffix sequences found in {path}")
    return seqs


def normalize_bcr_id(raw: object, known_suffixes: set[str]) -> str:
    s = str(raw).strip().strip('"')
    if s.startswith("BCRseq_"):
        s = s[len("BCRseq_"):]
    parts = s.split("_")
    if parts and parts[-1].upper() in known_suffixes:
        parts = parts[:-1]
    return "_".join(parts)


def require_columns(df: pd.DataFrame, columns: list[str], label: str) -> None:
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required {label} column(s): " + ", ".join(missing))


def unique_join(values: pd.Series) -> str:
    out: list[str] = []
    seen: set[str] = set()
    for value in values.dropna().astype(str):
        for item in split_accessions(value):
            if item not in seen:
                seen.add(item)
                out.append(item)
    return "; ".join(out)


def summarize_filtered_psms(psm: pd.DataFrame) -> pd.DataFrame:
    psm = psm.copy()
    psm["psm_abundance_num"] = to_number(psm["psm_abundance"]).fillna(0.0)
    psm["pep_num"] = to_number(psm["Percolator PEP"]) if "Percolator PEP" in psm.columns else pd.NA
    psm["dm_num"] = to_number(psm["DeltaM [ppm]"]) if "DeltaM [ppm]" in psm.columns else pd.NA
    psm["scan_key"] = psm["Spectrum File"].astype(str) + ":" + psm["First Scan"].astype(str)

    abundance_by_injection = (
        psm.groupby(["peptide_sequence", "cluster_id", "Spectrum File"], dropna=False)["psm_abundance_num"]
        .sum()
        .reset_index()
    )
    abundance_by_injection["injection_pair"] = abundance_by_injection.apply(
        lambda r: f"{r['Spectrum File']}={r['psm_abundance_num']:.6g}", axis=1
    )
    inj_abundance_summary = (
        abundance_by_injection.groupby(["peptide_sequence", "cluster_id"], dropna=False)["injection_pair"]
        .agg(lambda xs: "; ".join(xs))
        .rename("peptide_abundance_by_injection")
        .reset_index()
    )
    avg_abundance_summary = (
        abundance_by_injection.groupby(["peptide_sequence", "cluster_id"], dropna=False)
        .agg(avg_abundance_per_injection=("psm_abundance_num", "mean"))
        .reset_index()
    )

    psm_count_by_injection = (
        psm.groupby(["peptide_sequence", "cluster_id", "Spectrum File"], dropna=False)
        .size()
        .reset_index(name="psm_count_injection")
    )
    psm_count_by_injection["injection_pair"] = psm_count_by_injection.apply(
        lambda r: f"{r['Spectrum File']}={int(r['psm_count_injection'])}", axis=1
    )
    inj_psm_summary = (
        psm_count_by_injection.groupby(["peptide_sequence", "cluster_id"], dropna=False)
        .agg(
            psm_count_by_injection=("injection_pair", lambda xs: "; ".join(xs)),
            n_injections_observed=("Spectrum File", "nunique"),
            avg_psm_per_injection=("psm_count_injection", "mean"),
        )
        .reset_index()
    )

    grouped = (
        psm.groupby(["peptide_sequence", "peptide_sequence_il_norm", "cluster_id"], dropna=False)
        .agg(
            peptide_abundance=("psm_abundance_num", "sum"),
            peptide_abundance_total=("psm_abundance_num", "sum"),
            psm_count=("psm_abundance_num", "size"),
            psm_count_total=("psm_abundance_num", "size"),
            scan_count=("scan_key", "nunique"),
            accession_list_str=("accession_list_str", unique_join),
            spectrum_files=("Spectrum File", lambda s: "; ".join(sorted(set(s.dropna().astype(str))))),
            mean_percolator_pep=("pep_num", "mean"),
            median_percolator_pep=("pep_num", "median"),
            mean_delta_m_ppm=("dm_num", "mean"),
        )
        .reset_index()
    )
    grouped = grouped.merge(inj_abundance_summary, on=["peptide_sequence", "cluster_id"], how="left")
    grouped["per_injection_abundance"] = grouped["peptide_abundance_by_injection"]
    grouped = grouped.merge(avg_abundance_summary, on=["peptide_sequence", "cluster_id"], how="left")
    grouped = grouped.merge(inj_psm_summary, on=["peptide_sequence", "cluster_id"], how="left")
    return grouped


def parse_bcrseq_data(bcrseq: pd.DataFrame, suffix: str, known_suffixes: set[str]) -> pd.DataFrame:
    df = bcrseq.copy()
    for name in REGIONS:
        if name not in df.columns:
            df[name] = ""
        df[name] = df[name].fillna("")

    def compute(row: pd.Series) -> pd.Series:
        positions: list[dict[str, Any]] = []
        start = 0
        seqs = []
        for region in REGIONS:
            seq = str(row[region]).upper()
            seqs.append(seq)
            end = start + len(seq)
            positions.append({"region": region, "start": start, "end": end})
            start = end
        positions.append({"region": "suffix", "start": start, "end": start + len(suffix)})
        full = "".join(seqs) + suffix
        return pd.Series({"full_bcr_seq": full, "region_positions": positions})

    df[["full_bcr_seq", "region_positions"]] = df.apply(compute, axis=1)
    df["normalized_full_bcr_seq"] = df["full_bcr_seq"].map(il_normalize)
    df["bcrseq_id_processed"] = df["bcrseq_id"].map(lambda s: normalize_bcr_id(s, known_suffixes))
    return df


def mapping_rows(peptides: pd.DataFrame, bcrseq: pd.DataFrame, known_suffixes: set[str]) -> pd.DataFrame:
    bcr_by_id = bcrseq.set_index("bcrseq_id_processed", drop=False)
    rows: list[dict[str, Any]] = []
    for _, pep in peptides.iterrows():
        accessions = split_accessions(pep["accession_list_str"])
        normalized_accessions = sorted({normalize_bcr_id(a, known_suffixes) for a in accessions})
        for norm_acc in normalized_accessions:
            if norm_acc not in bcr_by_id.index:
                rows.append({**pep.to_dict(), "bcrseq_id_processed": norm_acc, "match_found": False})
                continue
            matches = bcr_by_id.loc[[norm_acc]] if isinstance(bcr_by_id.loc[norm_acc], pd.DataFrame) else bcr_by_id.loc[[norm_acc]]
            for _, bcr in matches.iterrows():
                pep_norm = str(pep["peptide_sequence_il_norm"])
                text = str(bcr["normalized_full_bcr_seq"])
                starts = [m.start() for m in re.finditer(f"(?={re.escape(pep_norm)})", text)]
                if not starts:
                    rows.append({
                        **pep.to_dict(),
                        **bcr.to_dict(),
                        "match_found": False,
                        "matched_regions": "",
                        "region_coverage_percentages": "",
                        "match_start": pd.NA,
                        "match_end": pd.NA,
                        "cdr3_coverage_percentage": pd.NA,
                        "cdr3_cover_aa": pd.NA,
                    })
                    continue
                for start in starts:
                    end = start + len(pep_norm)
                    matched: list[str] = []
                    coverage: dict[str, float] = {}
                    cdr3_overlap = 0
                    cdr3_len = 0
                    for comp in bcr["region_positions"]:
                        rname, rs, re_ = comp["region"], comp["start"], comp["end"]
                        rlen = max(0, re_ - rs)
                        ov_s, ov_e = max(start, rs), min(end, re_)
                        if ov_s < ov_e:
                            matched.append(rname)
                            if rlen:
                                coverage[rname] = (ov_e - ov_s) / rlen * 100.0
                            if rname == "cdr3_aa":
                                cdr3_overlap = ov_e - ov_s
                                cdr3_len = rlen
                    rows.append({
                        **pep.to_dict(),
                        **bcr.to_dict(),
                        "match_found": True,
                        "matched_regions": ",".join(matched),
                        "region_coverage_percentages": ",".join(f"{k}:{v:.1f}%" for k, v in coverage.items()),
                        "match_start": start,
                        "match_end": end,
                        "cdr3_coverage_percentage": (cdr3_overlap / cdr3_len * 100.0) if cdr3_len else 0.0,
                        "cdr3_cover_aa": int(cdr3_overlap),
                    })
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Quantify filtered PSMs by peptide and map to exact BCRseq accessions."
    )
    parser.add_argument("filtered_psms", help="Filtered heavy single-lineage PSM TSV from filter_psms.py.")
    parser.add_argument("bcrseq_tsv", help="BCRseq annotation TSV with bcrseq_id.")
    parser.add_argument("suffix_txt", help="Suffix/extension text file, one AA sequence per line.")
    parser.add_argument("--out-dir", help="Output directory. Default: timestamped directory next to filtered PSMs.")
    args = parser.parse_args()

    psm_path = Path(args.filtered_psms).expanduser()
    bcr_path = Path(args.bcrseq_tsv).expanduser()
    suffix_path = Path(args.suffix_txt).expanduser()
    out_dir = make_out_dir(psm_path, args.out_dir)
    stem = psm_path.stem

    psm = read_table(psm_path)
    bcr = read_table(bcr_path)
    require_columns(psm, REQUIRED_PSM_COLUMNS, "PSM")
    require_columns(bcr, ["bcrseq_id", "ClusterID"] + REGIONS, "BCRseq")

    suffixes = read_suffixes(suffix_path)
    longest_suffix = max(suffixes, key=len)
    known_suffixes = set(suffixes)

    peptide_summary = summarize_filtered_psms(psm)
    bcr_prepped = parse_bcrseq_data(bcr, longest_suffix, known_suffixes)
    mapped = mapping_rows(peptide_summary, bcr_prepped, known_suffixes)

    peptide_path = out_dir / f"{stem}_quantified_peptides.tsv"
    mapped_path = out_dir / f"{stem}_mapped_peptides.tsv"
    peptide_summary.to_csv(peptide_path, sep="\t", index=False)
    mapped_to_write = mapped.drop(columns=["region_positions"], errors="ignore")
    mapped_to_write.to_csv(mapped_path, sep="\t", index=False)

    log = {
        "command": " ".join(sys.argv),
        "filtered_psms": str(psm_path),
        "bcrseq_tsv": str(bcr_path),
        "suffix_txt": str(suffix_path),
        "longest_suffix": longest_suffix,
        "input_psm_rows": int(len(psm)),
        "unique_spectrum_files": int(psm["Spectrum File"].dropna().astype(str).nunique()),
        "unique_peptides": int(len(peptide_summary)),
        "mapped_rows": int(len(mapped)),
        "matched_rows": int(mapped["match_found"].sum()) if "match_found" in mapped.columns else 0,
        "peptide_abundance_sum": float(peptide_summary["peptide_abundance"].sum()) if "peptide_abundance" in peptide_summary.columns else 0.0,
        "avg_abundance_per_injection_summary": {
            "min": float(peptide_summary["avg_abundance_per_injection"].min()) if "avg_abundance_per_injection" in peptide_summary.columns else 0.0,
            "median": float(peptide_summary["avg_abundance_per_injection"].median()) if "avg_abundance_per_injection" in peptide_summary.columns else 0.0,
            "max": float(peptide_summary["avg_abundance_per_injection"].max()) if "avg_abundance_per_injection" in peptide_summary.columns else 0.0,
        },
        "psm_count_total_sum": int(peptide_summary["psm_count_total"].sum()) if "psm_count_total" in peptide_summary.columns else 0,
        "avg_psm_per_injection_summary": {
            "min": float(peptide_summary["avg_psm_per_injection"].min()) if "avg_psm_per_injection" in peptide_summary.columns else 0.0,
            "median": float(peptide_summary["avg_psm_per_injection"].median()) if "avg_psm_per_injection" in peptide_summary.columns else 0.0,
            "max": float(peptide_summary["avg_psm_per_injection"].max()) if "avg_psm_per_injection" in peptide_summary.columns else 0.0,
        },
        "outputs": {"quantified_peptides": str(peptide_path), "mapped_peptides": str(mapped_path)},
    }
    log_path = out_dir / f"{stem}_quantify_map_log_{timestamp()}.json"
    log_path.write_text(json.dumps(log, indent=2))

    print(f"Wrote quantified peptides: {peptide_path}")
    print(f"Wrote mapped peptides: {mapped_path}")
    print(f"Wrote log: {log_path}")


if __name__ == "__main__":
    main()
