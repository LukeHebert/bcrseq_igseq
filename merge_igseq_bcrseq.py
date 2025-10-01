#!/usr/bin/env python3
"""
merge_igseq_bcrseq.py

Merge Proteome Discoverer peptide-spectrum matches (PSMs) with BCR cDNA/IgBLAST
(BCRseq) data, then locate where each peptide maps within the antibody variable
region (FWR/CDR), report regional coverage (including CDR3), and summarize
confident matches.

Highlights
----------
- Multiple PSM files:            -psm_files file1 file2 ...
- Multiple BCR chain types:      --bcrseq CHAIN:file [--bcrseq CHAIN:file ...]
- Per-chain extension/C1 files:  --suffix CHAIN:file [--suffix CHAIN:file ...]
- I/L equivalence during matching (I or L → X)
- Keeps peptides that map to exactly one BCR sequence
- Outputs a merged TSV and a brief text summary
"""

from __future__ import annotations

import argparse
import csv
import os
import re
from collections import defaultdict
from typing import Dict, Iterable, List, Sequence, Tuple

import pandas as pd


# -------------------------------
# Parsing and validation
# -------------------------------

def parse_chain_specs(specs: Sequence[str] | None, flag: str) -> Dict[str, List[str]]:
    """
    Parse values like ["IGH:fileA.tsv", "IGK:fileB.tsv", "IGK:fileC.tsv"] into
    a mapping {chain: [files...]}.

    Parameters
    ----------
    specs : sequence of CHAIN:file strings
    flag  : flag name for error messages (e.g., "--bcrseq", "--suffix")

    Returns
    -------
    dict mapping chain -> list of file paths
    """
    mapping: Dict[str, List[str]] = defaultdict(list)
    if not specs:
        return dict(mapping)
    for spec in specs:
        if ":" not in spec:
            raise ValueError(f"{flag} values must be CHAIN:FILE; got '{spec}'")
        chain, file = spec.split(":", 1)
        chain, file = chain.strip(), file.strip()
        if not chain or not file:
            raise ValueError(f"{flag} values must be CHAIN:FILE; got '{spec}'")
        mapping[chain].append(file)
    return dict(mapping)


def read_suffixes_per_chain(suffix_specs: Dict[str, List[str]]) -> Dict[str, List[str]]:
    """
    Read per-chain extension/C1 files (one AA sequence per line).

    Parameters
    ----------
    suffix_specs : dict chain -> list of files

    Returns
    -------
    dict chain -> list of suffix sequences (uppercased, stripped)
    """
    per_chain: Dict[str, List[str]] = {}
    for chain, files in suffix_specs.items():
        seqs: List[str] = []
        for path in files:
            with open(path) as fh:
                for line in fh:
                    s = line.strip()
                    if s:
                        seqs.append(s.upper())
        per_chain[chain] = seqs
    return per_chain


def require_suffixes_for_all_chains(
    bcrseq_chain_files: Dict[str, List[str]],
    chain_suffix_lists: Dict[str, List[str]],
) -> None:
    """
    Ensure that every chain provided via --bcrseq has at least one suffix entry
    provided via --suffix.

    Raises
    ------
    ValueError if any chain is missing or has zero suffix sequences.
    """
    missing = [ch for ch in bcrseq_chain_files if ch not in chain_suffix_lists]
    if missing:
        raise ValueError(
            "Missing --suffix for chain(s): " + ", ".join(sorted(missing))
        )
    empty = [ch for ch, seqs in chain_suffix_lists.items() if not seqs]
    if empty:
        raise ValueError(
            "Empty suffix list for chain(s) (each must contain ≥1 sequence): "
            + ", ".join(sorted(empty))
        )


def build_known_suffixes(chain_suffix_lists: Dict[str, List[str]]) -> set[str]:
    """
    Build a set of all suffix sequences observed across chains.

    Returns
    -------
    set of unique uppercase suffix sequences.
    """
    known: set[str] = set()
    for seqs in chain_suffix_lists.values():
        known.update(s.upper() for s in seqs if s)
    return known


# -------------------------------
# Normalization helpers
# -------------------------------

def normalize_bcr_id(raw: str, known_suffixes: set[str]) -> str:
    """
    Normalize a BCR identifier so PSM accessions and BCRseq IDs align.

    Rules
    -----
    - Strip leading 'BCRseq_' if present.
    - If the final underscore-delimited token equals a known suffix, drop it.

    Examples
    --------
    "BCRseq_XYZ_..._IGHJ4*01_STTA" -> "XYZ_..._IGHJ4*01"
    "BCRseq_XYZ_..._IGHJ4*01"      -> "XYZ_..._IGHJ4*01"
    """
    s = str(raw).strip().strip('"')
    if s.startswith("BCRseq_"):
        s = s[len("BCRseq_"):]
    parts = s.split("_")
    if parts and parts[-1].upper() in known_suffixes:
        parts = parts[:-1]
    return "_".join(parts)


# -------------------------------
# I/O
# -------------------------------

def load_data(
    psm_files: Sequence[str],
    bcrseq_chain_files: Dict[str, List[str]],
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load and concatenate PSM files and chain-annotated BCRseq files.

    Returns
    -------
    (pep_data, bcrseq_data) where
      pep_data   : concatenated PSMs
      bcrseq_data: concatenated BCR sequences with a 'chain_type' column
    """
    print("Loading data...")

    psm_frames = [
        pd.read_csv(path, sep="\t", quoting=csv.QUOTE_ALL, dtype=str)
        for path in psm_files
    ]
    for df in psm_frames:
        df.columns = df.columns.str.strip('"')
    pep_data = pd.concat(psm_frames, ignore_index=True)
    print(f"Combined PSM file row count: {pep_data.shape[0]}")

    bcr_frames: List[pd.DataFrame] = []
    for chain, files in bcrseq_chain_files.items():
        for path in files:
            df = pd.read_csv(path, sep="\t", dtype=str)
            df.columns = df.columns.str.strip('"')
            df["chain_type"] = chain
            bcr_frames.append(df)
    if not bcr_frames:
        raise ValueError("No BCRseq files provided.")
    bcrseq_data = pd.concat(bcr_frames, ignore_index=True, sort=False)
    print(f"Combined BCRseq file row count: {bcrseq_data.shape[0]}")

    return pep_data, bcrseq_data


# -------------------------------
# Core transforms
# -------------------------------

def filter_igseq_data(pep_data: pd.DataFrame, known_suffixes: set[str]) -> pd.DataFrame:
    """
    Keep PSM rows whose 'Protein Accessions' all begin with 'BCRseq_'.
    Explode per-accession and build the normalized merge key 'accession'.
    """
    pep = pep_data.copy()
    pep["Protein Accessions"] = pep["Protein Accessions"].astype(str)

    pep["accession_list"] = (
        pep["Protein Accessions"]
        .str.split(";")
        .apply(lambda xs: [acc.strip().strip('"') for acc in xs if acc.strip().strip('"')])
    )

    pep["valid_accessions"] = pep["accession_list"].apply(
        lambda accs: all(a.startswith("BCRseq_") for a in accs)
    )
    pep = pep[pep["valid_accessions"]].drop(columns=["valid_accessions"])
    print(
        "IgSeq data count after filtering out accessions not starting with 'BCRseq_':",
        pep.shape[0],
    )

    pep = pep.explode("accession_list", ignore_index=True)
    pep["accession"] = pep["accession_list"].apply(
        lambda s: normalize_bcr_id(s, known_suffixes)
    )
    return pep


def choose_longest_suffix_per_chain(
    chain_suffix_lists: Dict[str, List[str]]
) -> Dict[str, str]:
    """
    For each chain, choose the longest suffix string to append to the BCR full sequence.
    """
    return {chain: max(seqs, key=len) for chain, seqs in chain_suffix_lists.items()}


def parse_bcrseq_data(
    bcrseq_data: pd.DataFrame,
    chain_suffix_lists: Dict[str, List[str]],
    known_suffixes: set[str],
) -> pd.DataFrame:
    """
    Build a concatenated 'full_bcr_seq' per row from region AA fields and a
    chain-specific downstream extension. Record region boundaries and a version
    normalized for I/L equivalence. Produce 'bcrseq_id_processed' for merging.
    """
    df = bcrseq_data.copy()
    region_names = [
        "fwr1_aa", "cdr1_aa", "fwr2_aa", "cdr2_aa",
        "fwr3_aa", "cdr3_aa", "fwr4_aa",
    ]
    for name in region_names:
        if name not in df.columns:
            df[name] = ""
        df[name] = df[name].fillna("")

    longest_by_chain = choose_longest_suffix_per_chain(chain_suffix_lists)

    def compute_row(row: pd.Series) -> pd.Series:
        seqs = [str(row[name]).upper() for name in region_names]
        lengths = [len(s) for s in seqs]
        positions = []
        start = 0
        for name, ln in zip(region_names, lengths):
            end = start + ln
            positions.append({"region": name, "start": start, "end": end})
            start = end
        suffix = longest_by_chain[row["chain_type"]]
        positions.append({"region": "suffix", "start": start, "end": start + len(suffix)})
        full_seq = "".join(seqs) + suffix
        return pd.Series({"full_bcr_seq": full_seq, "region_positions": positions})

    df[["full_bcr_seq", "region_positions"]] = df.apply(compute_row, axis=1)

    print("Adding downstream peptide extensions (per chain):")
    for ch, suf in longest_by_chain.items():
        print(f"  chain={ch}: longest_suffix='{suf}'")
    print("BCRseq data count after adding chain-specific extension:", df.shape[0])

    df["normalized_full_bcr_seq"] = df["full_bcr_seq"].str.replace(r"[IL]", "X", regex=True)

    if "bcrseq_id" not in df.columns:
        raise ValueError("BCRseq input must include a 'bcrseq_id' column.")

    df["bcrseq_id_processed"] = df["bcrseq_id"].apply(
        lambda s: normalize_bcr_id(s, known_suffixes)
    )
    return df


def merge_datasets(pep_data: pd.DataFrame, bcrseq_data: pd.DataFrame) -> pd.DataFrame:
    """
    Inner-join peptides to BCR sequences on their normalized IDs.
    """
    print("Merging data on processed 'accession' and 'bcrseq_id'...")
    merged = pep_data.merge(
        bcrseq_data, left_on="accession", right_on="bcrseq_id_processed"
    )
    merged = merged.loc[:, ~merged.columns.duplicated()]  # ensure unique column names
    print(f"Merged data count: {merged.shape[0]}")
    return merged


def find_matched_regions(merged_data: pd.DataFrame) -> pd.DataFrame:
    """
    Extract peptide strings from 'Annotated Sequence', perform normalized
    substring search against normalized BCR sequences, and compute regional
    coverage (% of each BCR region overlapped by the peptide).
    """
    if merged_data.empty:
        print("No merged rows; skipping match search.")
        cols = [
            "matched_regions", "region_coverage_percentages", "match_start",
            "match_end", "cdr3_coverage_percentage", "match_found",
        ]
        return merged_data.assign(
            **{c: pd.Series(dtype=object) for c in cols}
        )[:0]

    df = merged_data.copy()

    if "Annotated Sequence" not in df.columns:
        raise ValueError("PSM input must include 'Annotated Sequence' column.")

    df["peptide_sequence"] = df["Annotated Sequence"].apply(
        lambda x: str(x).split(".")[1].upper() if "." in str(x) else str(x).upper()
    )

    df["normalized_peptide_sequence"] = df["peptide_sequence"].str.replace(r"[IL]", "X", regex=True)
    df["normalized_full_bcr_seq"] = df["full_bcr_seq"].str.replace(r"[IL]", "X", regex=True)

    def regions_for_row(row: pd.Series) -> pd.Series:
        text = row["normalized_full_bcr_seq"]
        pep = row["normalized_peptide_sequence"]
        positions = row["region_positions"]
        matches = [m.start() for m in re.finditer(f"(?={re.escape(pep)})", text)]
        if not matches:
            return pd.Series({
                "match_info_list": [],
                "match_found": False,
            })
        info: List[dict] = []
        pep_len = len(pep)
        for start in matches:
            end = start + pep_len
            matched_regions: List[str] = []
            coverage: Dict[str, float] = {}
            cdr3_overlap, cdr3_len = 0, 0
            for comp in positions:
                rname, rs, re_ = comp["region"], comp["start"], comp["end"]
                rlen = max(0, re_ - rs)
                ov_s, ov_e = max(start, rs), min(end, re_)
                if ov_s < ov_e:
                    matched_regions.append(rname)
                    ov_len = ov_e - ov_s
                    if rlen > 0:
                        coverage[rname] = (ov_len / rlen) * 100.0
                    if rname == "cdr3_aa":
                        cdr3_overlap, cdr3_len = ov_len, rlen
            cov_str = ",".join(f"{k}:{v:.1f}%" for k, v in coverage.items())
            cdr3_cov = (cdr3_overlap / cdr3_len) * 100.0 if cdr3_len > 0 else 0.0
            info.append({
                "matched_regions": ",".join(matched_regions),
                "region_coverage_percentages": cov_str,
                "match_start": start,
                "match_end": end,
                "cdr3_coverage_percentage": cdr3_cov,
            })
        return pd.Series({"match_info_list": info, "match_found": True})

    out = df.join(df.apply(regions_for_row, axis=1))

    # Expand per-match rows and extract fields
    out = out[out["match_found"]].explode("match_info_list")
    if out.empty:
        return out.assign(
            matched_regions="", region_coverage_percentages="",
            match_start=pd.NA, match_end=pd.NA, cdr3_coverage_percentage=pd.NA
        )

    get = lambda key: out["match_info_list"].apply(lambda d: d[key])
    out = out.assign(
        matched_regions=get("matched_regions"),
        region_coverage_percentages=get("region_coverage_percentages"),
        match_start=get("match_start"),
        match_end=get("match_end"),
        cdr3_coverage_percentage=get("cdr3_coverage_percentage"),
    ).drop(columns=["match_info_list"])

    print("Number of rows where peptide sequence was found in full_bcr_seq:", out.shape[0])
    return out


def filter_matches(matched_data: pd.DataFrame) -> pd.DataFrame:
    """
    Remove matches occurring only in framework/suffix regions; keep those
    overlapping any CDR region.
    """
    if matched_data.empty:
        return matched_data

    def fwr_or_suffix_only(s: str) -> bool:
        regions = [r for r in str(s).split(",") if r]
        return all(not r.startswith("cdr") for r in regions) if regions else True

    df = matched_data.copy()
    initial = df.shape[0]
    df["is_only_fwr_or_suffix"] = df["matched_regions"].map(fwr_or_suffix_only)
    df = df[~df["is_only_fwr_or_suffix"]].drop(columns=["is_only_fwr_or_suffix"])
    print(
        "Number of matches removed that occur only in FWR or suffix regions:",
        initial - df.shape[0],
    )
    print(
        "Number of matches remaining after removing FWR/suffix-only matches:",
        df.shape[0],
    )
    return df


def finalize_data(filtered_data: pd.DataFrame) -> pd.DataFrame:
    """
    Count distinct BCR sequences per peptide and keep those with exactly one.
    """
    if filtered_data.empty:
        print("No matches to finalize.")
        return filtered_data

    counts = (
        filtered_data.groupby("Annotated Sequence")["bcrseq_id"]
        .nunique()
        .rename("bcrseq_match_count")
        .reset_index()
    )
    print("Statistical distribution of peptide-BCRseq match counts before filtering:")
    print(counts["bcrseq_match_count"].describe())

    df = filtered_data.merge(counts, on="Annotated Sequence", how="left")
    final = df[df["bcrseq_match_count"] == 1]
    print(
        "Number of matches remaining after filtering peptides matching to only one BCRseq sequence:",
        final.shape[0],
    )
    print("Statistical distribution of peptide-BCRseq match counts after filtering:")
    if not final.empty:
        print(final.groupby("Annotated Sequence")["bcrseq_id"].nunique().describe())
    else:
        print("empty")
    return final


# -------------------------------
# Orchestration
# -------------------------------

def main(psm_files: Sequence[str], bcrseq_specs: Sequence[str], suffix_specs: Sequence[str]) -> None:
    """
    Orchestrate loading, normalization, merging, matching, filtering,
    and writing outputs.
    """
    bcrseq_chain_files = parse_chain_specs(bcrseq_specs, "--bcrseq")
    suffix_chain_files = parse_chain_specs(suffix_specs, "--suffix")
    chain_suffix_lists = read_suffixes_per_chain(suffix_chain_files)
    require_suffixes_for_all_chains(bcrseq_chain_files, chain_suffix_lists)
    known_suffixes = build_known_suffixes(chain_suffix_lists)

    pep_data, bcrseq_data = load_data(psm_files, bcrseq_chain_files)
    print(f"Initial IgSeq data count: {pep_data.shape[0]}")

    pep_data = filter_igseq_data(pep_data, known_suffixes)
    bcrseq_data = parse_bcrseq_data(bcrseq_data, chain_suffix_lists, known_suffixes)

    # Debug sample to help verify keys if needed (comment out in routine runs)
    # print("PSM normalized accession sample:", pep_data["accession"].dropna().unique()[:5])
    # print("BCR normalized id sample:", bcrseq_data["bcrseq_id_processed"].dropna().unique()[:5])

    merged_data = merge_datasets(pep_data, bcrseq_data)
    matched_data = find_matched_regions(merged_data)
    filtered_data = filter_matches(matched_data)
    final_data = finalize_data(filtered_data)

    base = os.path.splitext(psm_files[0])[0]
    out_tsv = base + "_merged.tsv"
    final_data.to_csv(out_tsv, sep="\t", index=False)
    print(f"Data saved to {out_tsv}")

    total_full = final_data["full_bcr_seq"].shape[0]
    uniq_full = final_data["full_bcr_seq"].nunique()
    uniq_bcr = final_data["bcrseq_id"].nunique()
    cdr3_only = final_data[final_data["matched_regions"] == "cdr3_aa"].shape[0]
    cdr3_plus = final_data[
        final_data["matched_regions"].str.contains("cdr3_aa", na=False)
        & (final_data["matched_regions"] != "cdr3_aa")
    ].shape[0]

    summary = base + "_match_summary.txt"
    with open(summary, "w") as f:
        f.write(f"Initial IgSeq data count: {pep_data.shape[0]}\n")
        f.write(
            "IgSeq data count after filtering out accessions not starting with 'BCRseq_': "
            f"{pep_data.shape[0]}\n"
        )
        f.write(f"Merged data count: {merged_data.shape[0]}\n")
        f.write(
            "Number of matches before removing FWR/suffix-only matches: "
            f"{matched_data.shape[0]}\n"
        )
        f.write(
            "Number of matches after removing FWR/suffix-only matches: "
            f"{filtered_data.shape[0]}\n"
        )
        f.write(
            "Number of matches remaining after filtering peptides matching to only one BCRseq sequence: "
            f"{final_data.shape[0]}\n"
        )
        f.write(f"Number of total full_bcr_seq values in final matched data: {total_full}\n")
        f.write(f"Number of unique full_bcr_seq values in final matched data: {uniq_full}\n")
        f.write(f"Number of unique bcrseq_id values in final matched data: {uniq_bcr}\n")
        f.write(f"Number of matches within cdr3_aa only: {cdr3_only}\n")
        f.write(f"Number of matches within cdr3_aa plus other regions: {cdr3_plus}\n")
    print(f"Summary saved to {summary}")


# -------------------------------
# CLI
# -------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Merge PSMs with BCRseq and compute peptide regional coverage."
    )
    parser.add_argument(
        "-psm_files",
        nargs="+",
        required=True,
        help="Proteome Discoverer PSM TXT/TSV file(s).",
    )
    parser.add_argument(
        "--bcrseq",
        dest="bcrseq_specs",
        nargs="+",
        metavar="CHAIN:FILE",
        required=True,
        help="Repeatable. One or more BCRseq TSV files, annotated with chain type, "
             "e.g., --bcrseq IGH:heavy.tsv IGK:kappa.tsv",
    )
    parser.add_argument(
        "--suffix",
        dest="suffix_specs",
        nargs="+",
        metavar="CHAIN:FILE",
        required=True,
        help="Repeatable. Map each chain type to an extension/C1 file (one AA sequence per line). "
             "The longest sequence per chain is appended when forming full_bcr_seq. "
             "Example: --suffix IGH:heavy_suffixes.txt IGK:kappa_suffixes.txt",
    )
    args = parser.parse_args()
    main(args.psm_files, args.bcrseq_specs, args.suffix_specs)
