#!/usr/bin/env python3
"""
merge_igseq_bcrseq.py

Merge Proteome Discoverer peptide-spectrum matches (PSMs) with BCR cDNA/IgBLAST
(BCRseq) data (if heavy chain, must already be clustered with ClusterID column).
Upon merging, the location where each peptide maps within BCRseq antibody regions
(FWR/CDR) is determined and region coverage (with CDR3 emphasis) is computed.
Chain-specific results are written to a new 'merge_igseq_bcrseq' directory next
to the input PSM file.

Features
--------
- Takes ≥1 PSM file(s):         -psm_files file1 file2 ...
- Takes ≥1 BCRseq files:        --bcrseq CHAIN:file [--bcrseq CHAIN:file ...]
- Takes extension/C seq file(s): --suffix CHAIN:file [--suffix CHAIN:file ...]
- I/L equivalence for peptide region mapping (I and L collapsed to 'X')
- Optional handling of mixed PSM accessions (default: allow mixed; use --no-mixed-accessions to require all-BCRseq)
- Chain-specific TSV outputs (+ chain-aware, timestamped log)
- Adds cdr3_cover_aa (AA count overlapped in CDR3) and psm_clusterid_list
- Writes all outputs to a subdirectory named 'merge_igseq_bcrseq' next to the first PSM
- All console prints are also captured into the log file
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from collections import defaultdict
from datetime import datetime
from typing import Dict, List, Sequence, Tuple

import pandas as pd


# -------------------------------
# Small print+log helper
# -------------------------------

_LOG_LINES: List[str] = []

def log_print(*args, sep=" ", end="\n"):
    """Print to stdout and also record the message for the run log."""
    msg = sep.join(str(a) for a in args)
    # Normalize newline handling: splitlines to capture multi-line prints faithfully
    lines = (msg + ("" if end == "" else end)).splitlines()
    for line in lines:
        if line == "" and end == "":
            # Exact print behavior: don't echo trailing blank when end=""
            continue
        _LOG_LINES.append(line)
    print(msg, end=end)


# -------------------------------
# Parsing and validation
# -------------------------------

def parse_chain_specs(specs: Sequence[str] | None, flag: str) -> Dict[str, List[str]]:
    """Parse values like ['IGH:fileA.tsv', 'IGK:fileB.tsv'] into {chain: [files...]}."""
    mapping: Dict[str, List[str]] = defaultdict(list)
    if not specs:
        return dict(mapping)
    for spec in specs:
        if ":" not in spec:
            raise ValueError(f"{flag} values must be CHAIN:FILE; got '{spec}'")
        chain, file = (s.strip() for s in spec.split(":", 1))
        if not chain or not file:
            raise ValueError(f"{flag} values must be CHAIN:FILE; got '{spec}'")
        mapping[chain].append(file)
    return dict(mapping)


def read_suffixes_per_chain(suffix_specs: Dict[str, List[str]]) -> Dict[str, List[str]]:
    """Read per-chain extension/C1 files (one AA sequence per line)."""
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
    """Ensure every chain in --bcrseq has at least one suffix sequence in --suffix."""
    missing = [ch for ch in bcrseq_chain_files if ch not in chain_suffix_lists]
    if missing:
        raise ValueError("Missing --suffix for chain(s): " + ", ".join(sorted(missing)))
    empty = [ch for ch, seqs in chain_suffix_lists.items() if not seqs]
    if empty:
        raise ValueError(
            "Empty suffix list for chain(s) (each must contain ≥1 sequence): "
            + ", ".join(sorted(empty))
        )


def build_known_suffixes(chain_suffix_lists: Dict[str, List[str]]) -> set[str]:
    """Build a set of all suffix sequences observed across chains (uppercased)."""
    known: set[str] = set()
    for seqs in chain_suffix_lists.values():
        known.update(s.upper() for s in seqs if s)
    return known


# -------------------------------
# Normalization helpers
# -------------------------------

def normalize_bcr_id(raw: str, known_suffixes: set[str]) -> str:
    """
    Normalize a BCR identifier so PSM accessions and BCRseq IDs align:
      - strip leading 'BCRseq_' if present
      - drop trailing '_<suffix>' if <suffix> is in known_suffixes
    """
    s = str(raw).strip().strip('"')
    if s.startswith("BCRseq_"):
        s = s[len("BCRseq_"):]
    parts = s.split("_")
    if parts and parts[-1].upper() in known_suffixes:
        parts = parts[:-1]
    return "_".join(parts)


def il_normalize(seq: str) -> str:
    """Replace I and L with X for equivalence matching."""
    return re.sub(r"[IL]", "X", str(seq))


def extract_peptide(annot: str) -> str:
    """Extract the peptide from PD’s 'Annotated Sequence' (prefix.PEPTIDE.suffix or raw)."""
    s = str(annot)
    parts = s.split(".")
    return (parts[1] if len(parts) >= 3 else s).upper()


# -------------------------------
# I/O
# -------------------------------

def load_data(
    psm_files: Sequence[str],
    bcrseq_chain_files: Dict[str, List[str]],
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load and concatenate PSM files and chain-annotated BCRseq files."""
    log_print("Loading data...")
    psm_frames = [
        pd.read_csv(path, sep="\t", quoting=csv.QUOTE_ALL, dtype=str) for path in psm_files
    ]
    for df in psm_frames:
        df.columns = df.columns.str.strip('"')
    pep_data = pd.concat(psm_frames, ignore_index=True)
    log_print(f"Combined PSM file row count: {pep_data.shape[0]}")

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
    log_print(f"Combined BCRseq file row count: {bcrseq_data.shape[0]}")
    return pep_data, bcrseq_data


# -------------------------------
# Core transforms
# -------------------------------

def prep_igseq_data(
    pep_data: pd.DataFrame,
    known_suffixes: set[str],
    allow_mixed_accessions: bool,
) -> pd.DataFrame:
    """
    Prepare PSM data:
      - split/trim Protein Accessions
      - optionally keep rows that mix BCRseq_ and non-BCRseq accessions
      - explode to one row per BCRseq accession
      - normalize 'accession' key for merging
      - de-dupe exploded rows by a stable PSM identity + normalized accession
    """
    pep = pep_data.copy()

    # Initial unique counts (pre-any processing)
    init_unique_psm = pep["PSMs Peptide ID"].nunique(dropna=True) if "PSMs Peptide ID" in pep.columns else 0
    init_unique_pep = pep["Annotated Sequence"].map(extract_peptide).nunique(dropna=True) if "Annotated Sequence" in pep.columns else 0
    log_print(f"Initial unique PSMs Peptide ID: {init_unique_psm}")
    log_print(f"Initial unique peptide_sequence: {init_unique_pep}")

    pep["Protein Accessions"] = pep["Protein Accessions"].astype(str)
    pep["accession_list"] = (
        pep["Protein Accessions"]
        .str.split(";")
        .apply(lambda xs: [acc.strip().strip('"') for acc in xs if acc.strip().strip('"')])
    )

    if allow_mixed_accessions:
        # Keep mixed rows, but drop non-BCRseq tokens for exploding/merging
        pep["accession_list"] = pep["accession_list"].apply(
            lambda accs: [a for a in accs if a.startswith("BCRseq_")]
        )
        pep = pep[pep["accession_list"].apply(len) > 0]
    else:
        # Require all tokens be BCRseq_
        pep = pep[pep["accession_list"].apply(lambda accs: all(a.startswith("BCRseq_") for a in accs))]

    log_print("IgSeq data count after PSM filtering:", pep.shape[0])

    # Explode to per-accession rows
    pep = pep.explode("accession_list", ignore_index=True)

    # Build normalized key used for merging
    pep["accession"] = pep["accession_list"].apply(lambda s: normalize_bcr_id(s, known_suffixes))

    # De-dupe exploded PSM rows that only differ by suffix
    psm_id = ["PSMs Peptide ID"] if "PSMs Peptide ID" in pep.columns else []
    ident_cols = list(psm_id)
    if "Annotated Sequence" in pep.columns and "Annotated Sequence" not in ident_cols:
        ident_cols.append("Annotated Sequence")
    subset_cols = ident_cols + ["accession"]
    pep = pep.drop_duplicates(subset=subset_cols, keep="first")

    return pep


def choose_longest_suffix_per_chain(
    chain_suffix_lists: Dict[str, List[str]]
) -> Dict[str, str]:
    """For each chain, choose the longest suffix string to append to the BCR full sequence."""
    return {chain: max(seqs, key=len) for chain, seqs in chain_suffix_lists.items()}


def parse_bcrseq_data(
    bcrseq_data: pd.DataFrame,
    chain_suffix_lists: Dict[str, List[str]],
    known_suffixes: set[str],
) -> pd.DataFrame:
    """
    Build 'full_bcr_seq' from region AA fields and chain-specific extension.
    Record region boundaries and normalized version (I/L collapsed).
    Produce 'bcrseq_id_processed' for merging.
    """
    df = bcrseq_data.copy()
    regions = ["fwr1_aa", "cdr1_aa", "fwr2_aa", "cdr2_aa", "fwr3_aa", "cdr3_aa", "fwr4_aa"]
    for name in regions:
        if name not in df.columns:
            df[name] = ""
        df[name] = df[name].fillna("")

    longest = choose_longest_suffix_per_chain(chain_suffix_lists)

    def compute_row(row: pd.Series) -> pd.Series:
        seqs = [str(row[name]).upper() for name in regions]
        lengths = [len(s) for s in seqs]
        positions, start = [], 0
        for name, ln in zip(regions, lengths):
            end = start + ln
            positions.append({"region": name, "start": start, "end": end})
            start = end
        suffix = longest[row["chain_type"]]
        positions.append({"region": "suffix", "start": start, "end": start + len(suffix)})
        full = "".join(seqs) + suffix
        return pd.Series({"full_bcr_seq": full, "region_positions": positions})

    df[["full_bcr_seq", "region_positions"]] = df.apply(compute_row, axis=1)
    log_print("Adding downstream peptide extensions (per chain):")
    for ch, suf in longest.items():
        log_print(f"  chain={ch}: longest_suffix='{suf}'")
    log_print("BCRseq data count after adding chain-specific extension:", df.shape[0])

    if "bcrseq_id" not in df.columns:
        raise ValueError("BCRseq input must include a 'bcrseq_id' column.")
    df["normalized_full_bcr_seq"] = df["full_bcr_seq"].map(il_normalize)
    df["bcrseq_id_processed"] = df["bcrseq_id"].apply(lambda s: normalize_bcr_id(s, known_suffixes))
    return df


def merge_datasets(pep_data: pd.DataFrame, bcrseq_data: pd.DataFrame) -> pd.DataFrame:
    """Inner-join peptides to BCR sequences on normalized IDs."""
    log_print("Merging data on processed 'accession' and 'bcrseq_id'...")
    merged = pep_data.merge(bcrseq_data, left_on="accession", right_on="bcrseq_id_processed")
    merged = merged.loc[:, ~merged.columns.duplicated()]
    log_print(f"Merged data count: {merged.shape[0]}")
    return merged


def find_matched_regions(merged_data: pd.DataFrame) -> pd.DataFrame:
    """
    Extract peptide string from 'Annotated Sequence', perform I/L-collapsed
    substring search vs 'normalized_full_bcr_seq', compute per-region coverage
    (percent of region overlapped by the peptide match). Also computes:
      - cdr3_coverage_percentage
      - cdr3_cover_aa (integer AA overlap length within CDR3)
    """
    if merged_data.empty:
        log_print("No merged rows; skipping match search.")
        cols = ["matched_regions", "region_coverage_percentages", "match_start",
                "match_end", "cdr3_coverage_percentage", "cdr3_cover_aa", "match_found"]
        return merged_data.assign(**{c: pd.Series(dtype=object) for c in cols})[:0]

    df = merged_data.copy()

    if "Annotated Sequence" not in df.columns:
        raise ValueError("PSM input must include 'Annotated Sequence' column.")

    df["peptide_sequence"] = df["Annotated Sequence"].apply(extract_peptide)
    df["normalized_peptide_sequence"] = df["peptide_sequence"].map(il_normalize)
    df["normalized_full_bcr_seq"] = df["full_bcr_seq"].map(il_normalize)

    def regions_for_row(row: pd.Series) -> pd.Series:
        text, pep = row["normalized_full_bcr_seq"], row["normalized_peptide_sequence"]
        positions = row["region_positions"]
        hits = [m.start() for m in re.finditer(f"(?={re.escape(pep)})", text)]
        if not hits:
            return pd.Series({"match_info_list": [], "match_found": False})
        info = []
        p_len = len(pep)
        for start in hits:
            end = start + p_len
            matched, coverage = [], {}
            cdr3_overlap = cdr3_len = 0
            for comp in positions:
                rname, rs, re_ = comp["region"], comp["start"], comp["end"]
                rlen = max(0, re_ - rs)
                ov_s, ov_e = max(start, rs), min(end, re_)
                if ov_s < ov_e:
                    matched.append(rname)
                    if rlen > 0:
                        coverage[rname] = (ov_e - ov_s) / rlen * 100.0
                    if rname == "cdr3_aa":
                        cdr3_overlap, cdr3_len = (ov_e - ov_s), rlen
            cov_str = ",".join(f"{k}:{v:.1f}%" for k, v in coverage.items())
            cdr3_cov = (cdr3_overlap / cdr3_len * 100.0) if cdr3_len > 0 else 0.0
            info.append({
                "matched_regions": ",".join(matched),
                "region_coverage_percentages": cov_str,
                "match_start": start,
                "match_end": end,
                "cdr3_coverage_percentage": cdr3_cov,
                "cdr3_cover_aa": int(cdr3_overlap),
            })
        return pd.Series({"match_info_list": info, "match_found": True})

    out = df.join(df.apply(regions_for_row, axis=1))
    out = out[out["match_found"]].explode("match_info_list")
    if out.empty:
        return out.assign(
            matched_regions="", region_coverage_percentages="",
            match_start=pd.NA, match_end=pd.NA,
            cdr3_coverage_percentage=pd.NA, cdr3_cover_aa=pd.NA
        )

    get = lambda key: out["match_info_list"].apply(lambda d: d[key])
    out = out.assign(
        matched_regions=get("matched_regions"),
        region_coverage_percentages=get("region_coverage_percentages"),
        match_start=get("match_start"),
        match_end=get("match_end"),
        cdr3_coverage_percentage=get("cdr3_coverage_percentage"),
        cdr3_cover_aa=get("cdr3_cover_aa"),
    ).drop(columns=["match_info_list"])

    log_print("Number of rows where peptide sequence was found in full_bcr_seq:", out.shape[0])
    return out


# -------------------------------
# Output helpers
# -------------------------------

NONVERBOSE_COLS = [
    # PSM
    "PSMs Peptide ID", "Annotated Sequence", "Protein Accessions",
    "Precursor Abundance", "DeltaCn", "DeltaM [ppm]", "Percolator PEP",
    # BCR (selection)
    "sequence_id", "sequence", "sequence_aa", "v_call", "d_call", "j_call", "c_call",
    "fwr1", "fwr1_aa", "cdr1", "cdr1_aa", "fwr2", "fwr2_aa", "cdr2", "cdr2_aa",
    "fwr3", "fwr3_aa", "fwr4", "fwr4_aa", "cdr3", "cdr3_aa",
    "c_sequence_alignment", "c_sequence_alignment_aa",
    "v_identity", "d_identity", "j_identity", "c_identity", "ClusterID", "Collapsed",
    # Derived
    "full_bcr_seq", "normalized_full_bcr_seq",
    "peptide_sequence", "normalized_peptide_sequence",
    "region_coverage_percentages", "cdr3_coverage_percentage", "cdr3_cover_aa",
    "psm_clusterid_list",
]

def ensure_columns(df: pd.DataFrame, cols: List[str]) -> pd.DataFrame:
    """Add any missing columns in `cols` as empty strings to allow column selection."""
    missing = [c for c in cols if c not in df.columns]
    if missing:
        df = df.assign(**{c: "" for c in missing})
    return df


def _parse_clusterids_from_accessions(acc_str: str) -> str:
    """
    From a semicolon-separated 'Protein Accessions' string, keep the 3rd '_' token
    (index 2) from each BCRseq_* entry, returning a semicolon-separated ClusterID list.
    """
    parts = [p.strip() for p in str(acc_str).split(";") if p.strip()]
    out: List[str] = []
    for p in parts:
        if p.startswith("BCRseq_"):
            toks = p.split("_")
            if len(toks) >= 3:
                out.append(toks[2])
    return "; ".join(out)


def write_chain_outputs(
    data_fullcols: pd.DataFrame,
    out_dir: str,
    base_stem: str,
    verbose: bool,
) -> Dict[str, str]:
    """Write one TSV per chain_type into `out_dir` and return {chain: out_path}."""
    if data_fullcols.empty:
        return {}
    paths: Dict[str, str] = {}
    for chain, sub in data_fullcols.groupby("chain_type", dropna=False):
        ch = str(chain) if pd.notna(chain) else "NA"
        out_path = os.path.join(out_dir, f"{base_stem}_merged_{ch}.tsv")
        to_write = sub if verbose else ensure_columns(sub, NONVERBOSE_COLS)[NONVERBOSE_COLS]
        to_write.to_csv(out_path, sep="\t", index=False)
        log_print(f"Wrote: {out_path} (rows={to_write.shape[0]})")
        paths[ch] = out_path
    return paths


def write_log(
    out_dir: str,
    base_stem: str,
    cmdline: str,
    pep_count_after_psm_filter: int,
    merged_data: pd.DataFrame,
    matched_data: pd.DataFrame,
    init_unique_psm: int,
    init_unique_pepseq: int,
    final_unique_psm: int,
    final_unique_pepseq: int,
) -> str:
    """
    Write a timestamped log file with counts. Includes:
      - exact command line
      - initial/final unique counts
      - merged/matched counts (per chain)
      - full console transcript captured via log_print
    """
    ts = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    log_path = os.path.join(out_dir, f"{base_stem}_log_{ts}.txt")

    def chain_counts(df: pd.DataFrame, label: str) -> List[str]:
        if df.empty:
            return [f"{label}: 0 total (no data)"]
        lines = [f"{label}: total={df.shape[0]}"]
        for ch, sub in df.groupby("chain_type", dropna=False):
            chs = str(ch) if pd.notna(ch) else "NA"
            lines.append(f"  {chs}: {sub.shape[0]}")
        return lines

    with open(log_path, "w") as f:
        f.write("Command line:\n")
        f.write(cmdline + "\n\n")
        f.write(f"Initial unique PSMs Peptide ID: {init_unique_psm}\n")
        f.write(f"Initial unique peptide_sequence: {init_unique_pepseq}\n\n")
        f.write(f"IgSeq data count after PSM filtering: {pep_count_after_psm_filter}\n")
        for line in chain_counts(merged_data, "Merged rows"):
            f.write(line + "\n")
        for line in chain_counts(matched_data, "Matched rows (all regions considered)"):
            f.write(line + "\n")
        f.write("\n")
        f.write(f"Final unique PSMs Peptide ID: {final_unique_psm}\n")
        f.write(f"Final unique peptide_sequence: {final_unique_pepseq}\n")
        f.write("\n--- Console output transcript ---\n")
        for line in _LOG_LINES:
            f.write(line + "\n")

    log_print(f"Log saved to {log_path}")
    return log_path


# -------------------------------
# Orchestration
# -------------------------------

def main(
    psm_files: Sequence[str],
    bcrseq_specs: Sequence[str],
    suffix_specs: Sequence[str],
    allow_mixed_accessions: bool,
    verbose: bool,
) -> None:
    """Load, normalize, merge, match, and write outputs."""
    # Output directory (next to first PSM)
    first_psm_dir = os.path.dirname(os.path.abspath(psm_files[0]))
    out_dir = os.path.join(first_psm_dir, "merge_igseq_bcrseq")
    os.makedirs(out_dir, exist_ok=True)

    bcrseq_chain_files = parse_chain_specs(bcrseq_specs, "--bcrseq")
    suffix_chain_files = parse_chain_specs(suffix_specs, "--suffix")
    chain_suffix_lists = read_suffixes_per_chain(suffix_chain_files)
    require_suffixes_for_all_chains(bcrseq_chain_files, chain_suffix_lists)
    known_suffixes = build_known_suffixes(chain_suffix_lists)

    pep_data, bcrseq_data = load_data(psm_files, bcrseq_chain_files)

    # Initial unique counts (before any prep)
    init_unique_psm = pep_data["PSMs Peptide ID"].nunique(dropna=True) if "PSMs Peptide ID" in pep_data.columns else 0
    init_unique_pepseq = pep_data["Annotated Sequence"].map(extract_peptide).nunique(dropna=True) if "Annotated Sequence" in pep_data.columns else 0
    log_print(f"Initial unique PSMs Peptide ID: {init_unique_psm}")
    log_print(f"Initial unique peptide_sequence: {init_unique_pepseq}")

    pep_data = prep_igseq_data(pep_data, known_suffixes, allow_mixed_accessions)
    pep_count_after_psm_filter = pep_data.shape[0]

    bcrseq_data = parse_bcrseq_data(bcrseq_data, chain_suffix_lists, known_suffixes)

    merged_data = merge_datasets(pep_data, bcrseq_data)

    # Parse cluster IDs from the original PSM accessions (per-row)
    merged_data["psm_clusterid_list"] = merged_data["Protein Accessions"].map(_parse_clusterids_from_accessions)

    matched_data = find_matched_regions(merged_data)

    # Final uniques (post-match)
    final_unique_psm = matched_data["PSMs Peptide ID"].nunique(dropna=True) if "PSMs Peptide ID" in matched_data.columns else 0
    final_unique_pepseq = matched_data["peptide_sequence"].nunique(dropna=True) if "peptide_sequence" in matched_data.columns else 0
    log_print(f"Final unique PSMs Peptide ID: {final_unique_psm}")
    log_print(f"Final unique peptide_sequence: {final_unique_pepseq}")

    # Write per-chain outputs
    base_stem = os.path.splitext(os.path.basename(psm_files[0]))[0]
    outputs = write_chain_outputs(matched_data, out_dir, base_stem, verbose=verbose)

    # Log file with counts, command line, and full console transcript
    cmdline = " ".join(sys.argv)
    write_log(
        out_dir=out_dir,
        base_stem=base_stem,
        cmdline=cmdline,
        pep_count_after_psm_filter=pep_count_after_psm_filter,
        merged_data=merged_data,
        matched_data=matched_data,
        init_unique_psm=init_unique_psm,
        init_unique_pepseq=init_unique_pepseq,
        final_unique_psm=final_unique_psm,
        final_unique_pepseq=final_unique_pepseq,
    )

    if not outputs:
        log_print("No chain-specific TSVs were written (no matched rows).")


# -------------------------------
# CLI
# -------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Merge PSMs with BCRseq and compute peptide regional coverage."
    )
    parser.add_argument(
        "-psm_files", nargs="+", required=True,
        help="Proteome Discoverer PSM TXT/TSV file(s)."
    )
    parser.add_argument(
        "--bcrseq", dest="bcrseq_specs", nargs="+", metavar="CHAIN:FILE", required=True,
        help="Repeatable. One or more BCRseq TSV files, annotated with chain type, "
             "e.g., --bcrseq IGH:heavy.tsv IGK:kappa.tsv"
    )
    parser.add_argument(
        "--suffix", dest="suffix_specs", nargs="+", metavar="CHAIN:FILE", required=True,
        help="Repeatable. Map each chain type to an extension/C1 file (one AA sequence per line). "
             "The longest sequence per chain is appended when forming full_bcr_seq."
    )
    parser.add_argument(
        "--no-mixed-accessions", dest="allow_mixed_accessions",
        action="store_false",
        help="Require PSM rows to have ONLY BCRseq_ accessions. "
             "Default: mixed-accession PSM rows are allowed; non-BCRseq tokens are dropped before merging."
    )
    parser.add_argument(
        "--verbose", action="store_true",
        help="Write full merged columns (default writes a concise subset)."
    )

    args = parser.parse_args()
    main(
        psm_files=args.psm_files,
        bcrseq_specs=args.bcrseq_specs,
        suffix_specs=args.suffix_specs,
        allow_mixed_accessions=args.allow_mixed_accessions,  # default True; set False by --no-mixed-accessions
        verbose=args.verbose,
    )
