#!/usr/bin/env python3
"""
bcrseq_pipeline.py

Terminal-run orchestrator for a BCRseq pipeline across many samples stored as
subdirectories of a parent directory.

It can run any/all of these stages (in order):
1) trim_merge.py
2) identify_genes.py
3) filter_collapse.py
4) cluster.py
5) make_searchable.py

Inputs:
- parent directory containing sample subdirectories
- a config file specifying how to run each stage
- optional stage to start from
- optional flag to annotate filtered TSVs with isotype + tissue_type parsed
  from subdirectory names

Isotype & tissue parsing rules (only if --annotate_isotype_tissue is set):
- isotype: underscore-delimited token starting with "Ig" (e.g. "..._IgG_...")
- tissue:  underscore-delimited token starting with "tt" (e.g. "..._ttblood_...")
  tissue_type value is token after "tt" (e.g. "ttblood" -> "blood")

Example:
  mydata_IgG_ttblood/
  -> isotype="IgG", tissue_type="blood"

Behavior change when --annotate_isotype_tissue is set:
- isotype/tissue_type columns are added after filter_collapse (to the filtered TSV)
- all (annotated) filtered TSVs are concatenated into one combined dataset
- clustering is run once on the combined dataset (cluster.py by default)
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import subprocess
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

try:
    import yaml  # type: ignore
except Exception:
    yaml = None

try:
    import pandas as pd  # type: ignore
except Exception:
    pd = None

try:
    from tqdm import tqdm  # type: ignore
except Exception:
    tqdm = None


STAGES_ORDER = ["trim_merge", "identify_genes", "filter_collapse", "cluster", "make_searchable"]


@dataclass
class StageResult:
    stage: str
    sample_dir: Path
    ok: bool
    stdout_path: Path
    stderr_path: Path
    returncode: int
    command: str
    produced: Optional[Path] = None


def eprint(*args: object) -> None:
    print(*args, file=sys.stderr)


def now_stamp() -> str:
    return datetime.now().strftime("%Y-%m-%d_%H-%M-%S")


def load_config(config_path: Path) -> Dict[str, Any]:
    if not config_path.exists():
        raise FileNotFoundError(f"Config not found: {config_path}")

    text = config_path.read_text()
    suffix = config_path.suffix.lower()

    if suffix in [".json"]:
        return json.loads(text)

    if suffix in [".yml", ".yaml"]:
        if yaml is None:
            raise RuntimeError(
                "YAML config requested but PyYAML is not installed. "
                "Install pyyaml or use JSON."
            )
        return yaml.safe_load(text)

    # Try JSON first, then YAML
    try:
        return json.loads(text)
    except Exception:
        if yaml is None:
            raise RuntimeError(
                "Config is not valid JSON, and YAML parsing is unavailable. "
                "Install pyyaml or use JSON."
            )
        return yaml.safe_load(text)


def ensure_bool(x: Any, default: bool = False) -> bool:
    if isinstance(x, bool):
        return x
    if x is None:
        return default
    if isinstance(x, str):
        return x.strip().lower() in ["1", "true", "t", "yes", "y", "on"]
    if isinstance(x, (int, float)):
        return bool(x)
    return default


def stage_enabled(cfg: Dict[str, Any], stage: str) -> bool:
    stage_cfg = cfg.get("stages", {}).get(stage, {})
    return ensure_bool(stage_cfg.get("enabled", True), True)


def get_stage_args(cfg: Dict[str, Any], stage: str) -> List[str]:
    stage_cfg = cfg.get("stages", {}).get(stage, {})
    args = stage_cfg.get("args", [])
    if args is None:
        return []
    if not isinstance(args, list):
        raise ValueError(f'Config "stages.{stage}.args" must be a list of strings.')
    return [str(a) for a in args]


def get_stage_script(cfg: Dict[str, Any], stage: str) -> Path:
    scripts = cfg.get("scripts", {})
    p = scripts.get(stage)
    if not p:
        raise ValueError(f'Config missing scripts["{stage}"] path.')
    return Path(p).expanduser()


def list_sample_dirs(parent_dir: Path) -> List[Path]:
    return sorted([p for p in parent_dir.iterdir() if p.is_dir()])


def parse_isotype_tissue(sample_dir_name: str) -> Tuple[Optional[str], Optional[str]]:
    # underscore-delimited tokens
    tokens = [t for t in sample_dir_name.split("_") if t]
    isotype = None
    tissue = None

    for t in tokens:
        if isotype is None and t.startswith("Ig") and len(t) >= 3:
            isotype = t

    for t in tokens:
        if tissue is None and t.startswith("tt") and len(t) >= 3:
            tissue = t[2:]  # strip leading "tt"
    return isotype, tissue


def find_r1_r2(sample_dir: Path, patterns: Dict[str, str]) -> Tuple[Path, Path]:
    """
    Find R1 and R2 inputs for trim_merge stage.
    patterns:
      r1_regex: regex to match R1 file names
      r2_regex: regex to match R2 file names
    """
    r1_re = re.compile(patterns.get("r1_regex", r".*(_R1|R1).*\.fastq\.gz$",), re.IGNORECASE)
    r2_re = re.compile(patterns.get("r2_regex", r".*(_R2|R2).*\.fastq\.gz$",), re.IGNORECASE)

    files = [p for p in sample_dir.iterdir() if p.is_file()]
    r1s = [p for p in files if r1_re.match(p.name)]
    r2s = [p for p in files if r2_re.match(p.name)]

    if len(r1s) == 0 or len(r2s) == 0:
        raise FileNotFoundError(
            f"Could not find R1/R2 fastq.gz in {sample_dir} using regex "
            f"r1={r1_re.pattern} r2={r2_re.pattern}"
        )
    if len(r1s) > 1 or len(r2s) > 1:
        raise RuntimeError(
            f"Ambiguous R1/R2 in {sample_dir}. Found R1={len(r1s)} R2={len(r2s)}. "
            f"Adjust config r1_regex/r2_regex."
        )
    return r1s[0], r2s[0]


def find_single_input(sample_dir: Path, glob_pattern: str) -> Path:
    matches = sorted(sample_dir.glob(glob_pattern))
    if not matches:
        raise FileNotFoundError(f"No matches in {sample_dir} for glob: {glob_pattern}")
    if len(matches) > 1:
        raise RuntimeError(
            f"Ambiguous input in {sample_dir} for glob: {glob_pattern}. "
            f"Found {len(matches)} matches. Refine the glob."
        )
    return matches[0]


def run_command(
    cmd: List[str],
    workdir: Path,
    logs_dir: Path,
    stage: str,
    sample_name: str,
) -> Tuple[int, Path, Path, str]:
    logs_dir.mkdir(parents=True, exist_ok=True)
    stdout_path = logs_dir / f"{sample_name}.{stage}.stdout.txt"
    stderr_path = logs_dir / f"{sample_name}.{stage}.stderr.txt"
    cmd_str = " ".join([sh_quote(c) for c in cmd])

    with open(stdout_path, "w") as out_f, open(stderr_path, "w") as err_f:
        p = subprocess.run(cmd, cwd=str(workdir), stdout=out_f, stderr=err_f)
    return p.returncode, stdout_path, stderr_path, cmd_str


def sh_quote(s: str) -> str:
    # minimal shell quoting for logs, not used for execution
    if re.search(r"[^\w@%+=:,./-]", s):
        return "'" + s.replace("'", "'\"'\"'") + "'"
    return s


def annotate_tsv_with_isotype_tissue(tsv_path: Path, isotype: str, tissue: str) -> Path:
    """
    Add 'isotype' and 'tissue_type' columns to a TSV.
    Writes a new file next to the original with suffix _annotated.tsv.
    """
    if pd is None:
        raise RuntimeError("pandas is required for annotation. Install pandas or disable annotation.")

    df = pd.read_csv(tsv_path, sep="\t", dtype=str)
    df["isotype"] = isotype
    df["tissue_type"] = tissue
    out_path = tsv_path.with_name(tsv_path.stem + "_annotated.tsv")
    df.to_csv(out_path, sep="\t", index=False)
    return out_path


def concatenate_tsvs(tsv_paths: List[Path], out_path: Path) -> Path:
    if pd is None:
        raise RuntimeError("pandas is required for concatenation. Install pandas or disable concatenation.")
    if not tsv_paths:
        raise RuntimeError("No TSVs provided for concatenation.")

    frames = []
    for p in tsv_paths:
        frames.append(pd.read_csv(p, sep="\t", dtype=str))
    combined = pd.concat(frames, axis=0, ignore_index=True, sort=False)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(out_path, sep="\t", index=False)
    return out_path


def compute_stage_start_index(start_stage: str) -> int:
    if start_stage not in STAGES_ORDER:
        raise ValueError(f"Invalid --start_stage {start_stage}. Must be one of: {', '.join(STAGES_ORDER)}")
    return STAGES_ORDER.index(start_stage)


def describe_plan(cfg: Dict[str, Any], start_idx: int) -> List[str]:
    planned = []
    for stage in STAGES_ORDER[start_idx:]:
        if stage_enabled(cfg, stage):
            planned.append(stage)
    return planned


def require_tqdm_if_progress() -> None:
    if tqdm is None:
        raise RuntimeError("tqdm is required for progress display. Install with: pip install tqdm")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run a multi-stage BCRseq pipeline across subdirectories of a parent directory."
    )
    parser.add_argument("parent_dir", type=str, help="Parent directory containing sample subdirectories.")
    parser.add_argument("config", type=str, help="Config file path (.json, .yml, .yaml).")
    parser.add_argument(
        "--start_stage",
        default="trim_merge",
        choices=STAGES_ORDER,
        help="Which pipeline stage to start from.",
    )
    parser.add_argument(
        "--annotate_isotype_tissue",
        action="store_true",
        help="Annotate filtered TSV outputs with isotype and tissue_type parsed from subdirectory names.",
    )
    parser.add_argument(
        "--continue_on_error",
        action="store_true",
        help="If set, keep processing other samples/stages after an error.",
    )
    parser.add_argument(
        "--dry_run",
        action="store_true",
        help="Print commands that would run but do not execute.",
    )
    args = parser.parse_args()

    require_tqdm_if_progress()

    parent_dir = Path(args.parent_dir).expanduser()
    config_path = Path(args.config).expanduser()

    cfg = load_config(config_path)

    logs_root = Path(cfg.get("logs_dir", parent_dir / "pipeline_logs")).expanduser()
    logs_root.mkdir(parents=True, exist_ok=True)

    sample_dirs = list_sample_dirs(parent_dir)
    if not sample_dirs:
        raise RuntimeError(f"No sample subdirectories found in {parent_dir}")

    start_idx = compute_stage_start_index(args.start_stage)
    planned_stages = describe_plan(cfg, start_idx)

    # Per-stage input discovery rules
    input_rules = cfg.get("inputs", {})
    r1r2_patterns = input_rules.get("trim_merge", {}) if isinstance(input_rules.get("trim_merge", {}), dict) else {}

    assembled_glob = str(input_rules.get("identify_genes", {}).get("assembled_glob", "*.assembled.fastq"))
    annotated_glob = str(input_rules.get("filter_collapse", {}).get("annotated_tsv_glob", "*_IgBLAST.tsv"))
    filtered_glob = str(input_rules.get("cluster", {}).get("filtered_tsv_glob", "*_filtered.tsv"))
    clustered_glob = str(input_rules.get("make_searchable", {}).get("clustered_tsv_glob", "*_clustered.tsv"))
    extensions_txt = cfg.get("make_searchable_extensions_txt", None)

    all_results: List[StageResult] = []

    # Resolve scripts once
    scripts = {stage: get_stage_script(cfg, stage) for stage in STAGES_ORDER}

    # Progress units: per-sample stages up to filter_collapse run per sample, then global stages once
    per_sample_stages = []
    for stage in STAGES_ORDER[start_idx:]:
        if not stage_enabled(cfg, stage):
            continue
        if stage in ["trim_merge", "identify_genes", "filter_collapse"]:
            per_sample_stages.append(stage)

    global_stages = []
    for stage in STAGES_ORDER[start_idx:]:
        if not stage_enabled(cfg, stage):
            continue
        if stage in ["cluster", "make_searchable"]:
            global_stages.append(stage)

    total_units = len(sample_dirs) * len(per_sample_stages) + len(global_stages)
    total_units = max(1, total_units)
    pbar = tqdm(total=total_units, unit="stage", dynamic_ncols=True)

    produced_filtered_by_sample: Dict[str, Path] = {}

    try:
        # Per-sample phases (up to filter_collapse)
        for sample_dir in sample_dirs:
            sample_name = sample_dir.name
            sample_logs_dir = logs_root / sample_name
            sample_logs_dir.mkdir(parents=True, exist_ok=True)

            isotype = None
            tissue = None
            if args.annotate_isotype_tissue:
                isotype, tissue = parse_isotype_tissue(sample_name)
                if isotype is None or tissue is None:
                    raise RuntimeError(
                        f"--annotate_isotype_tissue was set but could not parse isotype/tissue from "
                        f"directory name: {sample_name}. Expected tokens like '_IgG_' and '_ttblood_'."
                    )

            produced_assembled: Optional[Path] = None
            produced_annotated: Optional[Path] = None
            produced_filtered: Optional[Path] = None

            for stage in STAGES_ORDER[start_idx:]:
                if stage not in ["trim_merge", "identify_genes", "filter_collapse"]:
                    continue
                if not stage_enabled(cfg, stage):
                    continue

                pbar.set_description(f"{sample_name}")
                pbar.set_postfix_str(stage)

                try:
                    if stage == "trim_merge":
                        r1, r2 = find_r1_r2(sample_dir, r1r2_patterns)
                        cmd = [sys.executable, str(scripts["trim_merge"]), str(r1), str(r2)] + get_stage_args(
                            cfg, "trim_merge"
                        )
                        if args.dry_run:
                            print(f"[DRY RUN] {sample_name} trim_merge: {' '.join(map(sh_quote, cmd))}")
                            rc = 0
                            stdout_path = sample_logs_dir / f"{sample_name}.trim_merge.stdout.txt"
                            stderr_path = sample_logs_dir / f"{sample_name}.trim_merge.stderr.txt"
                            cmd_str = " ".join(map(sh_quote, cmd))
                        else:
                            rc, stdout_path, stderr_path, cmd_str = run_command(
                                cmd=cmd,
                                workdir=sample_dir,
                                logs_dir=sample_logs_dir,
                                stage="trim_merge",
                                sample_name=sample_name,
                            )
                        if rc == 0:
                            matches = sorted(sample_dir.glob("*.assembled.fastq"))
                            produced_assembled = matches[0] if matches else None
                        all_results.append(
                            StageResult("trim_merge", sample_dir, rc == 0, stdout_path, stderr_path, rc, cmd_str, produced_assembled)
                        )
                        if rc != 0 and not args.continue_on_error:
                            raise RuntimeError(f"trim_merge failed for {sample_name} (rc={rc})")

                    elif stage == "identify_genes":
                        if produced_assembled is None:
                            produced_assembled = find_single_input(sample_dir, assembled_glob)
                        cmd = [sys.executable, str(scripts["identify_genes"]), str(produced_assembled)] + get_stage_args(
                            cfg, "identify_genes"
                        )
                        if args.dry_run:
                            print(f"[DRY RUN] {sample_name} identify_genes: {' '.join(map(sh_quote, cmd))}")
                            rc = 0
                            stdout_path = sample_logs_dir / f"{sample_name}.identify_genes.stdout.txt"
                            stderr_path = sample_logs_dir / f"{sample_name}.identify_genes.stderr.txt"
                            cmd_str = " ".join(map(sh_quote, cmd))
                        else:
                            rc, stdout_path, stderr_path, cmd_str = run_command(
                                cmd=cmd,
                                workdir=sample_dir,
                                logs_dir=sample_logs_dir,
                                stage="identify_genes",
                                sample_name=sample_name,
                            )
                        if rc == 0:
                            matches = sorted(sample_dir.glob("*_IgBLAST.tsv"))
                            produced_annotated = matches[0] if matches else None
                            if produced_annotated is None:
                                matches2 = sorted(sample_dir.glob("*MiXCR.tsv"))
                                produced_annotated = matches2[0] if matches2 else None
                        all_results.append(
                            StageResult(
                                "identify_genes",
                                sample_dir,
                                rc == 0,
                                stdout_path,
                                stderr_path,
                                rc,
                                cmd_str,
                                produced_annotated,
                            )
                        )
                        if rc != 0 and not args.continue_on_error:
                            raise RuntimeError(f"identify_genes failed for {sample_name} (rc={rc})")

                    elif stage == "filter_collapse":
                        if produced_annotated is None:
                            produced_annotated = find_single_input(sample_dir, annotated_glob)

                        cmd = [sys.executable, str(scripts["filter_collapse"]), str(produced_annotated)] + get_stage_args(
                            cfg, "filter_collapse"
                        )
                        if args.dry_run:
                            print(f"[DRY RUN] {sample_name} filter_collapse: {' '.join(map(sh_quote, cmd))}")
                            rc = 0
                            stdout_path = sample_logs_dir / f"{sample_name}.filter_collapse.stdout.txt"
                            stderr_path = sample_logs_dir / f"{sample_name}.filter_collapse.stderr.txt"
                            cmd_str = " ".join(map(sh_quote, cmd))
                        else:
                            rc, stdout_path, stderr_path, cmd_str = run_command(
                                cmd=cmd,
                                workdir=sample_dir,
                                logs_dir=sample_logs_dir,
                                stage="filter_collapse",
                                sample_name=sample_name,
                            )
                        if rc == 0:
                            produced_filtered = Path(str(produced_annotated).replace(".tsv", "_filtered.tsv"))
                            if not produced_filtered.exists():
                                matches = sorted(sample_dir.glob(filtered_glob))
                                produced_filtered = matches[0] if matches else None

                            if produced_filtered is not None and args.annotate_isotype_tissue:
                                produced_filtered = annotate_tsv_with_isotype_tissue(
                                    produced_filtered,
                                    isotype=isotype,  # type: ignore[arg-type]
                                    tissue=tissue,  # type: ignore[arg-type]
                                )

                            if produced_filtered is not None:
                                produced_filtered_by_sample[sample_name] = produced_filtered

                        all_results.append(
                            StageResult(
                                "filter_collapse",
                                sample_dir,
                                rc == 0,
                                stdout_path,
                                stderr_path,
                                rc,
                                cmd_str,
                                produced_filtered,
                            )
                        )
                        if rc != 0 and not args.continue_on_error:
                            raise RuntimeError(f"filter_collapse failed for {sample_name} (rc={rc})")

                    else:
                        raise RuntimeError(f"Unknown stage: {stage}")

                except Exception as ex:
                    eprint(f"[ERROR] {sample_name} {stage}: {ex}")
                    if not args.continue_on_error:
                        raise
                finally:
                    pbar.update(1)

        # Global phases (cluster, then make_searchable) run once on concatenated data
        combined_filtered: Optional[Path] = None
        combined_clustered: Optional[Path] = None

        if "cluster" in global_stages:
            stage = "cluster"
            sample_name = "ALL_SAMPLES"
            sample_dir = parent_dir
            sample_logs_dir = logs_root / sample_name
            sample_logs_dir.mkdir(parents=True, exist_ok=True)

            pbar.set_description(sample_name)
            pbar.set_postfix_str(stage)

            try:
                filtered_inputs = list(produced_filtered_by_sample.values())
                if not filtered_inputs:
                    raise RuntimeError("No filtered TSVs were produced; cannot run global clustering.")

                combined_filtered = logs_root / f"combined_filtered_{now_stamp()}.tsv"
                if args.dry_run:
                    print(
                        f"[DRY RUN] {sample_name} concatenate: "
                        f"{len(filtered_inputs)} TSVs -> {combined_filtered}"
                    )
                else:
                    concatenate_tsvs(filtered_inputs, combined_filtered)

                cluster_args = get_stage_args(cfg, "cluster")
                if not cluster_args:
                    raise ValueError(
                        'Config "stages.cluster.args" must include at least the clustering_executable path '
                        "(the second positional arg to cluster.py)."
                    )
                clustering_executable = cluster_args[0]
                extra_cluster_args = cluster_args[1:]

                cmd = [sys.executable, str(scripts["cluster"]), str(combined_filtered), str(clustering_executable)] + extra_cluster_args

                if args.dry_run:
                    print(f"[DRY RUN] {sample_name} cluster: {' '.join(map(sh_quote, cmd))}")
                    rc = 0
                    stdout_path = sample_logs_dir / f"{sample_name}.cluster.stdout.txt"
                    stderr_path = sample_logs_dir / f"{sample_name}.cluster.stderr.txt"
                    cmd_str = " ".join(map(sh_quote, cmd))
                else:
                    rc, stdout_path, stderr_path, cmd_str = run_command(
                        cmd=cmd,
                        workdir=logs_root,
                        logs_dir=sample_logs_dir,
                        stage="cluster",
                        sample_name=sample_name,
                    )

                if rc == 0:
                    # cluster.py writes <input>_clustered.tsv
                    combined_clustered = combined_filtered.with_suffix("").with_name(combined_filtered.stem + "_clustered.tsv")
                    if not combined_clustered.exists():
                        matches = sorted(logs_root.glob(combined_filtered.stem + "*_clustered.tsv"))
                        combined_clustered = matches[0] if matches else None

                all_results.append(
                    StageResult("cluster", sample_dir, rc == 0, stdout_path, stderr_path, rc, cmd_str, combined_clustered)
                )

                if rc != 0 and not args.continue_on_error:
                    raise RuntimeError(f"cluster failed (rc={rc})")

            except Exception as ex:
                eprint(f"[ERROR] {sample_name} {stage}: {ex}")
                if not args.continue_on_error:
                    raise
            finally:
                pbar.update(1)

        if "make_searchable" in global_stages:
            stage = "make_searchable"
            sample_name = "ALL_SAMPLES"
            sample_dir = parent_dir
            sample_logs_dir = logs_root / sample_name
            sample_logs_dir.mkdir(parents=True, exist_ok=True)

            pbar.set_description(sample_name)
            pbar.set_postfix_str(stage)

            try:
                if extensions_txt is None:
                    raise ValueError('Config must include "make_searchable_extensions_txt" path to the extensions .txt file.')

                ext_path = Path(extensions_txt).expanduser()
                if not ext_path.exists():
                    raise FileNotFoundError(f"Extensions txt not found: {ext_path}")

                if combined_clustered is None:
                    # Fall back to finding a clustered file in logs_root (useful if starting at make_searchable)
                    matches = sorted(logs_root.glob(clustered_glob))
                    combined_clustered = matches[0] if matches else None

                if combined_clustered is None:
                    raise RuntimeError("No combined clustered TSV found; cannot run make_searchable.")

                cmd = [sys.executable, str(scripts["make_searchable"]), str(combined_clustered), str(ext_path)] + get_stage_args(
                    cfg, "make_searchable"
                )

                if args.dry_run:
                    print(f"[DRY RUN] {sample_name} make_searchable: {' '.join(map(sh_quote, cmd))}")
                    rc = 0
                    stdout_path = sample_logs_dir / f"{sample_name}.make_searchable.stdout.txt"
                    stderr_path = sample_logs_dir / f"{sample_name}.make_searchable.stderr.txt"
                    cmd_str = " ".join(map(sh_quote, cmd))
                else:
                    rc, stdout_path, stderr_path, cmd_str = run_command(
                        cmd=cmd,
                        workdir=logs_root,
                        logs_dir=sample_logs_dir,
                        stage="make_searchable",
                        sample_name=sample_name,
                    )

                all_results.append(StageResult("make_searchable", sample_dir, rc == 0, stdout_path, stderr_path, rc, cmd_str, None))

                if rc != 0 and not args.continue_on_error:
                    raise RuntimeError(f"make_searchable failed (rc={rc})")

            except Exception as ex:
                eprint(f"[ERROR] {sample_name} {stage}: {ex}")
                if not args.continue_on_error:
                    raise
            finally:
                pbar.update(1)

    finally:
        pbar.close()

    # Summary
    ok = sum(1 for r in all_results if r.ok)
    bad = sum(1 for r in all_results if not r.ok)

    summary_path = logs_root / f"summary_{now_stamp()}.tsv"
    with open(summary_path, "w") as f:
        f.write("\t".join(["sample", "stage", "ok", "returncode", "produced", "stdout", "stderr", "command"]) + "\n")
        for r in all_results:
            f.write(
                "\t".join(
                    [
                        r.sample_dir.name,
                        r.stage,
                        "1" if r.ok else "0",
                        str(r.returncode),
                        str(r.produced) if r.produced else "",
                        str(r.stdout_path),
                        str(r.stderr_path),
                        r.command,
                    ]
                )
                + "\n"
            )

    print(f"Done. Stage runs ok={ok} failed={bad}")
    print(f"Summary: {summary_path}")
    if bad > 0:
        sys.exit(2)


if __name__ == "__main__":
    main()