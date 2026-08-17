"""Tests for CDR-derived lineage abundance comparisons."""

from __future__ import annotations

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd


SCRIPT_PATH = Path(__file__).parents[1] / "workflows" / "igseq_proteomics" / "compare_cdr_lineage_abundance.py"
SPEC = importlib.util.spec_from_file_location("compare_cdr_lineage_abundance", SCRIPT_PATH)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def mapped_row(
    cluster: str,
    peptide: str,
    start: int,
    end: int,
    abundance: float = 10.0,
    match_found: str = "True",
) -> dict[str, object]:
    return {
        "ClusterID": cluster,
        "peptide_sequence": peptide,
        "peptide_abundance_total": abundance,
        "match_found": match_found,
        "match_start": start,
        "match_end": end,
        "fwr1_aa": "AAAAA",
        "cdr1_aa": "BBB",
        "fwr2_aa": "CCCCC",
        "cdr2_aa": "DDD",
        "fwr3_aa": "EEEEE",
        "cdr3_aa": "FFF",
        "fwr4_aa": "GG",
    }


class CompareCdrLineageAbundanceTests(unittest.TestCase):
    def test_coordinates_and_minimum_overlap_classification(self) -> None:
        row = pd.Series(mapped_row("1", "pep", 5, 8))
        coordinates = MODULE.region_coordinates(row)
        self.assertEqual(coordinates["cdr1_aa"], (5, 8))
        self.assertEqual(coordinates["cdr2_aa"], (13, 16))
        self.assertEqual(MODULE.overlap_length(4, 8, 5, 8), 3)
        classified = MODULE.classify_mapped_rows(pd.DataFrame([mapped_row("1", "pep", 4, 8)]), 3)
        self.assertTrue(classified.iloc[0]["cdr1_aa_qualifies"])
        classified = MODULE.classify_mapped_rows(pd.DataFrame([mapped_row("1", "pep", 4, 7)]), 3)
        self.assertFalse(classified.iloc[0]["cdr1_aa_qualifies"])

    def test_selection_classes_and_combined_one_count(self) -> None:
        rows = [
            mapped_row("1", "cdr1", 5, 8, 10),
            mapped_row("1", "cdr2", 13, 16, 20),
            mapped_row("1", "both", 5, 16, 30),
            mapped_row("1", "cdr3", 21, 24, 40),
        ]
        collapsed = MODULE.collapse_quantitative_evidence(MODULE.classify_mapped_rows(pd.DataFrame(rows), 3), "peptide_abundance_total")
        classes = collapsed.set_index("peptide_sequence")
        self.assertTrue(classes.loc["cdr1", "cdr1"])
        self.assertTrue(classes.loc["cdr2", "cdr2"])
        self.assertTrue(classes.loc["both", "cdr1_cdr2"])
        self.assertFalse(classes.loc["both", "cdr1"])
        self.assertFalse(classes.loc["both", "cdr2"])
        self.assertTrue(classes.loc["cdr3", "cdr3"])
        abundance = MODULE.lineage_abundance_table(collapsed).iloc[0]
        self.assertEqual(abundance["cdr1_abundance"], 10)
        self.assertEqual(abundance["cdr2_abundance"], 20)
        self.assertEqual(abundance["cdr1_cdr2_abundance"], 60)

    def test_multiple_bcr_members_do_not_double_count(self) -> None:
        rows = [mapped_row("1", "shared", 5, 8, 25), mapped_row("1", "shared", 5, 8, 25)]
        collapsed = MODULE.collapse_quantitative_evidence(MODULE.classify_mapped_rows(pd.DataFrame(rows), 3), "peptide_abundance_total")
        self.assertEqual(len(collapsed), 1)
        self.assertEqual(collapsed.iloc[0]["abundance"], 25)

    def test_comparisons_top_membership_and_sparse_handling(self) -> None:
        rows = [
            mapped_row("1", "one_cdr1", 5, 8, 10),
            mapped_row("1", "one_cdr3", 21, 24, 20),
            mapped_row("2", "two_cdr1", 5, 8, 20),
            mapped_row("2", "two_cdr3", 21, 24, 10),
            mapped_row("3", "three_cdr1", 5, 8, 5),
        ]
        collapsed = MODULE.collapse_quantitative_evidence(MODULE.classify_mapped_rows(pd.DataFrame(rows), 3), "peptide_abundance_total")
        lineages = MODULE.lineage_abundance_table(collapsed)
        pair, membership, metrics = MODULE.comparison_tables(lineages, "cdr1", "cdr3", 2)
        self.assertEqual(metrics["shared_detection_count"], 2)
        self.assertEqual(metrics["top_n_overlap_count"], 2)
        self.assertEqual(metrics["top_n_jaccard_similarity"], 1.0)
        self.assertEqual(set(pair["detection_status"]), {"both", "left_only"})
        self.assertTrue(membership["in_intersection"].all())
        _, _, sparse = MODULE.comparison_tables(lineages, "cdr2", "cdr3", 2)
        self.assertEqual(sparse["shared_detection_count"], 0)
        self.assertIsNone(sparse["spearman_rho"])

    def test_run_writes_seven_comparisons_and_uses_abundance_fallback(self) -> None:
        rows = [
            mapped_row("1", "one_cdr1", 5, 8, 10),
            mapped_row("1", "one_cdr3", 21, 24, 20),
            mapped_row("2", "two_cdr1", 5, 8, 20),
            mapped_row("2", "two_cdr3", 21, 24, 10),
        ]
        for row in rows:
            row["peptide_abundance"] = row.pop("peptide_abundance_total")
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            input_path = root / "mapped.tsv"
            pd.DataFrame(rows).to_csv(input_path, sep="\t", index=False)
            outputs = MODULE.run(input_path, root / "results", 2, 3)
            summary = pd.read_csv(outputs["comparison_summary"], sep="\t")
            self.assertEqual(len(summary), 7)
            self.assertTrue(outputs["cdr1_vs_cdr3_scatter"].exists())
            self.assertFalse((root / "results" / "mapped_cdr2_vs_cdr3_scatter.png").exists())
            self.assertEqual(len(list((root / "results").glob("*_lineages.tsv"))), 7)
            self.assertEqual(len(list((root / "results").glob("*_membership.tsv"))), 7)
            log = json.loads(outputs["run_log"].read_text())
            self.assertEqual(log["outputs"]["run_log"], str(outputs["run_log"]))


if __name__ == "__main__":
    unittest.main()
