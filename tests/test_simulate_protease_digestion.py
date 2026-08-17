"""Tests for AIRR/IgBLAST protease-digestion simulation."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
import tempfile
import unittest

import pandas as pd


MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "workflows"
    / "igseq_proteomics"
    / "simulate_protease_digestion.py"
)
SPEC = importlib.util.spec_from_file_location("simulate_protease_digestion", MODULE_PATH)
SIMULATOR = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = SIMULATOR
SPEC.loader.exec_module(SIMULATOR)


def rule(**positions: str) -> object:
    """Create a rule with unspecified positions left unrestricted."""
    return SIMULATOR.ProteaseRule(
        {position: positions.get(position, "") for position in ("P4", "P3", "P2", "P1", "P1'", "P2'")}
    )


class ProteaseDigestionTests(unittest.TestCase):
    """Verify cleavage semantics and AIRR input validation."""

    def test_p1_after_cleavage_and_no_duplicate_residue(self) -> None:
        boundaries = SIMULATOR.find_cleavage_boundaries("AKCD", [rule(P1="K")])
        self.assertEqual(boundaries, [2])
        self.assertEqual(SIMULATOR.peptides_from_boundaries("AKCD", boundaries), ["AK", "CD"])

    def test_p1_prime_cleavage_and_terminal_cleavage(self) -> None:
        boundaries = SIMULATOR.find_cleavage_boundaries("ADK", [rule(**{"P1'": "D"}), rule(P1="K")])
        self.assertEqual(boundaries, [1, 3])
        self.assertEqual(SIMULATOR.peptides_from_boundaries("ADK", boundaries), ["A", "DK"])

    def test_duplicate_matching_rules_produce_one_boundary(self) -> None:
        boundaries = SIMULATOR.find_cleavage_boundaries("AKC", [rule(P1="K"), rule(P1="K")])
        self.assertEqual(boundaries, [2])

    def test_cdr3_metrics_and_weighted_summary(self) -> None:
        row = pd.Series(
            {
                "fwr3_aa": "AAAK",
                "cdr3_aa": "CDR",
                "fwr4_aa": "RAAA",
                "nt_seq_count": 5,
            },
            name=0,
        )
        detail = SIMULATOR.digest_record(row, "Trypsin", [rule(P1="K"), rule(P1="R")])
        self.assertFalse(detail["cuts_inside_cdr3"])
        self.assertEqual(detail["cdr3_spanning_peptide"], "CDR")
        self.assertEqual(detail["cdr3_spanning_peptide_length"], 3)
        self.assertEqual(detail["upstream_cut_distance_aa"], 0)
        self.assertEqual(detail["downstream_cut_distance_aa"], 0)

        summary = SIMULATOR.create_summary(pd.DataFrame([detail]))
        self.assertEqual(summary.loc[0, "weighted_sequence_count"], 5)
        self.assertEqual(summary.loc[0, "weighted_cdr3_spanning_count"], 5)

    def test_rejects_missing_columns_invalid_counts_and_empty_regions(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "input.tsv"
            pd.DataFrame({"fwr3_aa": ["AAA"]}).to_csv(path, sep="\t", index=False)
            with self.assertRaisesRegex(ValueError, "missing required"):
                SIMULATOR.read_annotation_table(path)

            pd.DataFrame(
                {
                    "fwr3_aa": ["AAA", "AAA"],
                    "cdr3_aa": ["CAR", "CAR"],
                    "fwr4_aa": ["AAA", "AAA"],
                    "nt_seq_count": ["0", "1"],
                }
            ).to_csv(path, sep="\t", index=False)
            with self.assertRaisesRegex(ValueError, "positive integer"):
                SIMULATOR.read_annotation_table(path)

            pd.DataFrame(
                {
                    "fwr3_aa": ["AAA"],
                    "cdr3_aa": [""],
                    "fwr4_aa": ["AAA"],
                    "nt_seq_count": ["1"],
                }
            ).to_csv(path, sep="\t", index=False)
            with self.assertRaisesRegex(ValueError, "cdr3_aa is empty"):
                SIMULATOR.read_annotation_table(path)


if __name__ == "__main__":
    unittest.main()
