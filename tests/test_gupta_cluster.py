"""Tests for efficient, row-preserving BCR CDR3 clustering."""

from __future__ import annotations

import importlib.util
import io
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd


MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "workflows"
    / "bcrseq_transcript"
    / "gupta_cluster.py"
)
SPEC = importlib.util.spec_from_file_location("gupta_cluster", MODULE_PATH)
CLUSTERER = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = CLUSTERER
SPEC.loader.exec_module(CLUSTERER)


def annotations() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sequence_id": ["nt1", "nt2", "nt3", "nt4", "nt5", "nt6"],
            "v_call": ["V1*01", "V1*01", "V1*01", "V1*01", "V2*01", "V2*01"],
            "j_call": ["J1*01", "J1*01", "J1*01", "J1*01", "J1*01", "J1*01"],
            "cdr3_aa": ["AAA", "AAA", "AAB", "BBB", "AAA", "AAA"],
        }
    )


class GuptaClusterTests(unittest.TestCase):
    def test_clustering_unique_cdr3s_preserves_all_rows_and_membership(self) -> None:
        original = CLUSTERER.add_grouping_columns(annotations())
        clustered = CLUSTERER.perform_clustering(original, threshold=0, log_file=io.StringIO())

        self.assertEqual(len(clustered), len(original))
        self.assertEqual(
            clustered.loc[clustered["sequence_id"].isin(["nt1", "nt2"]), "ClusterID"].nunique(),
            1,
        )
        self.assertNotEqual(
            clustered.loc[clustered["sequence_id"] == "nt1", "ClusterID"].item(),
            clustered.loc[clustered["sequence_id"] == "nt3", "ClusterID"].item(),
        )
        self.assertNotEqual(
            clustered.loc[clustered["sequence_id"] == "nt1", "ClusterID"].item(),
            clustered.loc[clustered["sequence_id"] == "nt5", "ClusterID"].item(),
        )

    def test_auto_threshold_uses_unique_cdr3_sequences(self) -> None:
        frame = pd.DataFrame(
            {
                "v_call": ["V1*01"] * 7,
                "j_call": ["J1*01"] * 7,
                "cdr3_aa": ["AAA"] * 5 + ["AAB", "ACC"],
            }
        )
        unique = CLUSTERER.unique_cdr3_sequences(CLUSTERER.add_grouping_columns(frame))

        with tempfile.TemporaryDirectory() as temp_dir:
            with patch.object(CLUSTERER, "hamming_distance", wraps=CLUSTERER.hamming_distance) as distance:
                CLUSTERER.compute_auto_threshold(unique, temp_dir, io.StringIO())

        self.assertEqual(distance.call_count, 6)

    def test_pairwise_work_estimate_reflects_deduplication(self) -> None:
        grouped = CLUSTERER.add_grouping_columns(annotations())
        unique = CLUSTERER.unique_cdr3_sequences(grouped)

        self.assertEqual(CLUSTERER.pairwise_comparisons(grouped), 7)
        self.assertEqual(CLUSTERER.pairwise_comparisons(unique), 3)


if __name__ == "__main__":
    unittest.main()
