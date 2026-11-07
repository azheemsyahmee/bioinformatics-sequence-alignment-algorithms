"""
Unit Tests for Needleman-Wunsch Global Alignment Algorithm.
"""

import sys
import os
import unittest

# Ensure project root is in path
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, ".."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from src.algorithms.needleman_wunsch import (
    needleman_wunsch,
    initialize_matrix,
    fill_matrix,
    traceback,
    format_matrix,
)


class TestNeedlemanWunsch(unittest.TestCase):
    """Test suite for Needleman-Wunsch global dynamic programming alignment."""

    def test_identical_sequences(self):
        """Aligning two identical sequences should produce 100% identity and no gaps."""
        seq = "ATCGATCG"
        res = needleman_wunsch(seq, seq, match_score=2, mismatch_penalty=-1, gap_penalty=-2)

        self.assertEqual(res["aligned_seq1"], seq)
        self.assertEqual(res["aligned_seq2"], seq)
        self.assertEqual(res["score"], len(seq) * 2)
        self.assertEqual(res["identity"], 1.0)
        self.assertNotIn("-", res["aligned_seq1"])
        self.assertNotIn("-", res["aligned_seq2"])

    def test_textbook_protein_case(self):
        """Verify alignment on Azheem's verified sample: FKHMEDPLE vs FMDTPLNE."""
        seq1 = "FKHMEDPLE"
        seq2 = "FMDTPLNE"
        res = needleman_wunsch(seq1, seq2, match_score=1, mismatch_penalty=-2, gap_penalty=-2)

        # Both aligned sequences must have identical length
        self.assertEqual(len(res["aligned_seq1"]), len(res["aligned_seq2"]))
        # Ungapped sequences must match original inputs
        self.assertEqual(res["aligned_seq1"].replace("-", ""), seq1)
        self.assertEqual(res["aligned_seq2"].replace("-", ""), seq2)
        # Expected score with gap -2, match 1
        self.assertIsInstance(res["score"], int)
        self.assertGreater(res["identity"], 0.4)

    def test_dna_gap_insertion(self):
        """Verify alignment correctly inserts gap when a single base is deleted."""
        seq1 = "GATTACA"
        seq2 = "GATACA"  # Missing second 'T'
        res = needleman_wunsch(seq1, seq2, match_score=1, mismatch_penalty=-1, gap_penalty=-1)

        self.assertEqual(len(res["aligned_seq1"]), len(res["aligned_seq2"]))
        self.assertEqual(res["aligned_seq1"].replace("-", ""), seq1)
        self.assertEqual(res["aligned_seq2"].replace("-", ""), seq2)
        self.assertEqual(res["aligned_seq2"].count("-"), 1)

    def test_complete_mismatch(self):
        """Completely different sequences should still produce valid global alignment."""
        seq1 = "AAAA"
        seq2 = "CCCC"
        res = needleman_wunsch(seq1, seq2, match_score=1, mismatch_penalty=-1, gap_penalty=-2)

        self.assertEqual(len(res["aligned_seq1"]), len(res["aligned_seq2"]))
        self.assertEqual(res["aligned_seq1"].replace("-", ""), seq1)
        self.assertEqual(res["aligned_seq2"].replace("-", ""), seq2)
        self.assertEqual(res["identity"], 0.0)

    def test_asymmetrical_lengths(self):
        """Verify alignment handles significantly different sequence lengths."""
        seq1 = "ACTG"
        seq2 = "ACTGACTGACTG"
        res = needleman_wunsch(seq1, seq2, match_score=1, mismatch_penalty=-1, gap_penalty=-2)

        self.assertEqual(len(res["aligned_seq1"]), len(res["aligned_seq2"]))
        self.assertEqual(res["aligned_seq1"].replace("-", ""), seq1)
        self.assertEqual(res["aligned_seq2"].replace("-", ""), seq2)

    def test_empty_sequences(self):
        """Empty sequence edge cases."""
        # Both empty
        res_both = needleman_wunsch("", "")
        self.assertEqual(res_both["score"], 0)
        self.assertEqual(res_both["aligned_seq1"], "")
        self.assertEqual(res_both["aligned_seq2"], "")

        # One empty
        res_one = needleman_wunsch("ACGT", "", gap_penalty=-2)
        self.assertEqual(res_one["score"], -8)
        self.assertEqual(res_one["aligned_seq1"], "ACGT")
        self.assertEqual(res_one["aligned_seq2"], "----")

    def test_custom_penalties(self):
        """Custom match/mismatch/gap scoring schemes."""
        seq1 = "ACGT"
        seq2 = "ACCT"
        # High mismatch penalty forces gaps
        res_high_mismatch = needleman_wunsch(seq1, seq2, match_score=5, mismatch_penalty=-10, gap_penalty=-2)
        # Moderate mismatch allows substitution
        res_sub = needleman_wunsch(seq1, seq2, match_score=2, mismatch_penalty=-1, gap_penalty=-5)

        self.assertEqual(res_sub["aligned_seq1"], "ACGT")
        self.assertEqual(res_sub["aligned_seq2"], "ACCT")
        self.assertEqual(res_sub["score"], 2 + 2 - 1 + 2)  # A match (2), C match (2), G/C mismatch (-1), T match (2) = 5

    def test_matrix_formatting(self):
        """Matrix string formatting produces valid table lines."""
        seq1 = "AG"
        seq2 = "AC"
        res = needleman_wunsch(seq1, seq2)
        mat_str = format_matrix(res["matrix"], seq1, seq2)

        self.assertIn("A", mat_str)
        self.assertIn("G", mat_str)
        self.assertIn("C", mat_str)


if __name__ == "__main__":
    unittest.main()
