"""
Unit Tests for Smith-Waterman Local Alignment Algorithm.
"""

import sys
import os
import unittest

# Ensure project root is in path
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, ".."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from src.algorithms.smith_waterman import smith_waterman, format_matrix


class TestSmithWaterman(unittest.TestCase):
    """Test suite for Smith-Waterman local dynamic programming alignment."""

    def test_conserved_local_motif(self):
        """Verify algorithm isolates conserved motif embedded in unrelated flanking sequences."""
        # 'HELLOWORLD' embedded in both sequences
        seq1 = "XXXXXXHELLOWORLDYYYYYY"
        seq2 = "ZZZZZZHELLOWORLDWWWWWW"

        res = smith_waterman(seq1, seq2, match_score=3, mismatch_penalty=-2, gap_penalty=-3)

        self.assertEqual(res["aligned_seq1"], "HELLOWORLD")
        self.assertEqual(res["aligned_seq2"], "HELLOWORLD")
        self.assertEqual(res["score"], len("HELLOWORLD") * 3)
        self.assertEqual(res["identity"], 1.0)
        self.assertEqual(res["start1"], 7)
        self.assertEqual(res["end1"], 16)
        self.assertEqual(res["start2"], 7)
        self.assertEqual(res["end2"], 16)

    def test_verified_sample_acacacta(self):
        """Verify alignment on Azheem's sample: ACACACTA vs AGCACACA."""
        seq1 = "ACACACTA"
        seq2 = "AGCACACA"
        res = smith_waterman(seq1, seq2, match_score=3, mismatch_penalty=-1, gap_penalty=-2)

        self.assertGreater(res["score"], 0)
        self.assertEqual(len(res["aligned_seq1"]), len(res["aligned_seq2"]))
        # Gapped aligned segment must match substring
        ungapped1 = res["aligned_seq1"].replace("-", "")
        self.assertIn(ungapped1, seq1)

    def test_zero_floor_boundary(self):
        """Verify matrix values never drop below zero."""
        seq1 = "TTTT"
        seq2 = "AAAA"
        res = smith_waterman(seq1, seq2, match_score=2, mismatch_penalty=-5, gap_penalty=-5)

        for row in res["matrix"]:
            for val in row:
                self.assertGreaterEqual(val, 0)
        self.assertEqual(res["score"], 0)
        self.assertEqual(res["aligned_seq1"], "")

    def test_identical_sequences(self):
        """Identical sequences match in their entirety."""
        seq = "CGCATAG"
        res = smith_waterman(seq, seq, match_score=2, mismatch_penalty=-1, gap_penalty=-2)

        self.assertEqual(res["aligned_seq1"], seq)
        self.assertEqual(res["aligned_seq2"], seq)
        self.assertEqual(res["score"], len(seq) * 2)
        self.assertEqual(res["start1"], 1)
        self.assertEqual(res["end1"], len(seq))

    def test_empty_sequences(self):
        """Empty sequence edge cases."""
        res_both = smith_waterman("", "")
        self.assertEqual(res_both["score"], 0)
        self.assertEqual(res_both["aligned_seq1"], "")

        res_one = smith_waterman("ACGT", "")
        self.assertEqual(res_one["score"], 0)
        self.assertEqual(res_one["aligned_seq1"], "")

    def test_custom_penalties(self):
        """Custom match and mismatch weights affect motif extension."""
        seq1 = "AAACGTTT"
        seq2 = "TTTCGTAAA"

        # Match 'CGT'
        res = smith_waterman(seq1, seq2, match_score=4, mismatch_penalty=-2, gap_penalty=-2)
        self.assertEqual(res["aligned_seq1"], "CGT")
        self.assertEqual(res["aligned_seq2"], "CGT")
        self.assertEqual(res["score"], 12)

    def test_matrix_formatting(self):
        """Matrix formatting outputs expected rows."""
        seq1 = "ACG"
        seq2 = "ACT"
        res = smith_waterman(seq1, seq2)
        out = format_matrix(res["matrix"], seq1, seq2)
        self.assertIn("A", out)
        self.assertIn("C", out)


if __name__ == "__main__":
    unittest.main()
