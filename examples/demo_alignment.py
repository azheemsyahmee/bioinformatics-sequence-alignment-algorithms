"""
Command-Line Alignment Demonstration Script
===========================================

Executes both Needleman-Wunsch global alignment and Smith-Waterman local
alignment on standard benchmark sequence pairs, printing formatted alignment
outputs and dynamic programming matrices directly to the console.

Usage:
    python examples/demo_alignment.py
"""

import sys
import os

# Ensure package root is in path
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, ".."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from src.algorithms.needleman_wunsch import needleman_wunsch, format_matrix as format_nw
from src.algorithms.smith_waterman import smith_waterman, format_matrix as format_sw


def run_demo():
    print("=" * 75)
    print(" BIOINFORMATICS SEQUENCE ALIGNMENT SUITE - CLI DEMONSTRATION")
    print("=" * 75)

    # 1. NEEDLEMAN-WUNSCH GLOBAL ALIGNMENT DEMO
    print("\n[1] DEMO: Needleman-Wunsch Global Alignment (End-to-End)")
    print("-" * 75)
    seq1_nw = "FKHMEDPLE"
    seq2_nw = "FMDTPLNE"
    match_nw = 1
    mismatch_nw = -2
    gap_nw = -2

    print(f"Input Sequence 1 : {seq1_nw}")
    print(f"Input Sequence 2 : {seq2_nw}")
    print(f"Scoring Scheme   : Match = +{match_nw}, Mismatch = {mismatch_nw}, Gap = {gap_nw}\n")

    res_nw = needleman_wunsch(seq1_nw, seq2_nw, match_score=match_nw, mismatch_penalty=mismatch_nw, gap_penalty=gap_nw)

    print("Optimal Global Alignment:")
    print(f"  Seq1 : {res_nw['aligned_seq1']}")
    print(f"         {res_nw['similarity_str']}")
    print(f"  Seq2 : {res_nw['aligned_seq2']}")
    print(f"Total Score : {res_nw['score']}")
    print(f"Identity    : {res_nw['identity'] * 100:.1f}%\n")
    print("Full Dynamic Programming Matrix:")
    print(format_nw(res_nw["matrix"], seq1_nw, seq2_nw))

    # 2. SMITH-WATERMAN LOCAL ALIGNMENT DEMO
    print("\n" + "=" * 75)
    print("[2] DEMO: Smith-Waterman Local Alignment (Conserved Motif Search)")
    print("-" * 75)
    seq1_sw = "TTTAAACACACTAGGGCCC"
    seq2_sw = "CCCGGAGCACACATTTTTT"
    match_sw = 3
    mismatch_sw = -1
    gap_sw = -2

    print(f"Input Sequence 1 (Target) : {seq1_sw}")
    print(f"Input Sequence 2 (Query)  : {seq2_sw}")
    print(f"Scoring Scheme            : Match = +{match_sw}, Mismatch = {mismatch_sw}, Gap = {gap_sw}\n")

    res_sw = smith_waterman(seq1_sw, seq2_sw, match_score=match_sw, mismatch_penalty=mismatch_sw, gap_penalty=gap_sw)

    print("Optimal Local Conserved Motif:")
    print(f"  Seq1 [{res_sw['start1']}..{res_sw['end1']}] : {res_sw['aligned_seq1']}")
    print(f"                 {res_sw['similarity_str']}")
    print(f"  Seq2 [{res_sw['start2']}..{res_sw['end2']}] : {res_sw['aligned_seq2']}")
    print(f"Maximum Local Score : {res_sw['score']}")
    print(f"Motif Identity      : {res_sw['identity'] * 100:.1f}%")
    print("=" * 75)


if __name__ == "__main__":
    run_demo()
