"""
Smith-Waterman Local Sequence Alignment Algorithm
================================================

A pure Python implementation of the Smith-Waterman dynamic programming
algorithm for identifying optimal local regions of similarity between two
biological sequences (DNA, RNA, or Protein).

Features:
- Dynamic programming matrix with zero-floor boundary condition (no negative propagation).
- Dynamic tracking of global maximum coordinate and score.
- Traceback terminating at positive boundaries to isolate conserved local motifs.
- Customizable scoring parameters: match reward, mismatch penalty, and gap penalty.
- Coordinate tracking: returns 1-based start and end indices of the local alignment.
- Handles edge cases: empty sequences, complete divergence (score = 0), and variable lengths.

Author: Azheem Syahmee
License: MIT
"""

from typing import List, Tuple, Dict, Any, Optional


def smith_waterman(
    seq1: str,
    seq2: str,
    match_score: int = 3,
    mismatch_penalty: int = -1,
    gap_penalty: int = -2,
) -> Dict[str, Any]:
    """
    Execute Smith-Waterman local alignment on two biological sequences.

    Args:
        seq1: First sequence string (target/reference).
        seq2: Second sequence string (query).
        match_score: Reward score for identity match (default: 3).
        mismatch_penalty: Penalty score for mismatch (default: -1).
        gap_penalty: Penalty score for gap insertion/deletion (default: -2).

    Returns:
        Dictionary containing:
        - 'aligned_seq1': Best locally aligned segment from seq1 with gaps.
        - 'aligned_seq2': Best locally aligned segment from seq2 with gaps.
        - 'score': Maximum local alignment score.
        - 'matrix': Full dynamic programming scoring matrix.
        - 'start1', 'end1': 1-based start and end coordinates in seq1.
        - 'start2', 'end2': 1-based start and end coordinates in seq2.
        - 'identity': Fraction of identical aligned positions in the motif.
        - 'similarity_str': Match visualization string ('|' match, '.' mismatch, ' ' gap).
    """
    seq1_clean = seq1.strip().upper()
    seq2_clean = seq2.strip().upper()

    if not seq1_clean or not seq2_clean:
        return {
            "aligned_seq1": "",
            "aligned_seq2": "",
            "score": 0,
            "matrix": [[0]],
            "start1": 0,
            "end1": 0,
            "start2": 0,
            "end2": 0,
            "identity": 0.0,
            "similarity_str": "",
        }

    n = len(seq1_clean) + 1
    m = len(seq2_clean) + 1
    score_matrix = [[0] * m for _ in range(n)]
    max_score = 0
    max_position = (0, 0)

    # Matrix filling with zero-floor local alignment recurrence
    for i in range(1, n):
        for j in range(1, m):
            is_match = seq1_clean[i - 1] == seq2_clean[j - 1]
            diag = score_matrix[i - 1][j - 1] + (match_score if is_match else mismatch_penalty)
            del_val = score_matrix[i - 1][j] + gap_penalty
            ins_val = score_matrix[i][j - 1] + gap_penalty

            cell_score = max(0, diag, del_val, ins_val)
            score_matrix[i][j] = cell_score

            if cell_score > max_score:
                max_score = cell_score
                max_position = (i, j)

    if max_score == 0:
        return {
            "aligned_seq1": "",
            "aligned_seq2": "",
            "score": 0,
            "matrix": score_matrix,
            "start1": 0,
            "end1": 0,
            "start2": 0,
            "end2": 0,
            "identity": 0.0,
            "similarity_str": "",
        }

    # Traceback starting from global maximum until reaching a zero cell
    align1 = []
    align2 = []
    i, j = max_position
    end1, end2 = i, j

    while i > 0 and j > 0 and score_matrix[i][j] > 0:
        current = score_matrix[i][j]
        diag = score_matrix[i - 1][j - 1]
        up = score_matrix[i - 1][j]
        left = score_matrix[i][j - 1]

        is_match = seq1_clean[i - 1] == seq2_clean[j - 1]
        expected_diag = diag + (match_score if is_match else mismatch_penalty)

        if current == expected_diag:
            align1.append(seq1_clean[i - 1])
            align2.append(seq2_clean[j - 1])
            i -= 1
            j -= 1
        elif current == up + gap_penalty:
            align1.append(seq1_clean[i - 1])
            align2.append("-")
            i -= 1
        elif current == left + gap_penalty:
            align1.append("-")
            align2.append(seq2_clean[j - 1])
            j -= 1
        else:
            break

    start1 = i + 1
    start2 = j + 1

    aligned_seq1 = "".join(reversed(align1))
    aligned_seq2 = "".join(reversed(align2))

    # Calculate identity and match symbols
    matches = 0
    sim_chars = []
    total_len = len(aligned_seq1)

    for c1, c2 in zip(aligned_seq1, aligned_seq2):
        if c1 == c2 and c1 != "-":
            matches += 1
            sim_chars.append("|")
        elif c1 == "-" or c2 == "-":
            sim_chars.append(" ")
        else:
            sim_chars.append(".")

    identity = (matches / total_len) if total_len > 0 else 0.0
    similarity_str = "".join(sim_chars)

    return {
        "aligned_seq1": aligned_seq1,
        "aligned_seq2": aligned_seq2,
        "score": max_score,
        "matrix": score_matrix,
        "start1": start1,
        "end1": end1,
        "start2": start2,
        "end2": end2,
        "identity": round(identity, 4),
        "similarity_str": similarity_str,
    }


def format_matrix(matrix: List[List[int]], seq1: str, seq2: str) -> str:
    """
    Format the dynamic programming matrix as a human-readable ASCII table.

    Args:
        matrix: 2D scoring matrix.
        seq1: First sequence (row labels).
        seq2: Second sequence (column labels).

    Returns:
        Multi-line formatted table string.
    """
    seq1_clean = seq1.strip().upper()
    seq2_clean = seq2.strip().upper()

    col_headers = [" ", " "] + list(seq2_clean)
    header_line = "  ".join(f"{ch:>4}" for ch in col_headers)
    lines = [header_line, "-" * len(header_line)]

    row_labels = [" "] + list(seq1_clean)
    for i, row in enumerate(matrix):
        label = row_labels[i] if i < len(row_labels) else " "
        row_str = f"{label:>4}  " + "  ".join(f"{val:>4}" for val in row)
        lines.append(row_str)

    return "\n".join(lines)


if __name__ == "__main__":
    # Test sample local alignment with embedded conserved motif
    sample_seq1 = "ACACACTA"
    sample_seq2 = "AGCACACA"

    res = smith_waterman(sample_seq1, sample_seq2, match_score=3, mismatch_penalty=-1, gap_penalty=-2)
    print("=== Smith-Waterman Local Alignment ===")
    print("Sequence 1:", sample_seq1)
    print("Sequence 2:", sample_seq2)
    print(f"Motif Span : Seq1 [{res['start1']}..{res['end1']}], Seq2 [{res['start2']}..{res['end2']}]")
    print("Aligned 1 :", res["aligned_seq1"])
    print("Match     :", res["similarity_str"])
    print("Aligned 2 :", res["aligned_seq2"])
    print(f"Max Score : {res['score']}")
    print(f"Identity  : {res['identity'] * 100:.1f}%\n")
    print("Dynamic Programming Matrix:")
    print(format_matrix(res["matrix"], sample_seq1, sample_seq2))
