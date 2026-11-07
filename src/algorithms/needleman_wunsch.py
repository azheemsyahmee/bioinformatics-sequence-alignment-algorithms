"""
Needleman-Wunsch Global Sequence Alignment Algorithm
===================================================

A pure Python implementation of the Needleman-Wunsch dynamic programming
algorithm for optimal global pairwise sequence alignment of biological
sequences (DNA, RNA, or Protein).

Features:
- Dynamic programming matrix initialization with linear gap penalties.
- Traceback logic supporting diagonal (match/mismatch), vertical (deletion/gap in seq2),
  and horizontal (insertion/gap in seq1) transitions.
- Customizable scoring scheme: match reward, mismatch penalty, and gap penalty.
- Formatted alignment and score matrix output for pedagogical and diagnostic use.
- Handles edge cases: empty sequences, asymmetrical sequence lengths, and complete mismatches.

Author: Azheem Syahmee
License: MIT
"""

from typing import List, Tuple, Dict, Any, Optional


def initialize_matrix(seq1: str, seq2: str, gap_penalty: int = -2) -> List[List[int]]:
    """
    Initialize the dynamic programming scoring matrix with gap penalty boundaries.

    Args:
        seq1: First sequence (rows, length n).
        seq2: Second sequence (columns, length m).
        gap_penalty: Penalty score for gap opening/extension (typically negative).

    Returns:
        (n + 1) x (m + 1) grid initialized with boundary gap penalties.
    """
    rows = len(seq1) + 1
    cols = len(seq2) + 1
    matrix = [[0 for _ in range(cols)] for _ in range(rows)]

    for i in range(1, rows):
        matrix[i][0] = i * gap_penalty
    for j in range(1, cols):
        matrix[0][j] = j * gap_penalty

    return matrix


def fill_matrix(
    matrix: List[List[int]],
    seq1: str,
    seq2: str,
    match_score: int = 1,
    mismatch_penalty: int = -1,
    gap_penalty: int = -2,
) -> None:
    """
    Fill the dynamic programming matrix using the Needleman-Wunsch recurrence relation.

    Args:
        matrix: Initialized scoring matrix modified in-place.
        seq1: First sequence string.
        seq2: Second sequence string.
        match_score: Reward for matching characters.
        mismatch_penalty: Penalty for mismatched characters.
        gap_penalty: Penalty for introducing a gap.
    """
    rows = len(seq1) + 1
    cols = len(seq2) + 1

    for i in range(1, rows):
        for j in range(1, cols):
            is_match = seq1[i - 1] == seq2[j - 1]
            score_diag = matrix[i - 1][j - 1] + (match_score if is_match else mismatch_penalty)
            score_up = matrix[i - 1][j] + gap_penalty
            score_left = matrix[i][j - 1] + gap_penalty
            matrix[i][j] = max(score_diag, score_up, score_left)


def traceback(
    matrix: List[List[int]],
    seq1: str,
    seq2: str,
    match_score: int = 1,
    mismatch_penalty: int = -1,
    gap_penalty: int = -2,
) -> Tuple[str, str]:
    """
    Trace back through the scoring matrix from (len(seq1), len(seq2)) to (0, 0)
    to reconstruct the optimal global alignment.

    Args:
        matrix: Filled dynamic programming scoring matrix.
        seq1: First sequence.
        seq2: Second sequence.
        match_score: Match reward.
        mismatch_penalty: Mismatch penalty.
        gap_penalty: Gap penalty.

    Returns:
        Tuple of (aligned_seq1, aligned_seq2).
    """
    aligned_seq1 = ""
    aligned_seq2 = ""
    i = len(seq1)
    j = len(seq2)

    while i > 0 and j > 0:
        current_score = matrix[i][j]
        diagonal_score = matrix[i - 1][j - 1]
        up_score = matrix[i - 1][j]
        left_score = matrix[i][j - 1]

        expected_diag = diagonal_score + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_penalty)

        if current_score == expected_diag:
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            i -= 1
            j -= 1
        elif current_score == up_score + gap_penalty:
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            i -= 1
        else:
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            j -= 1

    while i > 0:
        aligned_seq1 = seq1[i - 1] + aligned_seq1
        aligned_seq2 = "-" + aligned_seq2
        i -= 1

    while j > 0:
        aligned_seq1 = "-" + aligned_seq1
        aligned_seq2 = seq2[j - 1] + aligned_seq2
        j -= 1

    return aligned_seq1, aligned_seq2


def needleman_wunsch(
    seq1: str,
    seq2: str,
    match_score: int = 1,
    mismatch_penalty: int = -1,
    gap_penalty: int = -2,
) -> Dict[str, Any]:
    """
    Execute complete Needleman-Wunsch global alignment on two biological sequences.

    Args:
        seq1: First sequence string (e.g., DNA 'GATTACA' or protein 'FKHMEDPLE').
        seq2: Second sequence string.
        match_score: Reward score for identity match (default: 1).
        mismatch_penalty: Penalty score for mismatch (default: -1).
        gap_penalty: Penalty score for gap insertion/deletion (default: -2).

    Returns:
        Dictionary containing:
        - 'aligned_seq1': Aligned first sequence with gap characters.
        - 'aligned_seq2': Aligned second sequence with gap characters.
        - 'score': Total global alignment score.
        - 'matrix': Full dynamic programming scoring matrix.
        - 'identity': Fraction of identical aligned positions (0.0 to 1.0).
        - 'similarity_str': Match visualization string ('|' match, '.' mismatch, ' ' gap).
    """
    seq1_clean = seq1.strip().upper()
    seq2_clean = seq2.strip().upper()

    if not seq1_clean and not seq2_clean:
        return {
            "aligned_seq1": "",
            "aligned_seq2": "",
            "score": 0,
            "matrix": [[0]],
            "identity": 0.0,
            "similarity_str": "",
        }

    if not seq1_clean:
        return {
            "aligned_seq1": "-" * len(seq2_clean),
            "aligned_seq2": seq2_clean,
            "score": len(seq2_clean) * gap_penalty,
            "matrix": [[j * gap_penalty for j in range(len(seq2_clean) + 1)]],
            "identity": 0.0,
            "similarity_str": " " * len(seq2_clean),
        }

    if not seq2_clean:
        return {
            "aligned_seq1": seq1_clean,
            "aligned_seq2": "-" * len(seq1_clean),
            "score": len(seq1_clean) * gap_penalty,
            "matrix": [[i * gap_penalty] for i in range(len(seq1_clean) + 1)],
            "identity": 0.0,
            "similarity_str": " " * len(seq1_clean),
        }

    matrix = initialize_matrix(seq1_clean, seq2_clean, gap_penalty=gap_penalty)
    fill_matrix(
        matrix,
        seq1_clean,
        seq2_clean,
        match_score=match_score,
        mismatch_penalty=mismatch_penalty,
        gap_penalty=gap_penalty,
    )
    aligned_seq1, aligned_seq2 = traceback(
        matrix,
        seq1_clean,
        seq2_clean,
        match_score=match_score,
        mismatch_penalty=mismatch_penalty,
        gap_penalty=gap_penalty,
    )
    score = matrix[len(seq1_clean)][len(seq2_clean)]

    # Calculate identity and match representation
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
        "score": score,
        "matrix": matrix,
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
    # Test sample execution
    sample_seq1 = "FKHMEDPLE"
    sample_seq2 = "FMDTPLNE"

    res = needleman_wunsch(sample_seq1, sample_seq2, match_score=1, mismatch_penalty=-2, gap_penalty=-2)
    print("=== Needleman-Wunsch Global Alignment ===")
    print("Sequence 1:", sample_seq1)
    print("Sequence 2:", sample_seq2)
    print("Aligned 1 :", res["aligned_seq1"])
    print("Match     :", res["similarity_str"])
    print("Aligned 2 :", res["aligned_seq2"])
    print(f"Score     : {res['score']}")
    print(f"Identity  : {res['identity'] * 100:.1f}%\n")
    print("Dynamic Programming Matrix:")
    print(format_matrix(res["matrix"], sample_seq1, sample_seq2))
