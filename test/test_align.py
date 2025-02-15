# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    Unit test for the Needleman-Wunsch alignment function.
    
    Uses test_seq1.fa and test_seq2.fa with BLOSUM62 matrix,
    a gap open penalty of -10, and a gap extension penalty of -1.
    
    Asserts:
        - The first matrix element is initialized correctly.
        - The first gap matrix elements are set correctly.
        - The last value of the alignment matrix is finite.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    nw = NeedlemanWunsch("substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)
    
    nw.align(seq1, seq2)

    # Check initial conditions
    assert nw._align_matrix[0, 0] == 0
    assert nw._gapA_matrix[0, 0] == nw.gap_open
    assert nw._gapB_matrix[0, 0] == nw.gap_open

    # Check that the matrices are populated correctly
    assert np.isfinite(nw._align_matrix[-1, -1])  # Ensure last value is computed

    # Check gap initialization along first row/column
    assert nw._gapA_matrix[1, 0] == nw.gap_open + nw.gap_extend
    assert nw._gapB_matrix[0, 1] == nw.gap_open + nw.gap_extend

def test_nw_backtrace():
    """
    Unit test for the Needleman-Wunsch backtrace function.

    Uses test_seq3.fa and test_seq4.fa with BLOSUM62 matrix,
    a gap open penalty of -10, and a gap extension penalty of -1.

    Asserts:
        - The backtrace produces the expected alignment.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    nw = NeedlemanWunsch("substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)
    
    alignment_score, aligned_seq3, aligned_seq4 = nw.align(seq3, seq4)

    # Expected results (Modify these as needed)
    expected_seq3 = "M--QKL"
    expected_seq4 = "MGGQKL"
    
    assert aligned_seq3 == expected_seq3, f"Expected: {expected_seq3}, but got: {aligned_seq3}"
    assert aligned_seq4 == expected_seq4, f"Expected: {expected_seq4}, but got: {aligned_seq4}"
    assert np.isfinite(alignment_score)  # Ensure score is a finite number





