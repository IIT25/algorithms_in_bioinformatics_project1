import numpy as np    
from alignment_utils import *

def backtrack_from_score_matrix(i: int, j: int, sequences: list, score_matrix: list, cost_matrix: list, gapcost: int, alphabet: list, alignment1: str = '', alignment2: str = ''):
    if i > 0 and j > 0 and score_matrix[i][j] == score_matrix[i-1][j-1] + cost_matrix[alphabet.index(sequences[0][j-1])][alphabet.index(sequences[1][i-1])]:
        alignment1 += sequences[0][j-1]
        alignment2 += sequences[1][i-1]
        return backtrack_from_score_matrix(i-1, j-1, sequences, score_matrix, cost_matrix, gapcost, alphabet, alignment1, alignment2)
    elif i > 0 and j >= 0 and score_matrix[i][j] == score_matrix[i-1][j] + gapcost:
        alignment1 += '-'
        alignment2 += sequences[1][i-1]
        return backtrack_from_score_matrix(i-1, j, sequences, score_matrix, cost_matrix, gapcost, alphabet, alignment1, alignment2)
    elif i >= 0 and j > 0 and score_matrix[i][j] == score_matrix[i][j-1] + gapcost:
        alignment1 += sequences[0][j-1]
        alignment2 += '-'
        return backtrack_from_score_matrix(i, j-1, sequences, score_matrix, cost_matrix, gapcost, alphabet, alignment1, alignment2)
    return alignment1[::-1], alignment2[::-1] #The alignments were filled from back to front, so we must invert the string of each

def global_linear(sequences: list, alphabet: list, cost_matrix: list, gap_cost: int, backtracking: bool = False):
    """
    Function to get the global pairwise alignment using linear gap cost
    It receives:
        - sequences: a list of two sequences that we want to align
        - alphabet: a list of the current alphabet, aligned with the index of the cost matrix
        - cost_matrix: a matrix with the scores of each substitution
        - gap_cost: an integer representing the gap cost on the alignment
        - backtracking: a boolean variable to return an optimal alignment using backtracking, if desired 
    """
    n = len(sequences[0])
    m = len(sequences[1])
    score_matrix = [list(range(n+1))] + [[(i+1)* gap_cost] + [0] * n for i in range(m)]
    score_matrix[0] = [i*gap_cost for i in score_matrix[0]] 

    for i in range(1, m + 1):
        for j in range(1, n + 1):
                score_matrix[i][j] = min(score_matrix[i-1][j-1] + cost_matrix[alphabet.index(sequences[0][j-1])][alphabet.index(sequences[1][i-1])],
                                         score_matrix[i-1][j] + gap_cost,
                                         score_matrix[i][j-1] + gap_cost)
    if backtracking:
        alignment = backtrack_from_score_matrix(m,n,sequences, score_matrix, cost_matrix, gap_cost, alphabet)
        return score_matrix, alignment
    else:
        return score_matrix, ["", ""]


if __name__ == "__main__":
    tests_directory = 'tests/'
    filename = "test4.fasta"
    cost_matrix_filename = 'cost_matrix.txt'
    gapcost, alphabet, cost_matrix = read_cost_matrix(tests_directory + cost_matrix_filename)
    sequences = read_fasta(tests_directory + filename)
    backtracking = True
    check_sequences_validity(sequences, alphabet)
    score_matrix, alignment = global_linear(sequences, alphabet, cost_matrix, gapcost, backtracking)
    
    print(f"The cost of a global alignment is {score_matrix[-1][-1]}")
    write_alignment_in_fasta(alignment, tests_directory + "output_" + filename)