from global_affine import *
import numpy as np

"""
Compute the optimal score of an optimal alignment for each pair of the
5 sequences above using program global_affine with the score matrix M 
and gap cost g(k)=5+5*k using your program global_affine. The result
is a 5x5 table where entry (i,j) the optimal score of an alignment of 
seq_i and seq_j.
"""
    
def compute_optimal_score(sequences: list, alphabet: list, cost_matrix: list, gap_cost: int):
    n = len(sequences)
    result = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            score_mat, _ = global_affine_alignment([sequences[i], sequences[j]], alphabet, cost_matrix, a, b)
            score = score_mat[-1][-1]
            result[i][j] = score
            result[j][i] = score
    return result

if __name__ == "__main__":
    eval_directory = "eval/"
    cost_matrix_filename = 'cost_matrix.txt'
    gapcost, alphabet, cost_matrix = read_cost_matrix(eval_directory + cost_matrix_filename)
    a, b = gapcost
    sequences = read_fasta(eval_directory + "q1.fasta", 6)
    check_sequences_validity(sequences, alphabet)
    result = compute_optimal_score(sequences, alphabet, cost_matrix, gapcost)
    for row in result:
        print("".join(f"{val:>5}" for val in row))