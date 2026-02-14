from global_linear import *
import numpy as np

"""
Compute the optimal score of an optimal alignment for each pair of the
5 sequences above using global_linear with the score matrix M and gap
cost g(k)=5*k. The result is a 5x5 table where entry (i,j) the optimal
score of an alignment of seq_i and seq_j.
"""
    
def compute_optimal_score(sequences: list, alphabet: list, cost_matrix: list, gap_cost: int):
    n = len(sequences)
    result = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            score_mat, _ = global_linear([sequences[i], sequences[j]], alphabet, cost_matrix, gap_cost)
            score = score_mat[-1][-1]
            result[i][j] = score
            result[j][i] = score
    return result

if __name__ == "__main__":
    eval_directory = "eval/"
    cost_matrix_filename = 'cost_matrix.txt'
    gapcost, alphabet, cost_matrix = read_cost_matrix(eval_directory + cost_matrix_filename)
    sequences = []
    for i in range(1, 6):
        seq_from_file = read_fasta(eval_directory + "seq" + str(i) + ".fasta", 1)
        sequences.append(seq_from_file[0])
        
    check_sequences_validity(sequences, alphabet)
    result = compute_optimal_score(sequences, alphabet, cost_matrix, gapcost)
    for row in result:
        print("".join(f"{val:>5}" for val in row))