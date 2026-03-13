from global_general import *
from global_affine import *
import numpy as np
import time

"""
Compute the optimal score of an optimal alignment for each pair of the
5 sequences above using program global_general with the score matrix M 
and gap cost g(k)=5+5*k using your program global_affine. The result
is a 5x5 table where entry (i,j) the optimal score of an alignment of 
seq_i and seq_j.
"""
def compute_optimal_score(sequences: list, function: str):
    n = len(sequences)
    result = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            if function == "general":
                score_mat, _ = global_general_alignment([sequences[i], sequences[j]], alphabet, cost_matrix, gapcost, gapcost_params)
            elif function == "affine":
                score_mat, _ = global_affine_alignment([sequences[i], sequences[j]], alphabet, cost_matrix, a, b)
            score = score_mat[-1][-1]
            result[i][j] = score
            result[j][i] = score
    return result

if __name__ == "__main__":
    eval_directory = "eval/"
    cost_matrix_filename = 'cost_matrix.txt'
    gapcost_params, alphabet, cost_matrix = read_cost_matrix(eval_directory + cost_matrix_filename)
    a, b = gapcost_params
    gapcost = affine_gapcost
    sequences = read_fasta(eval_directory + "q1.fasta", 6)
        
    check_sequences_validity(sequences, alphabet)

    # --- Global General ---
    print("Global general function -> O(n^3)")
    start_general = time.time()
    result = compute_optimal_score(sequences, "general")
    end_general = time.time()
    for row in result:
        print("".join(f"{val:>5}" for val in row))
    print(f"It took {end_general - start_general:.4f} seconds.\n")  
    
    # --- Global Affine ---
    print("Global affine function -> O(n^2)")
    start_general = time.time()
    result = compute_optimal_score(sequences, "affine")
    end_general = time.time()
    for row in result:
        print("".join(f"{val:>5}" for val in row))
    print(f"It took {end_general - start_general:.4f} seconds.")  