import numpy as np
from alignment_utils import * 
from global_linear import *

def global_count(sequences: list, alphabet: list, cost_matrix: list, gap_cost: int):
    """
    Function that counts all possible optimal pairwise alignments using linear gap cost
    It recieves:
        - sequences: a list of two sequences that we want to align
        - alphabet: a list of the current alphabet, aligned with the index of the cost matrix
        - cost_matrix: a matrix with the scores of each substitution
        - gap_cost: an integer representing the gap cost on the alignment
    """
    n = len(sequences[0])
    m = len(sequences[1])  
    score_matrix, _ = global_linear(sequences, alphabet, cost_matrix, gap_cost, backtracking=False)
    count_matrix = np.zeros((m + 1, n + 1), dtype=int)
    count_matrix[0, :] = 1
    count_matrix[:, 0] = 1
    
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            current_score = score_matrix[i][j]
            
            if current_score == score_matrix[i-1][j-1] + cost_matrix[alphabet.index(sequences[0][j-1])][alphabet.index(sequences[1][i-1])]:
                count_matrix[i][j] += count_matrix[i-1][j-1]
                
            if current_score == score_matrix[i-1][j] + gap_cost:
                count_matrix[i][j] += count_matrix[i-1][j]
            
            if current_score == score_matrix[i][j-1] + gap_cost:
                count_matrix[i][j] += count_matrix[i][j-1]
                
    return count_matrix[m][n]    
            
            
if __name__ == "__main__":
    tests_directory = 'tests/'
    filename = "test1.fasta"
    cost_matrix_filename = 'cost_matrix.txt'
    gapcost, alphabet, cost_matrix = read_cost_matrix(tests_directory + cost_matrix_filename)
    sequences = read_fasta(tests_directory + filename)
    check_sequences_validity(sequences, alphabet)
    
    total_count = global_count(sequences, alphabet, cost_matrix, gapcost)
    print(f"There are {total_count} possible optimal alignments.")

