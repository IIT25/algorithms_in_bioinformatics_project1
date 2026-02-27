import numpy as np    
from alignment_utils import *

def affine_gapcost(k, a, b):
    """
    Function implementing the affine gapcost. It recieves
        - k: The length of the gap
        - a: affine cost parameter
        - b: affine cost parameter
    And returns the affine cost a+b*k
    """
    if k == 0:
        return 0
    return a+b*k


def global_affine_alignment(sequences: list, alphabet: list, cost_matrix: list, a: int, b: int, backtracking: bool = False):
    """
    Function to get the global pairwise alignment using linear gap cost
    It receives:
        - sequences: a list of two sequences that we want to align
        - alphabet: a list of the current alphabet, aligned with the index of the cost matrix
        - cost_matrix: a matrix with the scores of each substitution
        - a: parameter of the affine gap cost, where the affine function is a+b*k
        - b: parameter of the affine gap cost, where the affine function is a+b*k
        - backtracking: a boolean variable to return an optimal alignment using backtracking, if desired 
    """
    n = len(sequences[0])
    m = len(sequences[1])
    score_matrix = [[0 for j in range(n+1)] for i in range(m+1)]
    deletion_matrix = [[0 for j in range(n+1)] for i in range(m+1)]
    insertion_matrix = [[0 for j in range(n+1)] for i in range(m+1)]
    for i in range(m+1):
        for j in range(n+1):
            #I'm not sure if there is a better way to implement it
            #Maybe more efficient
            score_options = []
            deletion_options = []
            insertion_options = []
            if i == 0 and j == 0:
                deletion_options.append(0)
                insertion_options.append(0)
                score_options.append(0)
            #Fill the deletion matrix
            if i > 0 and j >= 0:
                deletion_options.append(score_matrix[i-1][j]+a+b)
            if i > 1 and j >= 0:
                deletion_options.append(deletion_matrix[i-1][j]+a)
            if i > 0:
                deletion_matrix[i][j] = min(deletion_options)

            #Fill the insertion matrix
            if i >= 0 and j > 0:
                insertion_options.append(score_matrix[i][j-1]+a+b)
            if i >= 0 and j > 1:
                insertion_options.append(insertion_matrix[i][j-1] + a)
            if j > 0:
                insertion_matrix[i][j] = min(insertion_options)

            #Fill the score matrix
            if i > 0 and j > 0:
                score_options.append(score_matrix[i-1][j-1] + cost_matrix[alphabet.index(sequences[0][j-1])][alphabet.index(sequences[1][i-1])])
            if i > 0 and j >= 0:
                score_options.append(deletion_matrix[i][j])
            if i >= 0 and j > 0:
                score_options.append(insertion_matrix[i][j]) 
            score_matrix[i][j] = min(score_options)
    print(np.array(score_matrix))

tests_directory = 'tests/'
filename = "test1.fasta"
cost_matrix_filename = 'cost_matrix.txt'
gapcost_params, alphabet, cost_matrix = read_cost_matrix(tests_directory + cost_matrix_filename)
a, b = gapcost_params
sequences = read_fasta(tests_directory + filename)
backtracking = True
check_sequences_validity(sequences, alphabet)
score_matrix, alignment = global_affine_alignment(sequences, alphabet, cost_matrix, a, b, backtracking)

print(f"The cost of a global alignment is {score_matrix[-1][-1]}")
write_alignment_in_fasta(alignment, tests_directory + "output_" + filename)

print(np.array(score_matrix))
print(alignment)