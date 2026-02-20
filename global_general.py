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

def backtrack_from_score_matrix(sequences: list, alphabet: list, score_matrix: list, cost_matrix: list, gapcost: int, gapcost_params):
    """
    This function implements an iterative backtracking for general gapcost functions in cubic time
    To do so, it looks first for replacements and then evaluates gaps going up or to the left in the score matrix
    It recieves:
    
        - sequences: The two sequences to be aligned
        - alphabet: a list of the current alphabet, aligned with the index of the cost matrix
        - score_matrix: The score matrix resulting with the cost of the alignments
        - cost_matrix: The matrix with the scores of each substitution
        - gapcost: The gapcost function considered
        - gapcost_params: a list with the parameters [a,b] for the gapcost function, for the affine gapcost it has the form a+b*k
    
    Returns the two sequences aligned
    """
    n = len(sequences[0])
    m = len(sequences[1])
    i = m
    j = n
    alignment1 = '' 
    alignment2 = ''
    #print(i,j)
    while i > 0 or j > 0:
        if i > 0 and j > 0 and score_matrix[i][j] == score_matrix[i-1][j-1] + cost_matrix[alphabet.index(sequences[0][j-1])][alphabet.index(sequences[1][i-1])]:
            alignment1 += sequences[0][j-1]
            alignment2 += sequences[1][i-1]
            i -= 1
            j -= 1
        else:
            k = 1
            while 1:
                if i >= k and score_matrix[i][j] == score_matrix[i-k][j] + gapcost(k, *gapcost_params):
                    alignment1 += "-"*k
                    alignment2 += sequences[1][i-k:i]
                    i -= k
                    break
                elif j >= k and score_matrix[i][j] == score_matrix[i][j-k] + gapcost(k, *gapcost_params):
                    alignment1 += sequences[0][j-k: j]
                    alignment2 += '-'*k
                    j -= k
                    break
                else:
                    k += 1
    return alignment1[::-1], alignment2[::-1] #The alignments were filled from back to front, so we must invert the string of each


def global_general_alignment(sequences: list, alphabet: list, cost_matrix: list, gap_cost, gapcost_params: list, backtracking: bool = False):
    """
    Function to get the global pairwise alignment using linear gap cost
    It recieves:
        - sequences: a list of two sequences that we want to align
        - alphabet: a list of the current alphabet, aligned with the index of the cost matrix
        - cost_matrix: a matrix with the scores of each substitution
        - gap_cost: a function that will be called to calculate the gapcost 
        - gapcost_params: a list with the parameters [a,b] for the gapcost function, for the affine gapcost it has the form a+b*k
        - backtracking: a boolean variable to return an optimal alignment using backtracking, if desired 
    """
    n = len(sequences[0])
    m = len(sequences[1])
    score_matrix = [[gap_cost(i, *gapcost_params)] + [0] * n for i in range(m+1)]
    score_matrix[0] = [gap_cost(i, *gapcost_params) for i in range(n+1)] 
    #Columns [1:n] seq 0 - j
    #Rows [1:m] seq 1 - i
    for i in range(1, m + 1):
        for j in range(1, n + 1):
                score_matrix[i][j] = min(score_matrix[i-1][j-1] + cost_matrix[alphabet.index(sequences[0][j-1])][alphabet.index(sequences[1][i-1])],
                                         min([score_matrix[i-k][j] + gap_cost(k, *gapcost_params) for k in range(1,i+1)]),
                                         min([score_matrix[i][j-k] + gap_cost(k, *gapcost_params) for k in range(1,j+1)]))
    #print(sequences[0])
    #print(sequences[1])
    #print(np.array(score_matrix)) #Convert to numpy just to view it pretty in terminal
    if backtracking:
        alignment = backtrack_from_score_matrix(sequences, alphabet, score_matrix, cost_matrix, gap_cost, gapcost_params)
        return score_matrix, alignment
    else:
        return score_matrix, ["", ""]


tests_directory = 'tests/'
filename = "test4.fasta"
cost_matrix_filename = 'cost_matrix.txt'
gapcost_params, alphabet, cost_matrix = read_cost_matrix(tests_directory + cost_matrix_filename)
gap_cost = affine_gapcost
sequences = read_fasta(tests_directory + filename)
backtracking = True
check_sequences_validity(sequences, alphabet)
score_matrix, alignment = global_general_alignment(sequences, alphabet, cost_matrix, gap_cost, gapcost_params, backtracking)

print(f"The cost of a global alignment is {score_matrix[-1][-1]}")
write_alignment_in_fasta(alignment, tests_directory + "output_" + filename)

print(alignment)