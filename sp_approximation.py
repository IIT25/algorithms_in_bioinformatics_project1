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
        return score_matrix, [list(seq) for seq in alignment]
    else:
        return score_matrix, ["", ""]
    
def compute_optimal_score(sequences: list, alphabet: list, cost_matrix: list, gap_cost: int):
    n = len(sequences)
    result = np.zeros((n,n))
    for i in range(n):
        for j in range(i, n):
            score_mat, _ = global_linear([sequences[i], sequences[j]], alphabet, cost_matrix, gap_cost)
            score = score_mat[-1][-1]
            result[i][j] = score
            result[j][i] = score
    return result

def merge_alignments(multiple_alignment, pairwise_alignment):
    """
    Function to merge a multiple alignment to a pairwise alignment
    It recieves
      - multiple alignment: A matrix of dimensions (k,n)
      - pairwise_alignment: A matrix of dimensions (2,m)

    This method will assume that the first sequence of the multiple
    alignment is the same as the first row of the pairwise alignment
    And will traverse them following the logic we saw in class
    """
    if multiple_alignment == []:
        return pairwise_alignment
    n = len(multiple_alignment[0])
    m = len(pairwise_alignment[0])
    k = len(multiple_alignment)
    i = 0 
    j = 0
    final_multiple_alignment = [[] for _ in range(k+1)]
    while i < n or j < m:
        if multiple_alignment[0][i] == pairwise_alignment[0][j]:
            for l in range(k):
                final_multiple_alignment[l].append(multiple_alignment[l][i])
            final_multiple_alignment[-1].append(pairwise_alignment[-1][j])
            i+=1
            j+=1
        elif multiple_alignment[0][i] == '-':
            for l in range(k):
                final_multiple_alignment[l].append(multiple_alignment[l][i])
            final_multiple_alignment[-1].append('-')
            i+=1
        else:
            for l in range(k):
                final_multiple_alignment[l].append('-')
            final_multiple_alignment[-1].append(pairwise_alignment[-1][j])
            j+=1
    return final_multiple_alignment

def two_sp_approximation(sequences, alphabet, cost_matrix, gap_cost):
    """
    Method that creates a 2-sp approximation for a multiple sequence alignment
    It creates it by using a star-like guide tree, finding a central string 
    and then running a pairwise alignment of the central string to every other string

    It recieves:
        - sequences: A list containing all the sequences to be aligned
        - alphabet: a list of the current alphabet, aligned with the index of the cost matrix
        - cost_matrix: a matrix with the scores of each substitution
        - gap_cost: an integer representing the gap cost on the alignment
    """
    result = compute_optimal_score(sequences, alphabet, cost_matrix, gap_cost)
    central_sequence = np.argmax(result.sum(axis = 0))
    multiple_alignment = []
    for i in range(len(sequences)):
        if i != central_sequence:
            _, pairwise_alignment = global_linear([sequences[central_sequence], sequences[i]], alphabet, cost_matrix, gap_cost, True)
            multiple_alignment = merge_alignments(multiple_alignment, pairwise_alignment)
    for i in range(len(multiple_alignment)):
        multiple_alignment[i] = "".join(multiple_alignment[i])
    
    return multiple_alignment

if __name__ == "__main__":
    eval_directory = "eval/"
    cost_matrix_filename = 'cost_matrix.txt'
    gapcost, alphabet, cost_matrix = read_cost_matrix(eval_directory + cost_matrix_filename)
    sequences = []
    for i in range(1, 6):
        seq_from_file = read_fasta(eval_directory + "seq" + str(i) + ".fasta", 1)
        sequences.append(seq_from_file[0])
        
    check_sequences_validity(sequences, alphabet)
    multiple_alignment = two_sp_approximation(sequences, alphabet, cost_matrix, gapcost)
    for row in multiple_alignment:
        print(row)

    sequences = read_fasta("tests/brca1-testseqs.fasta", num_sequences=8)
    eval_directory = "tests/"
    cost_matrix_filename = 'cost_matrix_capital.txt'
    gapcost, alphabet, cost_matrix = read_cost_matrix(eval_directory + cost_matrix_filename)
    check_sequences_validity(sequences, alphabet)
    multiple_alignment = two_sp_approximation(sequences, alphabet, cost_matrix, gapcost)
    for row in multiple_alignment:
        print(row)