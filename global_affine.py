import numpy as np    
from alignment_utils import *
from math import inf

def backtrack_from_score_matrix(S: list, D: list, I: list, sequences: list, alphabet: list, cost_matrix: list, a: int, b: int):
    align1 = ""
    align2 = ""
    A = sequences[0]
    B = sequences[1]
    n = len(A)
    m = len(B)
    i = n
    j = m
    # Check starting position
    if S[n][m] < D[n][m] and S[n][m] < I[n][m]:
        state = "S"
    elif D[n][m] < I[n][m]:
        state = "D"
    else:
        state = "I"
        
    while i > 0 or j > 0:
        if i == 0:
            align1 += "-"
            align2 += B[j-1]
            j -= 1
            continue
        if j == 0:
            align1 += A[i-1]
            align2 += "-"
            i -= 1
            continue
               
        if state == "S":
            cost = cost_matrix[alphabet.index(A[i-1])][alphabet.index(B[j-1])]
            if i > 0 and j > 0 and S[i][j] == S[i-1][j-1] + cost:
                align1 += A[i-1]
                align2 += B[j-1]
                i -= 1    
                j -= 1   
                state = "S" 
            elif S[i][j] == D[i][j]:
                state = "D"
            elif S[i][j] == I[i][j]:
                state = "I"
                     
        elif state == "D":
            if i > 0 and D[i][j] == D[i-1][j] + b:
                state = "D"
            else:
                state = "S"         
            align1 += A[i-1]
            align2 += "-"
            i -= 1
            
        elif state == "I":
            if j > 0 and I[i][j] == I[i][j-1] + b:
                state = "I"
            else:
                state = "S"
            align1 += "-"
            align2 += B[j-1]
            j -= 1

    return align1[::-1], align2[::-1] #The alignments were filled from back to front, so we must invert the string of each


def global_affine_alignment(sequences: list, alphabet: list, cost_matrix: list, a: int, b: int, backtracking: bool = False):
    """
    Function to compute global alignment score with affine gap costs.
    It receives:
        - sequences: a list of two sequences that we want to align
        - alphabet: a list of the current alphabet, aligned with the index of the cost matrix
        - cost_matrix: a matrix with the scores of each substitution
        - a: gap opening penalty
        - b: gap extension penalty
        - backtracking: a boolean variable to return an optimal alignment using backtracking, if desired 
    """
    A = sequences[0]
    B = sequences[1]
    n = len(A)
    m = len(B)
    S = [[inf] * (m+1) for _ in range(n+1)] # score matrix
    D = [[inf] * (m+1) for _ in range(n+1)] # deletion matrix
    I = [[inf] * (m+1) for _ in range(n+1)] # insertion matrix
    
    # Initialization
    S[0][0] = 0
    # Gaps in Seq 1 (Vertical)
    for i in range(1, n + 1):
        D[i][0] = a + b * i
        S[i][0] = D[i][0]
    # Gaps in Seq 2 (Horizontal)
    for j in range(1, m + 1):
        I[0][j] = a + b * j
        S[0][j] = I[0][j]
    
    for i in range(1, n+1):
        for j in range(1, m+1):
            if i == 0 and j == 0:
                continue
            
            # Update D[i][j]
            # We look at the row above (i-1)
            if i > 0:
                D[i][j] = min(S[i-1][j] + a + b, D[i-1][j] + b)
            
            # Update I[i][j]
            # We look at the column to the left (j-1)
            if j > 0:
                I[i][j] = min(S[i][j-1] + a + b, I[i][j-1] + b)
                
            # Update S[i][j]
            if i > 0 and j > 0:
                cost = cost_matrix[alphabet.index(A[i-1])][alphabet.index(B[j-1])] 
                match_score = S[i-1][j-1] + cost
                S[i][j] = min(match_score, D[i][j], I[i][j])
    
    if backtracking:
        alignment = backtrack_from_score_matrix(S, D, I, sequences, alphabet, cost_matrix, a,  b)
        return S, alignment
    else:
        return S, ["", ""]
    
if __name__ == "__main__":
    tests_directory = 'tests/'
    filename = "test4.fasta"
    cost_matrix_filename = 'cost_matrix.txt'
    gapcost, alphabet, cost_matrix = read_cost_matrix(tests_directory + cost_matrix_filename)
    a, b = gapcost
    sequences = read_fasta(tests_directory + filename)
    backtracking = True
    check_sequences_validity(sequences, alphabet)
    score_matrix, alignment = global_affine_alignment(sequences, alphabet, cost_matrix, a, b, backtracking)
    
    print(f"The cost of a global alignment is {score_matrix[-1][-1]}")
    write_alignment_in_fasta(alignment, tests_directory + "output_" + filename)