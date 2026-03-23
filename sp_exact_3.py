import numpy as np    
from alignment_utils import *


def backtrack(seq1: str, seq2: str, seq3: str,  score_matrix: list, cost_matrix: list, gap_cost: int, alphabet: list):
    i, j, k = len(seq1), len(seq2), len(seq3)
    align1, align2, align3 = "", "", ""
    
    while i > 0 or j > 0 or k > 0:
        current_score = score_matrix[i, j, k]
        charA = seq1[i-1] if i > 0 else "-"
        charB = seq2[j-1] if j > 0 else "-"
        charC = seq3[k-1] if k > 0 else "-"
        
        # match/mismatch in all three
        if i > 0 and j > 0 and k > 0 and \
            current_score == score_matrix[i-1, j-1, k-1] + get_sp_score(charA, charB, charC, cost_matrix, gap_cost, alphabet):
                align1 += charA
                align2 += charB
                align3 += charC
                i -= 1; j -= 1; k -= 1
                
        # one gap
        elif i > 0 and j > 0 and \
            current_score == score_matrix[i-1, j-1, k] + get_sp_score(charA, charB, "-", cost_matrix, gap_cost, alphabet):
                align1 += charA
                align2 += charB
                align3 += "-"
                i -= 1; j -= 1
        elif i > 0 and k > 0 and \
            current_score == score_matrix[i-1, j, k-1] + get_sp_score(charA, "-", charC, cost_matrix, gap_cost, alphabet):
                align1 += charA
                align2 += "-"
                align3 += charC
                i -= 1; k -= 1
        elif j > 0 and k > 0 and \
            current_score == score_matrix[i, j-1, k-1] + get_sp_score("-", charB, charC, cost_matrix, gap_cost, alphabet):
                align1 += "-"
                align2 += charB
                align3 += charC
                j -= 1; k -= 1
                
        # two gaps
        elif i > 0 and current_score == score_matrix[i-1, j, k] + get_sp_score(charA, "-", "-", cost_matrix, gap_cost, alphabet):
            align1 += charA
            align2 += "-"
            align3 += "-"
            i -= 1
        elif j > 0 and current_score == score_matrix[i, j-1, k] + get_sp_score("-", charB, "-", cost_matrix, gap_cost, alphabet):
            align1 += "-"
            align2 += charB
            align3 += "-"
            j -= 1
        elif k > 0 and current_score == score_matrix[i, j, k-1] + get_sp_score("-", "-", charC, cost_matrix, gap_cost, alphabet):
            align1 += "-"
            align2 += "-"
            align3 += charC
            k -= 1
    return align1[::-1], align2[::-1], align3[::-1] #The alignments were filled from back to front, so we must invert the string of each

def get_sp_score(char1: str, char2: str, char3: str, cost_matrix: list, gap_cost: int, alphabet):
    def pair_cost(a, b):
        if a == "-" and b == "-":
            return 0
        if a == "-" or b == "-":
            return gap_cost
        return cost_matrix[alphabet.index(a)][alphabet.index(b)]
    return pair_cost(char1, char2) + pair_cost(char1, char3) + pair_cost(char2, char3)

def sp_exact_3(sequences: list, alphabet: list, cost_matrix: list, gap_cost: int, backtracking: bool = False):
    """
    Function to get the global pairwise alignment using linear gap cost
    It receives:
        - sequences: a list of three sequences that we want to align
        - alphabet: a list of the current alphabet, aligned with the index of the cost matrix
        - cost_matrix: a matrix with the scores of each substitution
        - gap_cost: an integer representing the gap cost on the alignment
        - backtracking: a boolean variable to return an optimal alignment using backtracking, if desired 
    """
    A = sequences["seq1"]
    B = sequences["seq2"]
    C = sequences["seq3"]
    n = len(A) 
    m = len(B)
    o = len(C)
    score_matrix = np.zeros((n+1, m+1, o+1))

    for i in range(n+1):
        for j in range(m+1):
            for k in range(o+1):
                if i == 0 and j == 0 and k == 0:
                    score_matrix[i, j, k] = 0
                    continue
                                
                v1 = v2 = v3 = v4 = v5 = v6 = v7 = float("inf")
                # We use i-1 because the matrix index i corresponds to sequence char i-1
                charA = A[i-1] if i > 0 else "-"
                charB = B[j-1] if j > 0 else "-"
                charC = C[k-1] if k > 0 else "-"
                
                # match/mismatch in all three
                if i > 0 and j > 0 and k > 0:
                    v1 = score_matrix[i-1, j-1, k-1] + get_sp_score(charA, charB, charC, cost_matrix, gap_cost, alphabet)
                # one gap
                if i > 0 and j > 0:
                    v2 = score_matrix[i-1, j-1, k] + get_sp_score(charA, charB, "-", cost_matrix, gap_cost, alphabet)
                if i > 0 and k > 0:
                    v3 = score_matrix[i-1, j, k-1] + get_sp_score(charA, "-", charC, cost_matrix, gap_cost, alphabet)
                if j > 0 and k > 0:
                    v4 = score_matrix[i, j-1, k-1] + get_sp_score("-", charB, charC, cost_matrix, gap_cost, alphabet) 
                # two gaps
                if i > 0:
                    v5 = score_matrix[i-1, j, k] + get_sp_score(charA, "-", "-", cost_matrix, gap_cost, alphabet)
                if j > 0:
                    v6 = score_matrix[i, j-1, k] + get_sp_score("-", charB, "-", cost_matrix, gap_cost, alphabet)
                if k > 0:
                    v7 = score_matrix[i, j, k-1] + get_sp_score("-", "-", charC, cost_matrix, gap_cost, alphabet)
                
                score_matrix[i, j, k] = min(v1, v2, v3, v4, v5, v6, v7)
                
    if backtracking:
        alignment = backtrack(A, B, C, score_matrix, cost_matrix, gap_cost, alphabet)
        return score_matrix, alignment
    else:
        return score_matrix, ["", "", ""]

if __name__ == "__main__":
    eval_directory = "tests/"
    filename = "testdata_long.fasta"
    cost_matrix_filename = 'cost_matrix_capital.txt'
    sequences = read_fasta(eval_directory+filename, num_sequences=3)
    gapcost, alphabet, cost_matrix = read_cost_matrix(eval_directory + cost_matrix_filename)
    check_sequences_validity(sequences.values(), alphabet)
    score_matrix, alignment = sp_exact_3(sequences, alphabet, cost_matrix, gapcost, True)
    print(score_matrix[-1][-1][-1])
    for row in alignment:
        print("".join(row))
    write_alignment_in_fasta(alignment, eval_directory + "output_" + filename)