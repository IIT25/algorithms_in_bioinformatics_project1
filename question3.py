from global_general import *
from global_affine import *

"""
Compute and show an optimal alignment of sequence 1 and 2 using gap 
cost g(k)=5+5*k using your programs general_affine and general_general.
"""
eval_directory = "eval/"
cost_matrix_filename = 'cost_matrix.txt'
gapcost_params, alphabet, cost_matrix = read_cost_matrix(eval_directory + cost_matrix_filename)
a, b = gapcost_params
gapcost = affine_gapcost


sequences = read_fasta(eval_directory + "q1.fasta", 2)

check_sequences_validity(sequences, alphabet)

print("Global general function")
score_mat, alignment = global_general_alignment([sequences[0], sequences[1]], alphabet, cost_matrix, gapcost, gapcost_params, True)
print(">seq1")
print(alignment[0])
print(">seq2")
print(alignment[1])
print()
print("Global affine function")
score_mat, alignment = global_affine_alignment([sequences[0], sequences[1]], alphabet, cost_matrix, a, b, True)
print(">seq1")
print(alignment[0])
print(">seq2")
print(alignment[1])