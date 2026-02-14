from global_linear import *
from global_count import *

"""
Compute and show an optimal alignment of sequence 1 and 2 using gap 
cost g(k)=5*k
"""
eval_directory = "eval/"
cost_matrix_filename = 'cost_matrix.txt'
gapcost, alphabet, cost_matrix = read_cost_matrix(eval_directory + cost_matrix_filename)

sequences = []
for i in range(1, 3):
    seq_from_file = read_fasta(eval_directory + "seq" + str(i) + ".fasta", 1)
    sequences.append(seq_from_file[0])
check_sequences_validity(sequences, alphabet)
    
score_mat, alignment = global_linear([sequences[0], sequences[1]], alphabet, cost_matrix, gapcost, True)
print(">seq1")
print(alignment[0])
print(">seq2")
print(alignment[1])
#do we need to print all possible alignments ?? 