import numpy as np    
from alignment_utils import *
from sp_approximation import *
import sys
sys.setrecursionlimit(100000) 

if __name__ == "__main__":
    eval_directory = "presentation_sequences/"
    cost_matrix_filename = 'cost_matrix.txt'
    gapcost, alphabet, cost_matrix = read_cost_matrix(eval_directory +cost_matrix_filename)
    
    filename = f"brca1-full.fasta"
    sequences = read_fasta(eval_directory+filename, num_sequences=8)
    check_sequences_validity(sequences.values(), alphabet)

    multiple_alignment, sequence_names = two_sp_approximation(sequences, alphabet, cost_matrix, gapcost)
    write_alignment_in_fasta(multiple_alignment, eval_directory + "approximation_output_" + filename, sequence_names)