import numpy as np    
from alignment_utils import *

if __name__ == "__main__":
    eval_directory = "full_seq_alignment/"
    filename = "brca1-full.fasta"
    cost_matrix_filename = 'cost_matrix_capital.txt'
    sequences = read_fasta(eval_directory+filename, num_sequences=8)
    gapcost, alphabet, cost_matrix = read_cost_matrix(eval_directory + cost_matrix_filename)
    print(sequences)
    process_sequences(sequences, alphabet, eval_directory + filename)
