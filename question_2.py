from alignment_utils import *
from sp_approximation import *
import time

if __name__ == "__main__":
    eval_directory = "eval/"
    cost_matrix_filename = 'cost_matrix.txt'
    gapcost, alphabet, cost_matrix = read_cost_matrix(eval_directory +cost_matrix_filename)
    
    filename = "brca1-testseqs.fasta"
    sequences = read_fasta(eval_directory+filename, num_sequences=5)
    check_sequences_validity(sequences.values(), alphabet)

    start_time = time.time()
    multiple_alignment, sequence_names = two_sp_approximation(sequences, alphabet, cost_matrix, gapcost)
    end_time = time.time()

    write_alignment_in_fasta(multiple_alignment, eval_directory + "approximation_output_question_2.fasta" , sequence_names)
    #The cost of this alignment after looking at it with the msa_sp_score_3k.py is 4384