from alignment_utils import *
from sp_approximation import *
from sp_exact_3 import *
import time

if __name__ == "__main__":
    eval_directory = "testseqs/"
    cost_matrix_filename = 'cost_matrix.txt'
    gapcost, alphabet, cost_matrix = read_cost_matrix(eval_directory +cost_matrix_filename)
    
    for i in range(10,210,10):
        filename = f"testseqs_{i}_3.fasta"
        sequences = read_fasta(eval_directory+filename, num_sequences=3)
        check_sequences_validity(sequences.values(), alphabet)

        start_time = time.time()
        multiple_alignment, sequence_names = two_sp_approximation(sequences, alphabet, cost_matrix, gapcost)
        end_time = time.time()

        write_alignment_in_fasta(multiple_alignment, eval_directory + "approximation_output_"+filename , sequence_names, comment=f"time taken {end_time-start_time}")

        start_time = time.time()
        score_matrix, exact_alignment = sp_exact_3(sequences, alphabet, cost_matrix, gapcost, backtracking=True)
        score = score_matrix[-1][-1][-1]
        end_time = time.time()

        write_alignment_in_fasta(exact_alignment, eval_directory + "exact_output_"+filename , sequence_names, comment=f"time taken {end_time-start_time}")
