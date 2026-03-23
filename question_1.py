from alignment_utils import *
from sp_exact_3 import *

if __name__ == "__main__":
    eval_directory = "eval/"
    cost_matrix_filename = 'cost_matrix.txt'
    gapcost, alphabet, cost_matrix = read_cost_matrix(eval_directory + cost_matrix_filename)
    
    filename = "brca1-testseqs.fasta"
    sequences = read_fasta(eval_directory+filename, num_sequences=3)
    sequence_names = list(sequences.keys())
    check_sequences_validity(sequences.values(), alphabet)
    score_matrix, alignment = sp_exact_3(sequences, alphabet, cost_matrix, gapcost, True)
    print(score_matrix[-1][-1][-1]) # score is 790
    for row in alignment:
        print("".join(row))
    write_alignment_in_fasta(alignment, eval_directory + "exact_output_question_1.fasta", sequence_names)