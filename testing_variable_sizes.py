from alignment_utils import *
from global_linear import *
import numpy as np
from time import time
import pickle

def generate_random_sequences(n: int, alphabet: list):
    """
    function to generate a pair of random sequences to be aligned, 
    for simplicity and interpretability of the plot, we will assume that n=m

    Recieves:
        - n: The length of the sequences to be aligned
        - alphabet: A list with the alphabet considered to crete the sequences
    """
    seq1 = "".join([np.random.choice(alphabet) for _ in range(n)])
    seq2 = "".join([np.random.choice(alphabet) for _ in range(n)])
    return [seq1,seq2]

sequence_length = range(1,2000,10)
print(sequence_length)
n = []
t_n = []

tests_directory = 'tests/'
cost_matrix_filename = 'cost_matrix.txt'
gapcost, alphabet, cost_matrix = read_cost_matrix(tests_directory + cost_matrix_filename)
for length in sequence_length:
    for _ in range(5):
        t0 = time()
        sequences = generate_random_sequences(length, alphabet)
        backtracking = False
        score_matrix, alignment = global_linear(sequences, alphabet, cost_matrix, gapcost, backtracking)
        t1 = time()
        print(f"time taken for input of size {length}: {t1-t0}")
        n.append(length)
        t_n.append(t1-t0)

final_result = {"size": n, "time": t_n}
filehandler = open("tests/result.pkl","wb")
pickle.dump(final_result,filehandler)
filehandler.close()