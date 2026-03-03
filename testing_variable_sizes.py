from alignment_utils import *
from global_general import *
from global_affine import *
import numpy as np
from time import time
import pickle
def affine_gapcost(k: int, a: int, b:int):
    if k == 0:
        return 0
    return a + b * k

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

sequence_length = range(2,300,5)
print(sequence_length)
n = []
time_general = []


eval_directory = 'eval/'
cost_matrix_filename = 'cost_matrix.txt'
gapcost_params, alphabet, cost_matrix = read_cost_matrix(eval_directory + cost_matrix_filename)
a, b = gapcost_params
gapcost = affine_gapcost
backtracking = False

#for length in sequence_length:
#    for _ in range(5):
#        sequences = generate_random_sequences(length, alphabet)
#        t0 = time()
#        score_matrix, alignment = global_general_alignment(sequences, alphabet, cost_matrix, gapcost, gapcost_params, backtracking)
#        t1 = time()
#        print(f"time taken for general gaspcost function for input of size {length}: {t1-t0}")
#        time_general.append(t1-t0)
#        n.append(length)



#final_result = {"size": n, "time_general": time_general}
#filehandler = open("tests/result_general.pkl","wb")
#pickle.dump(final_result,filehandler)
#filehandler.close()


sequence_length = range(2,1000,5)
n = []
time_affine = []
for length in sequence_length:
    for _ in range(5):
        sequences = generate_random_sequences(length, alphabet)
        t0 = time()
        score_matrix, alignment = global_affine_alignment(sequences, alphabet, cost_matrix, a, b, backtracking)
        t1 = time()
        print(f"time taken for general affine function for input of size {length}: {t1-t0}")
        n.append(length)
        time_affine.append(t1-t0)

final_result = {"size": n, "time_affine": time_affine}
filehandler = open("tests/result_affine.pkl","wb")
pickle.dump(final_result,filehandler)
filehandler.close()