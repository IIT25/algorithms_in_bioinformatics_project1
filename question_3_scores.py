from msa_sp_score_3k import *
import matplotlib.pyplot as plt
from sp_exact_3 import *

# approximation_scores = [73, 144, 265, 334, 430, 500, 575, 650, 727, 756, 819, 936, 987, 1041, 1112, 1145, 1308, 1264, 1450, 1498]
# exact_stores = [70, 135, 231, 318, 385, 440, 516, 589, 628, 687, 754, 810, 895, 957, 1023, 1080, 1186, 1158, 1323, 1379]

if __name__ == "__main__":
    eval_directory = "testseqs/"
    cost_matrix_filename = 'cost_matrix.txt'
    gapcost, alphabet, cost_matrix = read_cost_matrix(eval_directory + cost_matrix_filename)
    
    exact_scores = []
    approximation_scores = []
    sizes = range(10,210,10)

    for i in sizes:
        approximation_scores.append(compute_sp_score(eval_directory + f"approximation_output_testseqs_{i}_3.fasta"))

        exact_scores.append(compute_sp_score(eval_directory + f"exact_output_testseqs_{i}_3.fasta"))
        
    print(exact_scores)
    print(approximation_scores)
    
    plt.plot(sizes,np.array(approximation_scores) / np.array(exact_scores))
    plt.show()
      
    