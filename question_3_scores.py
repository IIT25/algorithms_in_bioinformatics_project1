from msa_sp_score_3k import *
import matplotlib.pyplot as plt

"""
TODO
Add the exact alignment values to plot the ratios
"""


if __name__ == "__main__":
    alignments_directory = "testseqs/"
    

    approximation_scores = []
    sizes = range(10,210,10)

    for i in sizes:
        approximation_scores.append(compute_sp_score(alignments_directory + f"approximation_output_testseqs_{i}_3.fasta"))
    
    plt.plot(sizes,approximation_scores)
    plt.show()
    print(approximation_scores)