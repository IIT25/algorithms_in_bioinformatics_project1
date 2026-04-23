from neighbor_joining import *
from time import time
import os
import subprocess


results_folder = 'neighbor_joining_results/'
inputs_folder = 'unique_distance_matrices/'
for file in os.listdir(inputs_folder):
    print(file)
    taxa, distances = read_distance_matrix(inputs_folder + file)
    t0 = time()
    root = neighbour_joining(distances, taxa)
    t1 = time()
    print(file, t1-t0)
    newick_format = root.write_in_newick(taxa)
    write_newick_to_file(results_folder + file[:-4] + '.nwk', newick_format)
    #RapidNJ
    t0 = time()
    root =  subprocess.run(["./rapidNJ/rapidNJ-master/bin/rapidnj", (inputs_folder + file)], capture_output=True)
    t1 = time()
    print(file, " RapidNJ: ", t1-t0)

       

