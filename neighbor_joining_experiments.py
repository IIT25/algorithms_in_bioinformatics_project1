from neighbor_joining import *
from time import time
import os
import subprocess
import pandas as pd

results_folder = 'neighbor_joining_results/'
inputs_folder = 'unique_distance_matrices/'
quicktree_results = 'neighbor_joining_results_quicktree/'
nj_times = []
rapid_nj_times = []
quicktree_times = []
lengths = []

for file in os.listdir(inputs_folder):
    print(file)
    taxa, distances = read_distance_matrix(inputs_folder + file)
    input_path = os.path.join(inputs_folder, file)
    lengths.append(len(taxa))
    ### Our implementation
    t0 = time()
    root = neighbour_joining(distances, taxa)
    t1 = time()
    nj_times.append(t1-t0)
    print(file, t1-t0)
    newick_format = root.write_in_newick(taxa)
    write_newick_to_file(results_folder + file[:-4] + '.nwk', newick_format)

    ### RapidNJ
    t0 = time()
    rapid_nj_out =  subprocess.run(["./rapidNJ/rapidNJ-master/bin/rapidnj", (inputs_folder + file)], capture_output=True)
    t1 = time()
    rapid_nj_times.append(t1-t0)
    print(file, " RapidNJ: ", t1-t0)
    write_newick_to_file(results_folder[:-1] + '_rapid_nj/' + file[:-4] + '.nwk', str(rapid_nj_out.stdout))

    ### QuickTree
    t0 = time()
    quicktree_out = subprocess.run(["./quicktree/quicktree", "-in", "m", input_path],
                                   capture_output=True,
                                   text=True)
    t1 = time()
    print(f"{file} QuickTree time: {t1-t0:.4f}s")
    quicktree_times.append(t1-t0)
    output_filename = file[:-4] + '.nwk'
    output_path = os.path.join(quicktree_results, output_filename)
    write_newick_to_file(output_path, quicktree_out.stdout)

time_results = pd.DataFrame({"n": lengths, "neighbor_joining": nj_times, "rapid_nj": rapid_nj_times, "quicktree": quicktree_times})
time_results.to_csv("time_experiments_results.csv")

# Heello