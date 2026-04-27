import os
import io
import sys
from Bio import Phylo

file_names = os.listdir('neighbor_joining_results/')

# Get the tree results from different approaches
def get_list_of_trees(inputs_folder):
    trees_results = dict()
    for file in os.listdir(inputs_folder):
        #print(file)
        if file.endswith('.nwk'):
            file_path = os.path.join(inputs_folder, file)
            # Removing \n from rapid nj results
            with open(file_path, 'r') as f:
                content = f.read()
                if ';' in content:
                    content = content.split(';')[0] + ';'
                    #print(content)
            tree = Phylo.read(io.StringIO(content), "newick")
            trees_results[file] = tree
    return trees_results
                
NJ = get_list_of_trees('neighbor_joining_results/')
NJrapid = get_list_of_trees('neighbor_joining_results_rapid_nj/')
quicktree = get_list_of_trees('neighbor_joining_results_quicktree/')

# -------------------------------------------------------------------------------------------------
#
# rfdist.py <tree1> <tree2>
#
# Implementation of Day's algorihtm for computing the rf-distance between
# tree1 to tree2 over the same set of leaves, i.e. the number of splits not
# found in both trees.The two trees are read from the commandline and are
# assume to be in Newick-format. The Newick-parser from Biopyton
# (see https://biopython.org/wiki/Phylo) is used.
#
# Christian Storm Pedersen <cstorm@birc.au.dk>


def compare_trees(type1, type2, file_name):
    def dfs (node, splits):
        """
        performs a dfs traversal of a Phylo tree from 'node' and adds splits to the
        list 'splits' that form a consecutive interval cf. the naming of the leaves
        specified in the dictionary 'dfs_num' that maps leaf names to numbers
        """
        if node.is_terminal():
            minval = maxval = dfs_num[node.name]
            size = 1
        else:
            minval = sys.maxsize
            maxval = -sys.maxsize
            size = 0
            for child in node.clades:
                child_min, child_max, child_size = dfs(child, splits)
                if child_min < minval:
                    minval = child_min
                if child_max > maxval:
                    maxval = child_max
                size = size + child_size
            if size == maxval - minval + 1:
                splits.append((minval, maxval))
        return minval, maxval, size

    tree1 = type1[file_name]
    tree2 = type2[file_name]

    # Reroot the two trees by 'outgrouping' the same leaf.  
    root_leaf = tree1.get_terminals()[0].name
    tree1.root_with_outgroup(root_leaf)
    tree2.root_with_outgroup(root_leaf)

    # Make mapping from leaf names to dfs-numbers
    dfs_num = {}
    num = 1
    for leaf in tree1.find_clades("", True, "postorder"):
        dfs_num[leaf.name] = num
        num = num + 1

    # Collect all splits in tree1
    splits = []
    dfs(tree1.root, splits)

    # Remove the two trivial splits that occur due to the rooting of tree1
    splits = splits[:-2] 

    # Add potential shared splits form tree2
    dfs(tree2.root, splits)

    # Remove the two trivial splits that occur due to the rooting of tree2
    splits = splits[:-2]

    # The RF-distance is the number of unique splits in tree1 and tree2
    num_of_nontrivial_splits = len(tree1.get_nonterminals()) - 2 + len(tree2.get_nonterminals()) - 2
    num_of_shared_nontrivial_splits = len(splits) - len(list(set(splits)))
    rfdist = num_of_nontrivial_splits - 2 * num_of_shared_nontrivial_splits

    print(f"RF Distance for {file_name}: {rfdist}")

for file_name in file_names:
    #compare_trees(quicktree, NJ, file_name)
    #compare_trees(quicktree, NJrapid, file_name)
    compare_trees(NJrapid, NJ, file_name)

