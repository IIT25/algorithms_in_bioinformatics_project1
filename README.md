Project 1: Pairwise global alignment with linear gap cost

Group 10

Emma Réka Erdei

Martina Chabadová

Juan Nicolás Quintero Quintero

Introduction

We have implemented global optimal alignment using linear gap cost with a backtracking algorithm and a global count function to find the total number of optimal paths. Everything works as expected, and our results match the provided test cases.

Methods

To implement our solution, we use three different files, one as a helper and two containing the main functions to align the sequences and count the possible optimal alignments:

Helper Functions (alignment_utils.py)

read_fasta: Reads sequences from FASTA files using the Bio.SeqIO library.

read_cost_matrix: Reads text file with the cost matrix and the gapcost formated as following: n g where n is the size of the alphabet and g is the linear gapcost

next n lines: letter cost1 cost2 ... costn

Where the letter is the i'th letter of the alphabet and each row defines the cost of a substitution of alphabet[i] to the rest of the alphabet

check_sequences_validity: Checks the validity of all the sequences, verifies that they don't have characters outside of our alphabet

write_alignment_in_fasta: Writes the output alignment in a fasta file.

Global Linear Alignment (global_linear.py)

Function to get the global pairwise alignment using linear gap cost It receives:

    - sequences: a list of two sequences that we want to align
    - alphabet: a list of the current alphabet, aligned with the index of the cost matrix
    - cost_matrix: a matrix with the scores of each substitution
    - gap_cost: an integer representing the gap cost on the alignment
    - backtracking: a boolean variable to return an optimal alignment using backtracking, if desired 
We implemented the linear optimal cost finding algorithm using an iterative Dynamic Programming (DP) approach, where we fill the matrix S that we saw in class iteratively, row by row to avoid hitting max recursion limits.

Initialization: We create a DP table of dimensions . The first row and column are initialized with gap costs scaled by their index (), this reduces the number of iterations and comparisons, followed by the observation that the alignments in the first row and in the first column can only come from a succession of insertions or deletions.
Implementation: The table is filled iteratively (top-left to bottom-right). We chose an iterative approach over recursion to avoid maximum recursion depth issues with long sequences.
Backtracking: If requested, a separate function traverses the score matrix from  back to  to reconstruct the optimal alignment. Opposite to the global score algorithm, we decided to implement the bactracking function with a recursive approach, even if an iterative approach is possible, the recursive approach was more straightforward to implement.
Global Count (global_count.py)

Computes the total number of optimal alignments using a count_matrix with the score_matrix calculated in global_linear.

Logic: We set the first row and column of the count matrix to 1, this is because, as for the score matrix calculation, the alignments in the first row and first column can only come from a series of insertions or deletions.
Calculation: For each cell, the value is the sum of the count values of its neighbors (diagonal, top, or left) * if and only only if* that neighbor provided the optimal score for the current cell in the score matrix.
Result: The final value is located at the bottom-right corner of the count matrix.
