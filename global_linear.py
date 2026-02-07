from Bio import SeqIO


def read_fasta(filename: str, num_sequences: int = 2):
    with open(filename) as handle:
        sequences = []
        for (i, record) in enumerate(SeqIO.parse(handle, "fasta")):
            if i < num_sequences:
                sequences.append(record.seq)
    return sequences

def global_linear(filename: str, backtracking: bool = False):
    sequences = read_fasta(filename)
    #print(sequences)
    n = len(sequences[0])
    m = len(sequences[1])
    score_matrix = [list(range(n+1))] + [[i+1] + [0] * n for i in range(m)]
    #print(score_matrix)
    #TODO: Get these from input
    gap_cost = 5
    alaphabet = ["a", "c", "g", "t"]
    cost_matrix = [[0, 5, 2, 5], [5, 0, 5, 2], [2, 5, 0, 5], [5, 2, 5, 0]]
    for i in range(1, m + 1):
        for j in range(1, n + 1):
                score_matrix[i][j] = min(score_matrix[i-1][j-1] + cost_matrix[alaphabet.index(sequences[0][j-1])][alaphabet.index(sequences[1][i-1])],
                                         score_matrix[i-1][j] + gap_cost,
                                         score_matrix[i][j-1] + gap_cost)
    print(score_matrix)

global_linear("test1.fasta")