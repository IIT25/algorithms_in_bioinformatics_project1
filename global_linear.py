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
    print(sequences)
    pass

global_linear("test1.fasta")