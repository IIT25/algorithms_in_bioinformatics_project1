from Bio import SeqIO
from Bio.Seq import Seq           
from Bio.SeqRecord import SeqRecord
from itertools import product

def read_fasta(filename: str, num_sequences: int = 2):
    sequences = {}
    with open(filename) as handle:
        for i, record in enumerate(SeqIO.parse(handle, "fasta")):
            if i < num_sequences:
                
                sequences[record.description] = record.seq
    return sequences
       
def check_sequences_validity(sequences: list, alphabet: list):
    """
    Function to check the validity of all the sequences
    Verify that they don't have characters outside of our alphabet   
    """
    alphabet_set = set(alphabet)
    for i, sequence in enumerate(sequences):
        seq_str = str(sequence) 
        invalid_chars = set(seq_str) - alphabet_set
        if invalid_chars:
            raise Exception(f"Sequence {i+1} contains invalid characters: {invalid_chars}")
        
def read_cost_matrix(filename: str):
    """
    Function to read a file with the cost matrix and the gapcost, 
    The input file must be a .txt file formated as follow:
    1 line:
    n g
    where n is the size of the alphabet and g is the linear gapcost

    next n lines:
    letter cost1 cost2 ... costn
    
    Where the letter is the i'th letter of the alphabet and each row defines the cost of a substitution
    of alphabet[i] to the rest of the alphabet
    Notice that there is just one space between values of the same line
    """
    with open(filename) as f:
        lines = f.readlines()
        n, gapcost = map(int,lines[0].split(' '))
        alphabet = []
        cost_matrix = []
        for i in range(1, n+1):
            line = lines[i].split(' ')
            alphabet.append(line[0])
            cost_matrix.append(list(map(int, line[1:])))
    return gapcost, alphabet, cost_matrix
    
def write_alignment_in_fasta(alignment: list, filename: str, sequence_names: list = [], comment: str = ""):
    """
    Function to write the output alignment in a fasta file. Receives:
        - alignment: a list with the two sequences already aligned
        - filename: the name of the file where we will write our output
        - sequence_names: A list with the names of the sequences to be inserted in the file
    """
    if sequence_names == []:
        records = [SeqRecord(Seq(alignment[i]), id=f"seq{i}", description="align") for i in range(len(alignment))]
    else:
        records = [SeqRecord(Seq(alignment[i]), id=sequence_names[i], description="align") for i in range(len(alignment))]
    with open(filename, "w") as f:
        if comment:
            f.write(f"; {comment}\n")
        SeqIO.write(records, f, "fasta")
