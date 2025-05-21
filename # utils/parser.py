# utils/parser.py
from Bio import SeqIO

def load_fasta(fasta_file):
    """
    Loads a FASTA file and returns a dictionary of {record_id: sequence}.
    """
    records = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        records[record.id] = str(record.seq).upper()
    return records
