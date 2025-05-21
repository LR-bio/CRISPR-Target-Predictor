def load_fasta(file):
    """
    Parses a FASTA file-like object and returns a dictionary of sequences.
    Key = header line without '>', value = sequence string (uppercase).
    """
    sequences = {}
    header = None
    seq_lines = []

    for line in file:
        line = line.decode("utf-8") if isinstance(line, bytes) else line
        line = line.strip()
        if line.startswith(">"):
            if header:
                sequences[header] = "".join(seq_lines).upper()
            header = line[1:]
            seq_lines = []
        else:
            seq_lines.append(line)
    if header:
        sequences[header] = "".join(seq_lines).upper()

    return sequences
