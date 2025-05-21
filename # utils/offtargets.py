def hamming_distance(s1, s2):
    if len(s1) != len(s2):
        return max(len(s1), len(s2))
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def score_off_targets(gRNA, genome_seq, max_mismatches=3):
    """
    Scan the genome sequence for potential off-target sites.
    Return list of tuples: (position, sequence, mismatches).
    """
    gRNA = gRNA.upper()
    seq_len = len(gRNA)
    off_targets = []

    for i in range(len(genome_seq) - seq_len + 1):
        candidate = genome_seq[i:i+seq_len].upper()
        mismatches = hamming_distance(gRNA, candidate)
        if mismatches <= max_mismatches and mismatches > 0:
            off_targets.append((i, candidate, mismatches))

    return off_targets
