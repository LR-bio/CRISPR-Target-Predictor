# utils/offtargets.py

def hamming_distance(s1, s2):
    """
    Calculate the Hamming distance between two equal-length strings.
    """
    if len(s1) != len(s2):
        return max(len(s1), len(s2))  # Penalize different lengths heavily
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def score_off_targets(gRNA, reference_seq, max_mismatches=3):
    """
    Scan reference sequence for potential off-target sites
    allowing up to max_mismatches mismatches.
    Returns list of tuples: (position, off_target_seq, mismatch_count)
    """
    gRNA = gRNA.upper()
    seq_len = len(gRNA)
    off_targets = []

    for i in range(len(reference_seq) - seq_len + 1):
        candidate = reference_seq[i:i + seq_len].upper()
        dist = hamming_distance(gRNA, candidate)
        if dist <= max_mismatches:
            off_targets.append((i, candidate, dist))

    return off_targets
