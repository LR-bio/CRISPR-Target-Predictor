def is_valid_pam(pam_seq, pam_pattern):
    """
    Check if PAM sequence matches PAM pattern (e.g., NGG, TTTV).
    N = any base, V = A/C/G, T = T, etc.
    """
    if len(pam_seq) != len(pam_pattern):
        return False

    pam_map = {
        "N": "ACGT",
        "V": "ACG",
        "T": "T",
        "G": "G",
        "A": "A",
        "C": "C",
    }

    for base, pattern_char in zip(pam_seq, pam_pattern):
        if base not in pam_map.get(pattern_char, ""):
            return False
    return True

def find_targets(seq, pam="NGG", grna_length=20):
    """
    Find CRISPR targets on the forward strand.
    Returns list of dicts with Start, End, Strand, gRNA, PAM.
    """
    seq = seq.upper()
    targets = []

    total_len = grna_length + len(pam)

    for i in range(len(seq) - total_len + 1):
        candidate = seq[i:i+total_len]
        grna = candidate[:grna_length]
        pam_seq = candidate[grna_length:]
        if is_valid_pam(pam_seq, pam):
            targets.append({
                "Start": i,
                "End": i + grna_length - 1,
                "Strand": "+",
                "gRNA": grna,
                "PAM": pam_seq,
            })
    return targets

def find_targets_both_strands(seq, pam="NGG", grna_length=20):
    """
    Find CRISPR targets on both strands.
    Adjust reverse strand coordinates to original sequence.
    """
    from Bio.Seq import Seq

    fwd_targets = find_targets(seq, pam=pam, grna_length=grna_length)
    rev_seq = str(Seq(seq).reverse_complement())
    rev_targets = find_targets(rev_seq, pam=pam, grna_length=grna_length)

    seq_len = len(seq)
    for t in rev_targets:
        start = seq_len - t["End"] - 1
        end = seq_len - t["Start"] - 1
        t["Start"], t["End"] = start, end
        t["Strand"] = "-"
        # gRNA on reverse strand is reverse complement of found sequence
        t["gRNA"] = str(Seq(t["gRNA"]).reverse_complement())
    return fwd_targets + rev_targets
