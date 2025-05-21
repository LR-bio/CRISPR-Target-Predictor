# utils/finder.py
from Bio.Seq import Seq

# Constants for PAM recognition; extendable for other Cas systems
PAMS = {
    "Cas9": "NGG",
    # Add other PAM types here if needed, e.g. Cas12a ("TTTV"), etc.
}

def is_valid_pam(pam_seq, pam="NGG"):
    """
    Checks if the PAM sequence matches the PAM pattern (N = any base).
    Currently supports Cas9 NGG PAM.
    """
    if pam == "NGG":
        return len(pam_seq) == 3 and pam_seq[1:] == "GG" and pam_seq[0] in "ACGT"
    # Add other PAM patterns here
    return False

def calculate_gc_content(seq):
    gc = seq.count("G") + seq.count("C")
    return 100 * gc / len(seq) if len(seq) > 0 else 0

def find_targets_in_seq(seq, pam="NGG", grna_length=20):
    """
    Find valid CRISPR targets (gRNA + PAM) on a single DNA strand.
    Returns list of dicts with position, gRNA, PAM, and strand.
    """
    targets = []
    window_size = grna_length + len(pam)

    for i in range(len(seq) - window_size + 1):
        candidate = seq[i:i + window_size]
        gRNA = candidate[:grna_length]
        pam_seq = candidate[grna_length:]
        if is_valid_pam(pam_seq, pam):
            targets.append({
                "Start": i,
                "End": i + window_size,
                "gRNA": gRNA,
                "PAM": pam_seq,
                "Strand": "+"
            })
    return targets

def find_targets_both_strands(seq, pam="NGG", grna_length=20):
    """
    Finds targets on both strands.
    """
    targets = []

    # Forward strand
    fwd_targets = find_targets_in_seq(seq, pam, grna_length)
    targets.extend(fwd_targets)

    # Reverse strand
    rev_seq = str(Seq(seq).reverse_complement())
    rev_targets = find_targets_in_seq(rev_seq, pam, grna_length)
    seq_len = len(seq)
    # Adjust positions to original sequence coordinates
    for t in rev_targets:
        start = seq_len - t["End"]
        end = seq_len - t["Start"]
        targets.append({
            "Start": start,
            "End": end,
            "gRNA": t["gRNA"],
            "PAM": t["PAM"],
            "Strand": "-"
        })

    return targets
