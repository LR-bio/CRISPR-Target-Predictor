def predict_efficiency(gRNA):
    """
    Placeholder ML model for predicting CRISPR gRNA cutting efficiency.
    Here we use a simple heuristic:
    - Penalize extreme GC content
    - Penalize homopolymers >4
    - Reward balanced GC around 50%
    Returns score 0 to 1.
    """
    seq = gRNA.upper()
    gc_content = (seq.count("G") + seq.count("C")) / len(seq)
    gc_score = 1 - abs(0.5 - gc_content) * 2  # max 1 at 0.5 GC

    # Check homopolymers
    max_homopolymer = max(len(s) for s in ''.join(c if c == seq[i] else ' ' for i, c in enumerate(seq)).split())
    homopolymer_penalty = 0 if max_homopolymer <= 4 else (max_homopolymer - 4) * 0.1

    score = gc_score - homopolymer_penalty
    return max(0, min(1, score))
