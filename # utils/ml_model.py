# utils/ml_model.py

def predict_efficiency(gRNA):
    """
    Placeholder function for gRNA efficiency prediction.
    Returns a score between 0 and 1.
    Uses simple heuristics on GC content and nucleotide position preferences.
    Replace this stub with a real ML model as needed.
    """
    gRNA = gRNA.upper()

    # GC content contribution (ideal between 40-60%)
    gc = (gRNA.count("G") + gRNA.count("C")) / len(gRNA)
    gc_score = 1 - abs(gc - 0.5) * 2  # peaks at 0.5, decreases to 0 at extremes

    # Positional nucleotide preferences (made-up weights)
    weights = [0.1 if nt in "GC" else 0.05 for nt in gRNA]

    position_score = sum(weights) / len(gRNA)

    score = (gc_score + position_score) / 2
    return max(0, min(1, score))  # Clamp between 0 and 1
