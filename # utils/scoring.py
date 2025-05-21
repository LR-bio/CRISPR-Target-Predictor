# utils/scoring.py
import pandas as pd

def calculate_gc_content(seq):
    gc = seq.count("G") + seq.count("C")
    return 100 * gc / len(seq) if len(seq) > 0 else 0

def rank_targets(df):
    """
    Rank targets based on:
    - Number of off-target hits (fewer is better)
    - Predicted efficiency (higher is better, if available)
    Returns DataFrame sorted by rank_score.
    """
    # If no predicted efficacy, fill with 0
    if "PredictedEfficacy" not in df.columns:
        df["PredictedEfficacy"] = 0

    # Normalize off-target hits and efficacy for ranking
    df["OffTargetHitsNorm"] = (df["OffTargetHits"] - df["OffTargetHits"].min()) / (df["OffTargetHits"].max() - df["OffTargetHits"].min() + 1e-6)
    df["PredictedEfficacyNorm"] = (df["PredictedEfficacy"] - df["PredictedEfficacy"].min()) / (df["PredictedEfficacy"].max() - df["PredictedEfficacy"].min() + 1e-6)

    # Composite rank score: lower off-target and higher efficacy preferred
    df["RankScore"] = df["PredictedEfficacyNorm"] - df["OffTargetHitsNorm"]

    return df.sort_values(by="RankScore", ascending=False).reset_index(drop=True)
