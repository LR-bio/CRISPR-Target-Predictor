def calculate_gc_content(seq):
    """
    Calculate GC content percentage of sequence.
    """
    seq = seq.upper()
    gc_count = seq.count("G") + seq.count("C")
    if len(seq) == 0:
        return 0
    return round(100 * gc_count / len(seq), 2)

def rank_targets(df):
    """
    Rank targets based on composite score from GC content, off-target hits, and predicted efficacy.
    Lower off-target hits and GC between 40-60% favored.
    """
    import numpy as np

    # Normalize columns
    df = df.copy()
    max_off = df["OffTargetHits"].max() if df["OffTargetHits"].max() > 0 else 1
    df["OffTargetScore"] = 1 - (df["OffTargetHits"] / max_off)  # higher better

    # Ideal GC range 40-60%
    df["GCScore"] = df["GC%"].apply(lambda x: 1 - abs(50 - x) / 50)  # max 1 at 50%

    # Normalize predicted efficacy between 0-1 if exists
    if "PredictedEfficacy" in df.columns and df["PredictedEfficacy"].notnull().all():
        max_eff = df["PredictedEfficacy"].max()
        min_eff = df["PredictedEfficacy"].min()
        if max_eff > min_eff:
            df["EffScore"] = (df["PredictedEfficacy"] - min_eff) / (max_eff - min_eff)
        else:
            df["EffScore"] = 0.5
    else:
        df["EffScore"] = 0.5

    # Composite rank score weighted sum
    df["RankScore"] = (0.4 * df["EffScore"]) + (0.3 * df["OffTargetScore"]) + (0.3 * df["GCScore"])
    df = df.sort_values(by="RankScore", ascending=False)
    return df
