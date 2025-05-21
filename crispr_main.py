# crispr_main.py
# Master script for CRISPR Target Predictor with CLI and Streamlit entry points

from utils.parser import load_fasta
from utils.finder import find_targets_both_strands
from utils.offtargets import score_off_targets
from utils.scoring import calculate_gc_content, rank_targets
from utils.ml_model import predict_efficiency
import pandas as pd
import argparse
import os

def run_pipeline(fasta_path, output_path, max_mismatches=3, model_enabled=True):
    print("[+] Loading FASTA sequences...")
    records = load_fasta(fasta_path)

    all_results = []

    for rec_id, seq in records.items():
        print(f"[+] Finding CRISPR targets in: {rec_id}")
        targets = find_targets_both_strands(seq)

        for t in targets:
            t["GC%"] = calculate_gc_content(t["gRNA"])
            t["Record"] = rec_id

            # Off-target score (basic hamming-based placeholder)
            t["OffTargetHits"] = len(score_off_targets(t["gRNA"], seq, max_mismatches))

            # AI prediction (efficacy)
            if model_enabled:
                t["PredictedEfficacy"] = predict_efficiency(t["gRNA"])

        all_results.extend(targets)

    df = pd.DataFrame(all_results)

    if model_enabled:
        df = rank_targets(df)

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    df.to_csv(output_path, index=False)
    print(f"[✓] Saved results to {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Advanced CRISPR Target Predictor")
    parser.add_argument("--fasta", required=True, help="Input FASTA file")
    parser.add_argument("--out", default="reports/output.csv", help="Output CSV path")
    parser.add_argument("--no-ml", action="store_true", help="Disable ML-based prediction")
    parser.add_argument("--mismatches", type=int, default=3, help="Allowed mismatches for off-targets")

    args = parser.parse_args()
    run_pipeline(args.fasta, args.out, max_mismatches=args.mismatches, model_enabled=not args.no_ml)
