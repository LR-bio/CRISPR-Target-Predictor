import streamlit as st
import pandas as pd
from Bio import Entrez
from Bio.Seq import Seq
from utils.parser import load_fasta
from utils.finder import find_targets_both_strands
from utils.offtargets import score_off_targets
from utils.scoring import calculate_gc_content, rank_targets
from utils.ml_model import predict_efficiency
import plotly.express as px

st.set_page_config(layout="wide")
st.title("🔬 Advanced CRISPR Target Predictor")

# --- Enzyme & PAM support ---
CAS_ENZYMES = {
    "Cas9 (NGG)": {"pam": "NGG", "grna_len": 20},
    "Cas12a (TTTV)": {"pam": "TTTV", "grna_len": 23},
}

Entrez.email = "your_email@example.com"  # Replace with your email for NCBI queries

def fetch_sequence_from_ncbi(gene_name):
    try:
        handle = Entrez.esearch(db="nucleotide", term=gene_name, retmax=1)
        record = Entrez.read(handle)
        if record["IdList"]:
            seq_id = record["IdList"][0]
            handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta", retmode="text")
            fasta = handle.read()
            return fasta
        else:
            return None
    except Exception as e:
        st.error(f"NCBI fetch error: {e}")
        return None

def parse_fasta_string(fasta_str):
    lines = fasta_str.strip().split("\n")
    header = lines[0][1:]
    seq = "".join(lines[1:]).upper()
    return {header: seq}

# UI Elements
enzyme = st.selectbox("Select CRISPR enzyme", list(CAS_ENZYMES.keys()))
pam = CAS_ENZYMES[enzyme]["pam"]
grna_length = CAS_ENZYMES[enzyme]["grna_len"]

input_mode = st.radio("Input sequence or fetch by gene name?", ["Upload FASTA file", "Fetch from NCBI by gene"])

max_mismatches = st.slider("Max mismatches for off-target scanning", 0, 5, 3)
run_ml = st.checkbox("Enable efficiency prediction (ML model)", value=True)

sequence_data = {}

if input_mode == "Upload FASTA file":
    uploaded_file = st.file_uploader("Upload FASTA file", type=["fasta", "fa", "fna"])
    if uploaded_file:
        sequence_data = load_fasta(uploaded_file)

else:
    gene_name = st.text_input("Enter gene or accession to fetch from NCBI")
    if gene_name:
        with st.spinner("Fetching sequence from NCBI..."):
            fasta_str = fetch_sequence_from_ncbi(gene_name)
            if fasta_str:
                sequence_data = parse_fasta_string(fasta_str)
                st.success(f"Fetched sequence for {gene_name}")
            else:
                st.error("No sequence found.")

if sequence_data:
    for rec_id, seq in sequence_data.items():
        st.write(f"## Sequence: {rec_id} (Length: {len(seq)})")

        targets = find_targets_both_strands(seq, pam=pam, grna_length=grna_length)

        for t in targets:
            t["GC%"] = calculate_gc_content(t["gRNA"])
            t["Record"] = rec_id
            off_targets = score_off_targets(t["gRNA"], seq, max_mismatches)
            t["OffTargetHits"] = len(off_targets)
            t["OffTargetDetails"] = off_targets
            if run_ml:
                t["PredictedEfficacy"] = predict_efficiency(t["gRNA"])
            else:
                t["PredictedEfficacy"] = None

        df = pd.DataFrame(targets)

        if run_ml and not df.empty:
            df = rank_targets(df)

        if df.empty:
            st.warning("No CRISPR targets found.")
            continue

        tab1, tab2, tab3 = st.tabs(["Targets Table", "Off-target Details", "Summary & Visualization"])

        with tab1:
            st.write("### CRISPR Targets")
            display_cols = ["Start", "End", "Strand", "gRNA", "PAM", "GC%", "OffTargetHits", "PredictedEfficacy"]
            st.dataframe(df[display_cols].sort_values(by="RankScore" if "RankScore" in df.columns else "PredictedEfficacy", ascending=False))

            csv = df.to_csv(index=False).encode("utf-8")
            st.download_button(
                label="Download targets as CSV",
                data=csv,
                file_name=f"{rec_id}_crispr_targets.csv",
                mime="text/csv"
            )

        with tab2:
            st.write("### Off-target Details")
            for idx, row in df.iterrows():
                st.write(f"#### gRNA: {row['gRNA']} (Start: {row['Start']}, Strand: {row['Strand']})")
                off_targets = row["OffTargetDetails"]
                if off_targets:
                    off_df = pd.DataFrame(off_targets, columns=["Position", "Sequence", "Mismatches"])
                    st.dataframe(off_df)
                else:
                    st.write("No off-targets within mismatch threshold.")

        with tab3:
            st.write("### Summary Charts")

            counts = df["Strand"].value_counts().rename_axis("Strand").reset_index(name="Counts")
            fig1 = px.bar(counts, x="Strand", y="Counts", title="Target Count by Strand")
            st.plotly_chart(fig1, use_container_width=True)

            fig2 = px.histogram(df, x="GC%", nbins=20, title="GC% Distribution")
            st.plotly_chart(fig2, use_container_width=True)

            if run_ml:
                fig3 = px.histogram(df, x="PredictedEfficacy", nbins=20, title="Predicted Efficiency Distribution")
                st.plotly_chart(fig3, use_container_width=True)

            st.write("### Genome Browser - Target Positions")
            browser_df = df[["Start", "Strand", "gRNA"]].copy()
            browser_df["y"] = 1  # dummy y-axis
            fig4 = px.scatter(browser_df, x="Start", y="y", color="Strand",
                              hover_data=["gRNA"],
                              labels={"x": "Genome Position", "y": ""},
                              title="Target sites along sequence")
            fig4.update_yaxes(showticklabels=False)
            st.plotly_chart(fig4, use_container_width=True)
