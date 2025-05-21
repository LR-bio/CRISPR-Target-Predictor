# 🧬 CRISPR Target Predictor

A complete and advanced web-based CRISPR target site analysis tool built with Streamlit and Python. This tool identifies gRNA target sites for Cas9, checks off-targets, computes GC content, predicts editing efficiency, and ranks candidates — all from raw DNA or FASTA input.

🚀 Features
- PAM pattern support (default: NGG for Cas9)
- Full genome scan on both strands
- Off-target mismatch scoring
- GC content analysis
- Efficiency prediction (basic ML model)
- Ranking of targets by composite score
- Upload FASTA files or enter raw sequences
- Interactive, responsive UI via Streamlit

📁 File Structure
crispr_predictor/
├── app.py # Main Streamlit UI
└── utils/
├── init.py
├── parser.py # FASTA parser
├── finder.py # gRNA + PAM site finder
├── offtargets.py # Off-target analysis
├── scoring.py # GC content, ranking
└── ml_model.py # Efficiency prediction (simple model)

⚙️ Requirements
Install the following Python packages (Python 3.8+ recommended):
--------------------------------------
pip install streamlit pandas biopython
--------------------------------------

🧑‍💻 How to Run
Clone or download this repo, then from the project root:
--------------------------------------
streamlit run app.py
Your browser will open to http://localhost:8501 with the CRISPR web app.
--------------------------------------

📥 Inputs
You can provide input in two ways:
-Paste a raw DNA sequence
-Upload a .fasta file (multi-sequence supported)

📊 Outputs
The app will display:
-Valid CRISPR targets with gRNA, PAM, strand, start/end
-GC content of each gRNA
-Off-target site count (based on mismatches)
-Predicted efficacy score
-Composite rank score to pick best candidates
-CSV download of results

⚠️ Disclaimer
This tool is for educational and research purposes only and does not replace lab validation or
professional tools like Benchling, CHOPCHOP, or CRISPOR.

🧠 Future Ideas
-Integrate with CRISPR-Cas variants (Cas12a, CasX)
-Add support for chromatin accessibility scoring
-Improve ML model with training on experimental datasets
-Visualize off-target sites on chromosomes
-Add reverse-complement gRNA alignment viewer

