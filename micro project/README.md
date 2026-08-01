# AI-Powered Healthcare Decision Support System for Target Protein Identification, Protein Structure Analysis, and Drug Prediction Using Protein–Protein Interaction Networks

An end-to-end computational healthcare research and decision-support prototype that accepts a disease input, constructs and analyzes Protein–Protein Interaction (PPI) networks, ranks therapeutic target proteins, fetches 3D structures (PDB/AlphaFold DB), detects binding pockets, retrieves candidate drug molecules (ChEMBL/PubChem), executes/simulates molecular docking (AutoDock Vina), computes a multi-factor computational ranking score (0–100), and presents explainable recommendations in an interactive Streamlit dashboard.

---

## ⚠️ Important Scientific Disclaimer

> **The purpose of this system is NOT to clinically prescribe medicines.** It is a computational healthcare research and decision-support prototype for identifying and analyzing potential therapeutic candidates. All scores represent computational confidence/evidence strength and require laboratory experimental and clinical validation.

---

## 🧬 Scientific Pipeline Overview

```
Disease Input (Alzheimer's, T2 Diabetes, Breast Cancer)
        ↓
Phase 1: Disease-Associated Genes & Target Proteins
        ↓
Phase 2: Protein–Protein Interaction (PPI) Network (STRING DB + NetworkX Centralities)
        ↓
Phase 3: Therapeutic Target Protein Identification & Gene Information Mapping
        ↓
Phase 4: 3D Protein Structure Retrieval (RCSB PDB / AlphaFold DB)
        ↓
Phase 5 & 6: 3D Protein Visualization & Binding Pocket Analysis (Grid Center & Box)
        ↓
Phase 7: Candidate Drug Retrieval & Molecular Processing (ChEMBL / PubChem / RDKit)
        ↓
Phase 8: Molecular Docking (AutoDock Vina / Physics Simulation Engine)
        ↓
Phase 9: Multi-Factor Computational Drug Ranking Engine (0–100 Score)
        ↓
Phase 10: Explainable AI (XAI) Evidence Synthesis & Report Generation
        ↓
Phase 11: Interactive Streamlit Dashboard
```

---

## 📊 Multi-Factor Drug Scoring Formula

$$\text{Final Score} = w_{ppi} \cdot S_{PPI} + w_{disease} \cdot S_{Disease} + w_{bio} \cdot S_{Bioactivity} + w_{dock} \cdot S_{Docking} + w_{lit} \cdot S_{Literature}$$

Where:
- **PPI Target Centrality ($S_{PPI}$, 25%)**: NetworkX Degree, Betweenness, and Eigenvector centralities.
- **Disease Association Score ($S_{Disease}$, 20%)**: Open Targets / DisGeNET association strength.
- **Drug-Target Bioactivity ($S_{Bioactivity}$, 25%)**: ChEMBL IC50 / Ki bioactivity (nM).
- **Docking Affinity ($S_{Docking}$, 20%)**: AutoDock Vina binding free energy ($\Delta G$ in $\text{kcal/mol}$).
- **Literature Support ($S_{Literature}$, 10%)**: PubMed citation density.

---

## 🚀 How to Run the Dashboard

### 1. Install Dependencies
```bash
pip install -r requirements.txt
```

### 2. Run Streamlit Application
```bash
streamlit run src/app.py
```

---

## 📁 Project Directory Structure

```
micro project/
├── data/
│   ├── cache/                  # PPI networks & API caches
│   ├── diseases/               # Disease association datasets
│   ├── proteins/               # UniProt protein metadata
│   ├── structures/             # PDB and mmCIF 3D structures
│   ├── drugs/                  # SMILES & SDF drug files
│   ├── docking/                # AutoDock Vina pose outputs & logs
│   ├── references/             # Generated research reports (.md)
│   └── database.db             # Structured SQLite database
├── config/
│   └── settings.py             # System settings & scoring weights
├── src/
│   ├── disease_analysis.py     # Phase 1: Disease-gene association
│   ├── ppi_network.py          # Phase 2: PPI graph construction
│   ├── target_ranking.py       # Phase 2: Centrality target ranking
│   ├── gene_mapping.py         # Phase 3: Gene & protein metadata
│   ├── structure_retrieval.py  # Phase 4: RCSB PDB & AlphaFold fetcher
│   ├── pocket_detection.py     # Phase 6: Binding pocket detector
│   ├── drug_retrieval.py       # Phase 7: Candidate drug retriever
│   ├── molecule_processing.py  # Phase 8: RDKit 3D molecule preparation
│   ├── docking.py              # Phase 8: AutoDock Vina runner
│   ├── drug_ranking.py         # Phase 9: Multi-factor drug ranker
│   ├── explainability.py       # Phase 10: XAI evidence report engine
│   ├── report_generator.py     # Phase 10: Report exporter
│   └── app.py                  # Phase 11: Streamlit dashboard
├── tests/                      # Unit test suites
├── requirements.txt            # Python dependencies
└── README.md
```
