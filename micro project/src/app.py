import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import json
from pathlib import Path
import sys

# Ensure root directory is in sys.path
sys.path.append(str(Path(__file__).resolve().parent.parent))

from config.settings import SUPPORTED_DISEASES, STRUCTURE_DIR, DRUG_DIR, REFERENCE_DIR
from src.disease_analysis import DiseaseAnalyzer
from src.ppi_network import PPINetworkBuilder
from src.target_ranking import TargetRanker
from src.gene_mapping import GeneProteinMapper
from src.structure_retrieval import StructureRetriever
from src.pocket_detection import PocketDetector
from src.drug_retrieval import DrugRetriever
from src.molecule_processing import MoleculeProcessor
from src.docking import DockingEngine
from src.drug_ranking import DrugRanker
from src.explainability import ExplainabilityEngine
from src.report_generator import ReportGenerator

# Page Configuration
st.set_page_config(
    page_title="AI Healthcare Decision Support - Target & Drug Prediction",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom Styling (Glassmorphism & Modern Medical Palette)
st.markdown("""
<style>
    .main-header {
        font-size: 2.2rem;
        font-weight: 700;
        background: linear-gradient(135deg, #0284c7 0%, #0d9488 100%);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        margin-bottom: 0.5rem;
    }
    .sub-header {
        font-size: 1.1rem;
        color: #64748b;
        margin-bottom: 1.5rem;
    }
    .metric-card {
        background: #f8fafc;
        border: 1px solid #e2e8f0;
        border-radius: 10px;
        padding: 1.2rem;
        text-align: center;
        box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.05);
    }
    .disclaimer-box {
        background-color: #fffbe6;
        border-left: 5px solid #d97706;
        padding: 1rem;
        border-radius: 5px;
        font-size: 0.9rem;
        color: #92400e;
        margin-top: 1rem;
        margin-bottom: 1rem;
    }
    .badge-experimental {
        background-color: #dcfce7;
        color: #166534;
        padding: 3px 10px;
        border-radius: 12px;
        font-weight: 600;
        font-size: 0.85rem;
    }
    .badge-alphafold {
        background-color: #e0f2fe;
        color: #075985;
        padding: 3px 10px;
        border-radius: 12px;
        font-weight: 600;
        font-size: 0.85rem;
    }
</style>
""", unsafe_allow_html=True)

# Application Header
st.markdown('<div class="main-header">AI-Powered Healthcare Decision Support System</div>', unsafe_allow_html=True)
st.markdown('<div class="sub-header">Target Protein Identification, 3D Structure Analysis & Molecular Docking Drug Prediction</div>', unsafe_allow_html=True)

# Mandatory Disclaimer
st.markdown("""
<div class="disclaimer-box">
    <strong>⚠️ Scientific & Decision Support Disclaimer:</strong> 
    This platform is a computational research and decision-support prototype for target protein identification and candidate drug analysis. 
    Predicted targets, molecular interactions, and candidate drug scores are computational confidence metrics and require experimental and clinical validation.
</div>
""", unsafe_allow_html=True)

# Sidebar Disease Selector
st.sidebar.header("🔬 Pipeline Configuration")
disease_option = st.sidebar.selectbox(
    "Select Input Disease",
    options=list(SUPPORTED_DISEASES.keys()),
    format_func=lambda x: SUPPORTED_DISEASES[x]
)

disease_label = SUPPORTED_DISEASES[disease_option]
st.sidebar.success(f"Active Disease Model: **{disease_label}**")

# Instantiate Services
disease_analyzer = DiseaseAnalyzer()
target_ranker = TargetRanker()
gene_mapper = GeneProteinMapper()
struct_retriever = StructureRetriever()
pocket_detector = PocketDetector()
drug_retriever = DrugRetriever()
molecule_processor = MoleculeProcessor()
docking_engine = DockingEngine()
drug_ranker_engine = DrugRanker()
explainability_engine = ExplainabilityEngine()
report_generator = ReportGenerator()

# Load Disease Associations & Targets
try:
    disease_data = disease_analyzer.fetch_disease_associations(disease_option)
    top_targets_df, ppi_graph = target_ranker.rank_target_proteins(disease_option, top_n=6)
except Exception as e:
    st.error(f"Error initializing pipeline data: {e}")
    st.stop()

# Dashboard Tabs
tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "📊 1. Disease & PPI Network",
    "🧬 2. Target Proteins & Ranking",
    "🔮 3. 3D Protein Structure & Pocket",
    "🧪 4. Molecular Docking & Binding Pose",
    "💡 5. Drug Prediction & XAI Report"
])

# ==================== TAB 1: DISEASE & PPI NETWORK ====================
with tab1:
    st.subheader(f"Disease Genomics & PPI Network Analysis: {disease_label}")
    st.write(disease_data.get("description", ""))

    col1, col2 = st.columns([1, 2])
    
    with col1:
        st.markdown("### Disease Associated Genes")
        df_assoc = disease_analyzer.get_association_dataframe(disease_option)
        st.dataframe(
            df_assoc[["gene_symbol", "uniprot_id", "gene_name", "association_score"]],
            use_container_width=True,
            height=380
        )

    with col2:
        st.markdown("### Interactive Protein-Protein Interaction (PPI) Network")
        if ppi_graph.number_of_nodes() > 0:
            # Generate 2D Network Plot with Plotly
            import networkx as nx
            pos = nx.spring_layout(ppi_graph, seed=42)
            
            edge_x, edge_y = [], []
            for edge in ppi_graph.edges():
                x0, y0 = pos[edge[0]]
                x1, y1 = pos[edge[1]]
                edge_x.extend([x0, x1, None])
                edge_y.extend([y0, y1, None])

            edge_trace = go.Scatter(
                x=edge_x, y=edge_y,
                line=dict(width=1.5, color='#94a3b8'),
                hoverinfo='none',
                mode='lines'
            )

            node_x, node_y, node_text, node_size = [], [], [], []
            for node in ppi_graph.nodes():
                x, y = pos[node]
                node_x.append(x)
                node_y.append(y)
                node_text.append(f"Protein: {node}")
                deg = ppi_graph.degree(node)
                node_size.append(25 + deg * 8)

            node_trace = go.Scatter(
                x=node_x, y=node_y,
                mode='markers+text',
                hoverinfo='text',
                text=[node for node in ppi_graph.nodes()],
                textposition="top center",
                marker=dict(
                    showscale=True,
                    colorscale='YlGnBu',
                    color=node_size,
                    size=node_size,
                    colorbar=dict(title='Network Connectivity'),
                    line_width=2
                )
            )

            fig = go.Figure(data=[edge_trace, node_trace],
                 layout=go.Layout(
                    showlegend=False,
                    hovermode='closest',
                    margin=dict(b=20,l=5,r=5,t=20),
                    xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                    height=400
                ))
            st.plotly_chart(fig, use_container_width=True)

# ==================== TAB 2: TARGET RANKING ====================
with tab2:
    st.subheader("Therapeutic Target Protein Ranking (PPI Network Centrality)")
    st.write("Proteins are ranked based on Network Centrality (Degree, Betweenness, Eigenvector) and Disease Association Evidence.")

    st.dataframe(
        top_targets_df[["gene_symbol", "uniprot_id", "gene_name", "degree_centrality", "betweenness_centrality", "association_score", "target_importance_score"]].style.highlight_max(axis=0, color="#dcfce7"),
        use_container_width=True
    )

    # Bar chart of top targets
    fig_bar = px.bar(
        top_targets_df,
        x="gene_symbol",
        y="target_importance_score",
        color="association_score",
        title="Ranked Target Importance Score (PPI Centrality + Disease Association)",
        labels={"target_importance_score": "Composite Target Score", "gene_symbol": "Gene Symbol"},
        color_continuous_scale="Viridis"
    )
    st.plotly_chart(fig_bar, use_container_width=True)

# ==================== TAB 3: 3D PROTEIN STRUCTURE & POCKET ====================
with tab3:
    st.subheader("Target Protein 3D Structure & Binding Pocket Identification")

    target_genes = top_targets_df["gene_symbol"].tolist()
    selected_target = st.selectbox("Select Target Protein for 3D Analysis", options=target_genes, index=0)

    target_row = top_targets_df[top_targets_df["gene_symbol"] == selected_target].iloc[0]
    uniprot_id = target_row["uniprot_id"]

    # Fetch Detailed Protein Info & 3D Structure
    protein_info = gene_mapper.get_protein_info(uniprot_id, selected_target)
    struct_file, struct_meta = struct_retriever.get_structure(selected_target, uniprot_id)
    pockets = pocket_detector.detect_pockets(struct_file, selected_target)
    primary_pocket = pocket_detector.get_primary_pocket(struct_file, selected_target)

    c1, c2 = st.columns([1, 1])

    with c1:
        st.markdown(f"### Protein Info: {protein_info.get('protein_name')}")
        st.markdown(f"- **Gene Symbol:** `{selected_target}`")
        st.markdown(f"- **UniProt ID:** `{uniprot_id}`")
        st.markdown(f"- **Molecular Function:** {protein_info.get('molecular_function')}")
        st.markdown(f"- **Subcellular Location:** {protein_info.get('subcellular_location')}")
        st.markdown(f"- **Disease Role:** {protein_info.get('disease_association')}")

        st.markdown("---")
        st.markdown("### 3D Structure Provenance")
        source = struct_meta.get("structure_source", "")
        if "RCSB" in source:
            st.markdown(f'<span class="badge-experimental">Experimental Structure ({struct_meta.get("pdb_id")})</span>', unsafe_allow_html=True)
        else:
            st.markdown(f'<span class="badge-alphafold">AI-Predicted Structure ({struct_meta.get("pdb_id")})</span>', unsafe_allow_html=True)
        st.write(f"**Method:** {struct_meta.get('method')} | **Resolution/pLDDT:** {struct_meta.get('resolution')}")

    with c2:
        st.markdown("### Detected Binding Pocket Coordinates")
        st.write(f"**Pocket Name:** {primary_pocket.get('name')}")
        st.write(f"**Center (X, Y, Z):** `({primary_pocket['center']['x']}, {primary_pocket['center']['y']}, {primary_pocket['center']['z']})`")
        st.write(f"**Grid Size (Å):** `({primary_pocket['size']['x']}, {primary_pocket['size']['y']}, {primary_pocket['size']['z']})`")
        st.write(f"**Key Catalytic Residues:** `{', '.join(primary_pocket.get('key_residues', []))}`")

        # Display PDB file preview
        with open(struct_file, "r", encoding="utf-8") as pf:
            pdb_preview = pf.readlines()[:15]
        st.text_area("PDB Coordinate File Preview (Head)", "".join(pdb_preview), height=150)

# ==================== TAB 4: MOLECULAR DOCKING ====================
with tab4:
    st.subheader("Molecular Docking (AutoDock Vina Simulation)")

    candidate_drugs = drug_retriever.fetch_candidate_drugs(selected_target, disease_option)
    drug_names = [d["drug_name"] for d in candidate_drugs]
    
    selected_drug_name = st.selectbox("Select Candidate Drug Molecule for Docking", options=drug_names, index=0)
    selected_drug = next(d for d in candidate_drugs if d["drug_name"] == selected_drug_name)

    # Process molecule & run docking
    mol_info = molecule_processor.process_smiles(selected_drug["smiles"], selected_drug_name)
    dock_res = docking_engine.run_docking(selected_target, uniprot_id, selected_drug, primary_pocket, struct_file)

    dc1, dc2 = st.columns([1, 1])

    with dc1:
        st.markdown(f"### Candidate Drug Properties: {selected_drug_name}")
        st.markdown(f"- **ChEMBL ID:** `{selected_drug.get('chembl_id')}`")
        st.markdown(f"- **SMILES:** `{selected_drug.get('smiles')}`")
        st.markdown(f"- **Molecular Weight:** `{mol_info.get('mw')} g/mol`")
        st.markdown(f"- **LogP:** `{mol_info.get('logp')}`")
        st.markdown(f"- **H-Donors / H-Acceptors:** `{mol_info.get('hbd')} / {mol_info.get('hba')}`")
        st.markdown(f"- **Lipinski's Rule of 5:** `{mol_info.get('lipinski_rule_of_5')}`")

    with dc2:
        st.markdown("### Docking Simulation Output")
        st.metric("Predicted Binding Affinity (ΔG)", f"{dock_res.get('binding_affinity_kcal_mol')} kcal/mol")
        st.metric("Docking Pose Confidence Score", f"{dock_res.get('docking_score_confidence')}/100")
        st.write(f"**Docking Engine:** {dock_res.get('engine')}")
        
        st.markdown("#### Hydrogen Bonding Contacts:")
        for hb in dock_res.get("hydrogen_bonds", []):
            st.markdown(f"- Receptor Residue `{hb['receptor_residue']}` ↔ Ligand `{hb['ligand_atom']}` ({hb['distance_angstrom']} Å)")

# ==================== TAB 5: DRUG PREDICTION & XAI REPORT ====================
with tab5:
    st.subheader("Explainable Recommendation Engine & Research Report")

    # Dock all candidates for target
    all_dock_results = [
        docking_engine.run_docking(selected_target, uniprot_id, d, primary_pocket, struct_file)
        for d in candidate_drugs
    ]

    target_dict = {
        "gene_symbol": selected_target,
        "uniprot_id": uniprot_id,
        "protein_name": protein_info.get("protein_name"),
        "ppi_centrality_score": target_row.get("ppi_centrality_score", 0.90),
        "association_score": target_row.get("association_score", 0.90),
        "disease_role": protein_info.get("disease_association")
    }

    rankings_df = drug_ranker_engine.rank_candidate_list(target_dict, candidate_drugs, all_dock_results)

    st.markdown(f"### Final Candidate Drug Rankings for Target `{selected_target}`")
    st.dataframe(
        rankings_df[["rank", "drug_name", "computational_ranking_score", "docking_affinity_kcal_mol", "bioactivity_ic50_nm", "classification_label"]].style.highlight_max(subset=["computational_ranking_score"], color="#dcfce7"),
        use_container_width=True
    )

    # Top Candidate Explanation
    top_rank_record = rankings_df.iloc[0]
    top_drug_data = next(d for d in candidate_drugs if d["drug_name"] == top_rank_record["drug_name"])
    top_dock_res = next(r for r in all_dock_results if r["drug_name"] == top_rank_record["drug_name"])

    # Generate Explanation
    exp_output = explainability_engine.generate_explanation(
        disease_label,
        target_dict,
        struct_meta,
        primary_pocket,
        top_drug_data,
        top_dock_res,
        top_rank_record.to_dict()
    )

    st.markdown("---")
    st.markdown(f"### 📑 Explainable AI (XAI) Evidence Synthesis: {top_rank_record['drug_name']}")
    st.markdown(exp_output["narrative_explanation"])

    # Generate Downloadable Report
    md_content, r_path = report_generator.generate_markdown_report(exp_output)
    
    st.download_button(
        label="📥 Download Complete Healthcare Decision Support Report (.md)",
        data=md_content,
        file_name=f"Report_{disease_option}_{selected_target}.md",
        mime="text/markdown"
    )

st.markdown("---")
st.caption("AI-Powered Healthcare Decision Support System | Computer Science Capstone Prototype | Version 1.0")
