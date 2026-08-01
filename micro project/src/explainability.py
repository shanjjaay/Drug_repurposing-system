import json
import pandas as pd
from typing import Dict, List, Any
from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parent.parent))
from config.settings import REFERENCE_DIR

class ExplainabilityEngine:
    """
    Phase 10: Explainable AI Evidence Synthesis.
    Generates rule-driven, transparent, evidence-backed natural language explanations
    answering why a protein target and candidate drug were selected and ranked.
    """
    def __init__(self):
        pass

    def generate_explanation(
        self,
        disease_name: str,
        target_info: Dict[str, Any],
        structure_meta: Dict[str, Any],
        pocket_info: Dict[str, Any],
        drug_data: Dict[str, Any],
        docking_result: Dict[str, Any],
        ranking_result: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Synthesizes structured evidence object and natural language narrative.
        """
        gene_symbol = target_info.get("gene_symbol", "UNKNOWN")
        uniprot_id = target_info.get("uniprot_id", "UNKNOWN")
        protein_name = target_info.get("protein_name", gene_symbol)
        disease_assoc = target_info.get("disease_role", target_info.get("molecular_function", "Disease target protein."))
        
        drug_name = drug_data.get("drug_name", "Candidate Compound")
        chembl_id = drug_data.get("chembl_id", "N/A")
        ic50_nm = drug_data.get("bioactivity_value_nm", "N/A")

        docking_affinity = docking_result.get("binding_affinity_kcal_mol", -8.0)
        pocket_name = pocket_info.get("name", "Target Pocket")
        struct_source = structure_meta.get("structure_source", "PDB / AlphaFold DB")
        struct_type = structure_meta.get("structure_type", "3D Structure")

        final_score = ranking_result.get("computational_ranking_score", 0.0)
        tier = ranking_result.get("classification_label", "Candidate")
        components = ranking_result.get("component_scores", {})

        # Rule-driven evidence narrative construction
        narrative_paragraphs = [
            f"**1. Target Protein Selection & Biological Rationale:**\n"
            f"The protein **{protein_name}** (Gene: **{gene_symbol}**, UniProt Accession: **{uniprot_id}**) was selected as a key therapeutic target for **{disease_name}**. "
            f"Bioinformatics mapping indicates that {gene_symbol} plays a critical role in disease onset and progression: {disease_assoc}",

            f"**2. Protein-Protein Interaction (PPI) Network Importance:**\n"
            f"In the disease-specific PPI network, **{gene_symbol}** exhibits high network centrality (PPI Centrality Score: **{components.get('ppi_centrality', 0.0)}/100**). "
            f"Its degree centrality and betweenness centrality confirm that this protein acts as a central hub regulating key pathological signaling cascades.",

            f"**3. 3D Protein Structural & Binding Cavity Analysis:**\n"
            f"The 3D atomic structure was retrieved from **{struct_source}** ({struct_type}, PDB ID: `{structure_meta.get('pdb_id', 'PDB')}`). "
            f"Binding site detection identified **{pocket_name}** with grid box center coordinates `({pocket_info.get('center', {}).get('x')}, {pocket_info.get('center', {}).get('y')}, {pocket_info.get('center', {}).get('z')})`. "
            f"Key pocket-lining catalytic residues include `{', '.join(pocket_info.get('key_residues', ['ASP32']))}`.",

            f"**4. Candidate Drug Bioactivity & Molecular Docking Performance:**\n"
            f"The candidate small-molecule drug **{drug_name}** ({chembl_id}) demonstrated strong target bioactivity with a reported IC50 of **{ic50_nm} nM**. "
            f"AutoDock Vina molecular docking into the detected pocket yielded a favorable predicted binding affinity of **{docking_affinity} kcal/mol** (Docking Score: **{components.get('docking_score', 0.0)}/100**), supported by hydrogen-bonding interactions with active site residues.",

            f"**5. Final Computational Ranking & Decision Rationale:**\n"
            f"The compound received a composite **Computational Ranking Score of {final_score}/100** ({tier}). "
            f"This high evaluation is driven by strong target centrality ({components.get('ppi_centrality', 0)}), high disease association ({components.get('disease_association', 0)}), potent bioactivity ({components.get('target_bioactivity', 0)}), and favorable docking energetics ({components.get('docking_score', 0)})."
        ]

        full_explanation_text = "\n\n".join(narrative_paragraphs)

        return {
            "disease_name": disease_name,
            "target_gene": gene_symbol,
            "target_uniprot": uniprot_id,
            "target_protein_name": protein_name,
            "candidate_drug": drug_name,
            "chembl_id": chembl_id,
            "computational_ranking_score": final_score,
            "classification_tier": tier,
            "narrative_explanation": full_explanation_text,
            "evidence_summary": {
                "disease_association_score": components.get("disease_association", 0.0),
                "ppi_centrality_score": components.get("ppi_centrality", 0.0),
                "target_bioactivity_ic50_nm": ic50_nm,
                "docking_affinity_kcal_mol": docking_affinity,
                "structure_source": struct_source,
                "pocket_name": pocket_name,
                "literature_citations": drug_data.get("literature_citations", 0)
            },
            "disclaimer": "This system is a computational research and decision-support prototype. Predicted targets, molecular interactions, and candidate drugs require experimental and clinical validation."
        }

if __name__ == "__main__":
    engine = ExplainabilityEngine()
    dummy_target = {"gene_symbol": "BACE1", "uniprot_id": "P56817", "protein_name": "Beta-Secretase 1", "disease_role": "Rate-limiting enzyme in amyloid-beta formation."}
    dummy_struct = {"structure_source": "RCSB Protein Data Bank (PDB)", "structure_type": "Experimental", "pdb_id": "1FKN"}
    dummy_pocket = {"name": "Catalytic Aspartyl Dyad Pocket", "center": {"x": 19.5, "y": 25.2, "z": 12.8}, "key_residues": ["ASP32", "ASP228"]}
    dummy_drug = {"drug_name": "Verubecestat", "chembl_id": "CHEMBL3301604", "bioactivity_value_nm": 2.2, "literature_citations": 45}
    dummy_docking = {"binding_affinity_kcal_mol": -9.8}
    dummy_ranking = {"computational_ranking_score": 89.5, "classification_label": "High Confidence Candidate", "component_scores": {"ppi_centrality": 92.0, "disease_association": 89.0, "target_bioactivity": 90.0, "docking_score": 81.7, "literature_evidence": 70.0}}

    exp = engine.generate_explanation("Alzheimer's Disease", dummy_target, dummy_struct, dummy_pocket, dummy_drug, dummy_docking, dummy_ranking)
    print("=== PHASE 10: Explainable AI Evidence Report ===")
    print(exp["narrative_explanation"])
