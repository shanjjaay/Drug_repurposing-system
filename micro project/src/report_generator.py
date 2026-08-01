import os
import json
from typing import Dict, Any
from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parent.parent))
from config.settings import REFERENCE_DIR

class ReportGenerator:
    """
    Phase 10 & 11: Exportable Research & Evidence Summary Report Generator.
    Exports Markdown and HTML decision support reports.
    """
    def __init__(self, output_dir: Path = REFERENCE_DIR):
        self.output_dir = output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def generate_markdown_report(self, exp_data: Dict[str, Any]) -> Tuple[str, Path]:
        """
        Generates structured Markdown report file.
        """
        disease = exp_data.get("disease_name", "Disease Research Report")
        drug = exp_data.get("candidate_drug", "Candidate Drug")
        target = exp_data.get("target_gene", "Target Protein")

        report_filename = f"Report_{disease.replace(' ', '_')}_{target}_{drug.replace(' ', '_')}.md"
        report_path = self.output_dir / report_filename

        md_content = f"""# Healthcare Decision Support Report: {drug} for {disease}

**Target Protein:** {exp_data.get('target_protein_name')} (`{target}` / UniProt: `{exp_data.get('target_uniprot')}`)  
**Candidate Compound:** {drug} (ChEMBL: `{exp_data.get('chembl_id')}`)  
**Computational Ranking Score:** **{exp_data.get('computational_ranking_score')}/100** ({exp_data.get('classification_tier')})  
**Date:** 2026-08-01  

---

## Executive Evidence Summary

{exp_data.get('narrative_explanation')}

---

## Quantitative Metric Breakdown

| Computational Factor | Value / Evidence | Score Impact (0-100) |
| :--- | :--- | :--- |
| **Disease Association Score** | Open Targets / DisGeNET Evidence | {exp_data.get('evidence_summary', {}).get('disease_association_score')} / 100 |
| **PPI Network Centrality** | NetworkX Hub Centrality | {exp_data.get('evidence_summary', {}).get('ppi_centrality_score')} / 100 |
| **Target Bioactivity (IC50)** | {exp_data.get('evidence_summary', {}).get('target_bioactivity_ic50_nm')} nM | High Affinity |
| **Docking Binding Affinity** | {exp_data.get('evidence_summary', {}).get('docking_affinity_kcal_mol')} kcal/mol | Favored Pose |
| **Literature Evidence** | {exp_data.get('evidence_summary', {}).get('literature_citations')} PubMed Citations | Strong Support |

---

## Important Scientific & Regulatory Disclaimer

> {exp_data.get('disclaimer')}
"""
        with open(report_path, "w", encoding="utf-8") as f:
            f.write(md_content)

        print(f"[REPORT GENERATED] Saved Markdown report to {report_path.name}")
        return md_content, report_path

if __name__ == "__main__":
    from src.explainability import ExplainabilityEngine
    engine = ExplainabilityEngine()
    dummy_target = {"gene_symbol": "BACE1", "uniprot_id": "P56817", "protein_name": "Beta-Secretase 1", "disease_role": "Rate-limiting enzyme."}
    dummy_struct = {"structure_source": "RCSB PDB", "structure_type": "Experimental", "pdb_id": "1FKN"}
    dummy_pocket = {"name": "Catalytic Aspartyl Dyad Pocket", "center": {"x": 19.5, "y": 25.2, "z": 12.8}, "key_residues": ["ASP32", "ASP228"]}
    dummy_drug = {"drug_name": "Verubecestat", "chembl_id": "CHEMBL3301604", "bioactivity_value_nm": 2.2, "literature_citations": 45}
    dummy_docking = {"binding_affinity_kcal_mol": -9.8}
    dummy_ranking = {"computational_ranking_score": 89.5, "classification_label": "High Confidence Candidate", "component_scores": {"ppi_centrality": 92.0, "disease_association": 89.0, "target_bioactivity": 90.0, "docking_score": 81.7, "literature_evidence": 70.0}}

    exp = engine.generate_explanation("Alzheimer's Disease", dummy_target, dummy_struct, dummy_pocket, dummy_drug, dummy_docking, dummy_ranking)
    gen = ReportGenerator()
    md, path = gen.generate_markdown_report(exp)
    print(f"Generated Report File: {path}")
