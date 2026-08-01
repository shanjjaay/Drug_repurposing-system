import pandas as pd
import numpy as np
from typing import Dict, List, Any, Tuple
from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parent.parent))
from config.settings import SCORING_WEIGHTS

class DrugRanker:
    """
    Phase 9: Custom Multi-Factor Computational Drug Ranking Engine.
    Combines:
    1. PPI Target Centrality (25%)
    2. Disease-Gene Association Score (20%)
    3. Drug-Target Bioactivity (25%)
    4. AutoDock Vina Docking Affinity (20%)
    5. Literature Citation Support (10%)

    Outputs a normalized Computational Ranking Score between 0 and 100.
    """
    def __init__(self, weights: Dict[str, float] = SCORING_WEIGHTS):
        self.weights = weights

    def compute_drug_ranking(
        self,
        target_info: Dict[str, Any],
        drug_data: Dict[str, Any],
        docking_result: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Computes multi-factor score for a single candidate drug.
        """
        # 1. PPI Centrality Score (0-1 -> 0-100)
        ppi_score = float(target_info.get("ppi_centrality_score", 0.80)) * 100.0

        # 2. Disease Association Score (0-1 -> 0-100)
        disease_assoc_score = float(target_info.get("association_score", 0.85)) * 100.0

        # 3. Drug-Target Bioactivity Score (IC50 in nM -> 0-100 scale)
        # IC50 <= 1 nM -> 100, IC50 = 10 nM -> 85, IC50 = 100 nM -> 70, IC50 >= 1000 nM -> 40
        ic50_nm = float(drug_data.get("bioactivity_value_nm", 10.0))
        if ic50_nm <= 1.0:
            target_bioactivity_score = 100.0
        elif ic50_nm <= 10.0:
            target_bioactivity_score = 90.0 - (ic50_nm - 1.0) * 1.1
        elif ic50_nm <= 100.0:
            target_bioactivity_score = 80.0 - (ic50_nm - 10.0) * 0.22
        else:
            target_bioactivity_score = max(30.0, 60.0 - (ic50_nm - 100.0) * 0.03)

        # 4. Docking Binding Affinity Score (-kcal/mol -> 0-100 scale)
        # -12.0 kcal/mol -> 100, -9.0 kcal/mol -> 80, -6.0 kcal/mol -> 50, -4.0 kcal/mol -> 20
        affinity = float(docking_result.get("binding_affinity_kcal_mol", -8.5))
        docking_score = min(100.0, max(0.0, (-affinity / 12.0) * 100.0))

        # 5. Literature Evidence Score (PubMed Citation Count -> 0-100 scale)
        citations = int(drug_data.get("literature_citations", 25))
        if citations >= 500:
            literature_score = 100.0
        elif citations >= 100:
            literature_score = 85.0
        elif citations >= 20:
            literature_score = 70.0
        else:
            literature_score = 50.0

        # Composite Mathematical Formula
        final_ranking_score = (
            self.weights["ppi_centrality"] * ppi_score +
            self.weights["disease_association"] * disease_assoc_score +
            self.weights["target_bioactivity"] * target_bioactivity_score +
            self.weights["docking_score"] * docking_score +
            self.weights["literature_evidence"] * literature_score
        )
        final_ranking_score = round(min(100.0, max(0.0, final_ranking_score)), 1)

        return {
            "drug_name": drug_data["drug_name"],
            "target_gene": target_info.get("gene_symbol", "TARGET"),
            "target_uniprot": target_info.get("uniprot_id", "P00000"),
            "computational_ranking_score": final_ranking_score,
            "component_scores": {
                "ppi_centrality": round(ppi_score, 1),
                "disease_association": round(disease_assoc_score, 1),
                "target_bioactivity": round(target_bioactivity_score, 1),
                "docking_score": round(docking_score, 1),
                "literature_evidence": round(literature_score, 1)
            },
            "docking_affinity_kcal_mol": affinity,
            "bioactivity_ic50_nm": ic50_nm,
            "scoring_formula_weights": self.weights,
            "classification_label": self._get_ranking_tier(final_ranking_score)
        }

    def rank_candidate_list(
        self,
        target_info: Dict[str, Any],
        candidate_drugs: List[Dict[str, Any]],
        docking_results: List[Dict[str, Any]]
    ) -> pd.DataFrame:
        """
        Ranks a list of candidate drugs and outputs a sorted DataFrame.
        """
        rankings = []
        docking_map = {res["drug_name"]: res for res in docking_results}

        for drug in candidate_drugs:
            d_name = drug["drug_name"]
            dock_res = docking_map.get(d_name, {
                "binding_affinity_kcal_mol": -8.0,
                "docking_score_confidence": 65.0
            })
            score_card = self.compute_drug_ranking(target_info, drug, dock_res)
            rankings.append(score_card)

        df = pd.DataFrame(rankings)
        if not df.empty:
            df = df.sort_values(by="computational_ranking_score", ascending=False).reset_index(drop=True)
            df["rank"] = df.index + 1
        return df

    def _get_ranking_tier(self, score: float) -> str:
        if score >= 85.0:
            return "High Confidence Candidate (Top Tier)"
        elif score >= 70.0:
            return "Moderate Confidence Candidate (Second Tier)"
        else:
            return "Exploratory Candidate"

if __name__ == "__main__":
    ranker = DrugRanker()
    dummy_target = {"gene_symbol": "BACE1", "uniprot_id": "P56817", "ppi_centrality_score": 0.92, "association_score": 0.89}
    dummy_drug = {"drug_name": "Verubecestat", "bioactivity_value_nm": 2.2, "literature_citations": 45}
    dummy_docking = {"binding_affinity_kcal_mol": -9.8}

    res = ranker.compute_drug_ranking(dummy_target, dummy_drug, dummy_docking)
    print("=== PHASE 9: Multi-Factor Drug Ranking Engine ===")
    print(f"Candidate: {res['drug_name']}")
    print(f"Target: {res['target_gene']} ({res['target_uniprot']})")
    print(f"Final Computational Ranking Score: {res['computational_ranking_score']} / 100")
    print(f"Tier: {res['classification_label']}")
    print("Score Components:", res['component_scores'])
