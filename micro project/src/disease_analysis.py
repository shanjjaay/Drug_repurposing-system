import json
import os
import requests
import pandas as pd
from typing import Dict, List, Any, Optional
from pathlib import Path
import sys

# Ensure project root is in sys.path
sys.path.append(str(Path(__file__).resolve().parent.parent))
from config.settings import DISEASE_DIR, CACHE_DIR, OPEN_TARGETS_API, SUPPORTED_DISEASES

class DiseaseAnalyzer:
    """
    Phase 1: Disease-to-Gene-to-Protein Association Identification.
    Retrieves and maps disease-associated genes, proteins (UniProt IDs),
    association scores, and biological roles.
    """
    def __init__(self, cache_dir: Path = CACHE_DIR, disease_dir: Path = DISEASE_DIR):
        self.cache_dir = cache_dir
        self.disease_dir = disease_dir
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def fetch_disease_associations(self, disease_key: str = "alzheimers", use_live_api: bool = False) -> Dict[str, Any]:
        """
        Fetch disease associations for the given disease key.
        Checks local cache first, fallback to local curated JSON or Open Targets API.
        """
        disease_key_clean = disease_key.lower().replace(" ", "_").replace("'", "")
        cache_path = self.cache_dir / f"{disease_key_clean}_associations.json"

        # Step 1: Check cache file
        if cache_path.exists():
            with open(cache_path, "r", encoding="utf-8") as f:
                data = json.load(f)
                print(f"[CACHE HIT] Loaded {disease_key} disease associations from {cache_path.name}")
                return data

        # Step 2: Check local curated database
        local_db_path = self.disease_dir / f"{disease_key_clean}.json"
        if local_db_path.exists():
            with open(local_db_path, "r", encoding="utf-8") as f:
                data = json.load(f)
                # Cache it for consistency
                with open(cache_path, "w", encoding="utf-8") as cf:
                    json.dump(data, cf, indent=2)
                print(f"[LOCAL DB] Loaded {disease_key} associations from {local_db_path.name}")
                return data

        # Step 3: Optional live API query (Open Targets GraphQL API)
        if use_live_api:
            try:
                live_data = self._fetch_from_open_targets(disease_key)
                if live_data and "associations" in live_data and len(live_data["associations"]) > 0:
                    with open(cache_path, "w", encoding="utf-8") as cf:
                        json.dump(live_data, cf, indent=2)
                    print(f"[LIVE API] Successfully fetched from Open Targets for {disease_key}")
                    return live_data
            except Exception as e:
                print(f"[API ERROR] Open Targets query failed: {e}. Falling back to default.")

        raise FileNotFoundError(f"No association data found for disease key '{disease_key}' in {self.disease_dir}")

    def get_association_dataframe(self, disease_key: str = "alzheimers") -> pd.DataFrame:
        """
        Converts association dataset into a clean Pandas DataFrame.
        """
        data = self.fetch_disease_associations(disease_key)
        associations = data.get("associations", [])
        df = pd.DataFrame(associations)
        if not df.empty:
            df = df.sort_values(by="association_score", ascending=False).reset_index(drop=True)
        return df

    def get_target_uniprot_ids(self, disease_key: str = "alzheimers", min_score: float = 0.5) -> List[str]:
        """
        Extracts a list of UniProt IDs for proteins exceeding the minimum association score.
        """
        df = self.get_association_dataframe(disease_key)
        if df.empty:
            return []
        filtered_df = df[df["association_score"] >= min_score]
        return filtered_df["uniprot_id"].dropna().unique().tolist()

    def _fetch_from_open_targets(self, disease_name: str) -> Dict[str, Any]:
        """
        Executes a GraphQL query against Open Targets Platform API.
        """
        query = """
        query DiseaseTargets($queryString: String!) {
          search(queryString: $queryString, entityNames: ["disease"]) {
            hits {
              id
              name
            }
          }
        }
        """
        response = requests.post(OPEN_TARGETS_API, json={"query": query, "variables": {"queryString": disease_name}}, timeout=10)
        if response.status_code == 200:
            res_json = response.json()
            hits = res_json.get("data", {}).get("search", {}).get("hits", [])
            if hits:
                disease_id = hits[0]["id"]
                disease_label = hits[0]["name"]
                return {
                    "disease_id": disease_id,
                    "disease_name": disease_label,
                    "associations": []
                }
        return {}


if __name__ == "__main__":
    analyzer = DiseaseAnalyzer()
    print("=== PHASE 1: Disease -> Gene -> Protein Analysis ===")
    data = analyzer.fetch_disease_associations("alzheimers")
    print(f"\nDisease: {data.get('disease_name')} ({data.get('disease_id')})")
    print(f"Description: {data.get('description')}\n")
    
    df = analyzer.get_association_dataframe("alzheimers")
    print("Identified Genes & Target Proteins:")
    print(df[["gene_symbol", "uniprot_id", "gene_name", "association_score"]])
    
    uniprot_ids = analyzer.get_target_uniprot_ids("alzheimers", min_score=0.8)
    print(f"\nTarget UniProt IDs (score >= 0.8): {uniprot_ids}")
