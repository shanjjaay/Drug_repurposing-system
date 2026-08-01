import networkx as nx
import pandas as pd
import numpy as np
from typing import Dict, List, Any, Tuple
from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parent.parent))
from src.disease_analysis import DiseaseAnalyzer
from src.ppi_network import PPINetworkBuilder

class TargetRanker:
    """
    Phase 2: Target Protein Identification & Network Centrality Ranking.
    Calculates Degree, Betweenness, Eigenvector Centrality, and PageRank,
    normalizes values, and ranks potential therapeutic targets.
    """
    def __init__(self, w_degree: float = 0.35, w_betweenness: float = 0.35, w_eigenvector: float = 0.30):
        self.w_degree = w_degree
        self.w_betweenness = w_betweenness
        self.w_eigenvector = w_eigenvector

    def calculate_centralities(self, G: nx.Graph) -> pd.DataFrame:
        """
        Calculates Network Centralities for all protein nodes in graph G.
        """
        if G.number_of_nodes() == 0:
            return pd.DataFrame()

        # Centrality metrics
        deg_centrality = nx.degree_centrality(G)
        try:
            btw_centrality = nx.betweenness_centrality(G, weight="weight")
        except Exception:
            btw_centrality = {node: 0.0 for node in G.nodes()}

        try:
            eig_centrality = nx.eigenvector_centrality(G, weight="weight", max_iter=1000)
        except Exception:
            eig_centrality = {node: 1.0 / G.number_of_nodes() for node in G.nodes()}

        try:
            pagerank = nx.pagerank(G, weight="weight")
        except Exception:
            pagerank = {node: 1.0 / G.number_of_nodes() for node in G.nodes()}

        nodes_data = []
        for node in G.nodes(data=True):
            symbol = node[0]
            attr = node[1]
            nodes_data.append({
                "gene_symbol": symbol,
                "uniprot_id": attr.get("uniprot_id", symbol),
                "protein_name": attr.get("name", symbol),
                "degree_centrality": deg_centrality.get(symbol, 0.0),
                "betweenness_centrality": btw_centrality.get(symbol, 0.0),
                "eigenvector_centrality": eig_centrality.get(symbol, 0.0),
                "pagerank": pagerank.get(symbol, 0.0)
            })

        df = pd.DataFrame(nodes_data)

        # Min-Max Normalization (0.0 to 1.0)
        for col in ["degree_centrality", "betweenness_centrality", "eigenvector_centrality", "pagerank"]:
            min_val = df[col].min()
            max_val = df[col].max()
            if max_val > min_val:
                df[f"{col}_norm"] = (df[col] - min_val) / (max_val - min_val)
            else:
                df[f"{col}_norm"] = 1.0

        # Composite Centrality Score
        df["ppi_centrality_score"] = (
            self.w_degree * df["degree_centrality_norm"] +
            self.w_betweenness * df["betweenness_centrality_norm"] +
            self.w_eigenvector * df["eigenvector_centrality_norm"]
        )

        df = df.sort_values(by="ppi_centrality_score", ascending=False).reset_index(drop=True)
        return df

    def rank_target_proteins(self, disease_key: str = "alzheimers", top_n: int = 5) -> Tuple[pd.DataFrame, nx.Graph]:
        """
        Combines Disease Association Evidence + PPI Network Centrality to select top therapeutic targets.
        """
        # Step 1: Load Disease Associations
        disease_analyzer = DiseaseAnalyzer()
        disease_df = disease_analyzer.get_association_dataframe(disease_key)

        # Step 2: Build PPI Graph
        ppi_builder = PPINetworkBuilder()
        G, ppi_data = ppi_builder.build_ppi_graph(disease_key)

        # Step 3: Compute Centralities
        centrality_df = self.calculate_centralities(G)

        # Step 4: Merge Disease Association Score with Centrality
        if not disease_df.empty and not centrality_df.empty:
            merged_df = pd.merge(disease_df, centrality_df, on="gene_symbol", how="inner", suffixes=("_disease", "_ppi"))
            if "uniprot_id_disease" in merged_df.columns:
                merged_df["uniprot_id"] = merged_df["uniprot_id_disease"]

            # Combined Target Selection Score
            merged_df["target_importance_score"] = (
                0.60 * merged_df["ppi_centrality_score"] +
                0.40 * merged_df["association_score"]
            )
            merged_df = merged_df.sort_values(by="target_importance_score", ascending=False).reset_index(drop=True)
            top_targets = merged_df.head(top_n)
            return top_targets, G

        return centrality_df.head(top_n), G

if __name__ == "__main__":
    ranker = TargetRanker()
    top_df, G = ranker.rank_target_proteins("alzheimers", top_n=5)
    print("=== PHASE 2: Target Protein Network Ranking ===")
    print(top_df[["gene_symbol", "uniprot_id", "gene_name", "degree_centrality", "betweenness_centrality", "ppi_centrality_score", "target_importance_score"]])
