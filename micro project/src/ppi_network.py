import json
import requests
import networkx as nx
import pandas as pd
from typing import Dict, List, Any, Tuple
from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parent.parent))
from config.settings import CACHE_DIR, STRING_API_URL

class PPINetworkBuilder:
    """
    Phase 2: Protein-Protein Interaction (PPI) Network Construction.
    Fetches interactions from STRING DB API or local cache,
    and constructs a weighted NetworkX undirected graph.
    """
    def __init__(self, cache_dir: Path = CACHE_DIR):
        self.cache_dir = cache_dir
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def build_ppi_graph(self, disease_key: str = "alzheimers", protein_symbols: List[str] = None) -> Tuple[nx.Graph, Dict[str, Any]]:
        """
        Builds a NetworkX weighted undirected graph for given protein symbols or disease key.
        Nodes = Protein symbols / UniProt IDs.
        Edges = Interaction confidence scores (0.0 - 1.0).
        """
        ppi_data = self.fetch_ppi_data(disease_key, protein_symbols)
        G = nx.Graph()

        # Add Nodes
        for node in ppi_data.get("nodes", []):
            symbol = node["symbol"]
            G.add_node(symbol, uniprot_id=node.get("id"), name=node.get("name"))

        # Add Edges
        for edge in ppi_data.get("edges", []):
            src = edge["source"]
            tgt = edge["target"]
            weight = float(edge.get("combined_score", 0.5))
            if G.has_node(src) and G.has_node(tgt):
                G.add_edge(src, tgt, weight=weight, combined_score=weight)

        print(f"[PPI GRAPH] Built graph with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges.")
        return G, ppi_data

    def fetch_ppi_data(self, disease_key: str = "alzheimers", protein_symbols: List[str] = None) -> Dict[str, Any]:
        """
        Retrieves PPI network data from local cache or STRING REST API.
        """
        disease_clean = disease_key.lower().replace(" ", "_").replace("'", "")
        cache_path = self.cache_dir / f"{disease_clean}_ppi.json"

        if cache_path.exists():
            with open(cache_path, "r", encoding="utf-8") as f:
                print(f"[CACHE HIT] Loaded PPI network for {disease_key} from {cache_path.name}")
                return json.load(f)

        # Fallback to STRING DB REST API call if live internet is active
        if protein_symbols:
            try:
                string_data = self._query_string_api(protein_symbols)
                if string_data and "edges" in string_data and len(string_data["edges"]) > 0:
                    with open(cache_path, "w", encoding="utf-8") as f:
                        json.dump(string_data, f, indent=2)
                    return string_data
            except Exception as e:
                print(f"[STRING API WARNING] API query failed ({e}). Returning fallback template.")

        # Fallback empty structure
        return {"nodes": [], "edges": []}

    def _query_string_api(self, protein_list: List[str]) -> Dict[str, Any]:
        """
        Queries STRING DB network REST API for interaction edges.
        """
        url = f"{STRING_API_URL}/json/network"
        params = {
            "identifiers": "%0d".join(protein_list),
            "species": 9606,  # Human
            "required_score": 700  # High confidence > 0.7
        }
        resp = requests.post(url, data=params, timeout=10)
        if resp.status_code == 200:
            interactions = resp.json()
            nodes_set = set()
            edges = []
            for item in interactions:
                src = item.get("preferredName_A")
                tgt = item.get("preferredName_B")
                score = float(item.get("score", 0.7))
                nodes_set.add(src)
                nodes_set.add(tgt)
                edges.append({"source": src, "target": tgt, "combined_score": score})

            nodes = [{"symbol": symbol, "id": symbol, "name": symbol} for symbol in nodes_set]
            return {"disease": "Live API Query", "source": "STRING DB REST API", "nodes": nodes, "edges": edges}
        return {"nodes": [], "edges": []}

if __name__ == "__main__":
    builder = PPINetworkBuilder()
    G, data = builder.build_ppi_graph("alzheimers")
    print("\nGraph Summary:")
    print("Nodes:", list(G.nodes()))
    print("Edges:", [(u, v, d['weight']) for u, v, d in G.edges(data=True)])
