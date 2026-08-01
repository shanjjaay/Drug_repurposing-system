import unittest
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))
from src.target_ranking import TargetRanker
from src.ppi_network import PPINetworkBuilder

class TestTargetRanking(unittest.TestCase):
    def setUp(self):
        self.ranker = TargetRanker()
        self.ppi_builder = PPINetworkBuilder()

    def test_ppi_graph_construction(self):
        G, data = self.ppi_builder.build_ppi_graph("alzheimers")
        self.assertGreater(G.number_of_nodes(), 0)
        self.assertGreater(G.number_of_edges(), 0)

    def test_centrality_calculation(self):
        G, data = self.ppi_builder.build_ppi_graph("alzheimers")
        df = self.ranker.calculate_centralities(G)
        self.assertFalse(df.empty)
        self.assertIn("degree_centrality", df.columns)
        self.assertIn("betweenness_centrality", df.columns)
        self.assertIn("ppi_centrality_score", df.columns)

    def test_top_targets_selection(self):
        top_df, G = self.ranker.rank_target_proteins("alzheimers", top_n=3)
        self.assertEqual(len(top_df), 3)
        self.assertIn("target_importance_score", top_df.columns)
        top_symbols = top_df["gene_symbol"].tolist()
        self.assertTrue("APP" in top_symbols or "BACE1" in top_symbols or "MAPT" in top_symbols)

if __name__ == "__main__":
    unittest.main()
