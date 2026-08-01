import unittest
import sys
from pathlib import Path

# Add project root to path
sys.path.append(str(Path(__file__).resolve().parent.parent))

from src.disease_analysis import DiseaseAnalyzer

class TestDiseaseAnalysis(unittest.TestCase):
    def setUp(self):
        self.analyzer = DiseaseAnalyzer()

    def test_fetch_alzheimers_associations(self):
        data = self.analyzer.fetch_disease_associations("alzheimers")
        self.assertIsNotNone(data)
        self.assertIn("associations", data)
        self.assertEqual(data["disease_name"], "Alzheimer's Disease")
        self.assertGreater(len(data["associations"]), 0)

    def test_dataframe_generation(self):
        df = self.analyzer.get_association_dataframe("alzheimers")
        self.assertFalse(df.empty)
        self.assertIn("gene_symbol", df.columns)
        self.assertIn("uniprot_id", df.columns)
        self.assertIn("association_score", df.columns)

    def test_target_uniprot_ids(self):
        ids = self.analyzer.get_target_uniprot_ids("alzheimers", min_score=0.85)
        self.assertIsInstance(ids, list)
        self.assertIn("P05067", ids)  # APP
        self.assertIn("P10636", ids)  # MAPT
        self.assertIn("P56817", ids)  # BACE1

if __name__ == "__main__":
    unittest.main()
