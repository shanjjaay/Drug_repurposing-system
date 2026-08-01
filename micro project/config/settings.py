import os
from pathlib import Path

# Base Paths
BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
CACHE_DIR = DATA_DIR / "cache"
DISEASE_DIR = DATA_DIR / "diseases"
PROTEIN_DIR = DATA_DIR / "proteins"
STRUCTURE_DIR = DATA_DIR / "structures"
DRUG_DIR = DATA_DIR / "drugs"
DOCKING_DIR = DATA_DIR / "docking"
REFERENCE_DIR = DATA_DIR / "references"

# Ensure all directories exist
for path in [DATA_DIR, CACHE_DIR, DISEASE_DIR, PROTEIN_DIR, STRUCTURE_DIR, DRUG_DIR, DOCKING_DIR, REFERENCE_DIR]:
    path.mkdir(parents=True, exist_ok=True)

# Database
DB_PATH = DATA_DIR / "database.db"

# API Endpoints
OPEN_TARGETS_API = "https://api.platform.opentargets.org/api/v4/graphql"
STRING_API_URL = "https://string-db.org/api"
UNIPROT_API_URL = "https://rest.uniprot.org/uniprotkb"
RCSB_PDB_API_URL = "https://data.rcsb.org/rest/v1/core/entry"
ALPHAFOLD_API_URL = "https://alphafold.ebi.ac.uk/api/prediction"
CHEMBL_API_URL = "https://www.ebi.ac.uk/chembl/api/data"

# Default Supported Diseases
SUPPORTED_DISEASES = {
    "alzheimers": "Alzheimer's Disease",
    "diabetes_t2": "Type 2 Diabetes",
    "breast_cancer": "Breast Cancer"
}

# Multi-Factor Drug Scoring Weights (Total = 1.0 / 100%)
SCORING_WEIGHTS = {
    "ppi_centrality": 0.25,
    "disease_association": 0.20,
    "target_bioactivity": 0.25,
    "docking_score": 0.20,
    "literature_evidence": 0.10
}
