from pathlib import Path

DIRS = [
    "data/raw",
    "data/interim",
    "data/processed",
    "results/figures",
    "results/tables",
    "src",
    "notebooks",
    "reports",
]

FILES = {
    "README.md": """# Ophthalmology Transcriptomics (GEO): DEGs → Pathways → Drug Repurposing

End-to-end bioinformatics/translational project using a public GEO ophthalmology dataset.

Pipeline (step-by-step):
1) Download GEO series + build phenotype table
2) QC + normalization checks
3) Differential expression (case vs control)
4) Pathway enrichment (GO/KEGG/Reactome via Enrichr)
5) Drug repurposing (signature-based via Enrichr libraries)
6) Export tables/figures + reproducible report

> Research/education only (not clinical).
""",
    ".gitignore": """# Python
__pycache__/
*.py[cod]
.venv/
.env

# Jupyter
.ipynb_checkpoints/

# Data & outputs
data/raw/
data/interim/
data/processed/
results/
reports/*.html
reports/*.pdf

# OS
.DS_Store
Thumbs.db
""",
    "requirements.txt": """numpy
pandas
matplotlib
scipy
statsmodels
tqdm
requests
GEOparse
gseapy
openpyxl
jupyter
ipykernel
""",
    "src/config.py": """from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

DATA = ROOT / "data"
RAW = DATA / "raw"
INTERIM = DATA / "interim"
PROCESSED = DATA / "processed"

RESULTS = ROOT / "results"
FIGURES = RESULTS / "figures"
TABLES = RESULTS / "tables"

REPORTS = ROOT / "reports"
""",
    "src/utils.py": """import subprocess
from pathlib import Path

def run(cmd: list[str], cwd: Path | None = None) -> None:
    print("\\n>>", " ".join(cmd))
    subprocess.run(cmd, cwd=cwd, check=True)
""",
    "src/01_download_geo.py": """\"\"\"Step 1: Download GEO dataset + create phenotype table.\"\"\"""",
    "src/02_qc.py": """\"\"\"Step 2: QC checks + basic normalization diagnostics.\"\"\"""",
    "src/03_differential_expression.py": """\"\"\"Step 3: Differential expression + FDR correction.\"\"\"""",
    "src/04_enrichment.py": """\"\"\"Step 4: Pathway enrichment (Enrichr via gseapy).\"\"\"""",
    "src/05_drug_repurposing.py": """\"\"\"Step 5: Drug repurposing using gene signatures.\"\"\"""",
    "notebooks/01_quickstart.ipynb": "",
}

def main() -> None:
    root = Path.cwd()
    for d in DIRS:
        (root / d).mkdir(parents=True, exist_ok=True)

    for rel, content in FILES.items():
        path = root / rel
        if path.exists():
            continue
        if rel.endswith(".ipynb"):
            path.write_text(content, encoding="utf-8")
        else:
            path.write_text(content, encoding="utf-8")

    print("✅ Project skeleton created.")

if __name__ == "__main__":
    main()
