# pce26 — Protein Cage Engineering 2026

Research repository for computational protein cage design and docking studies.

---

## Repository rules

These rules keep the repo clean, reproducible, and easy to navigate for everyone.

### What belongs here

| Directory | Allowed content |
|---|---|
| `data/` | Numeric/tabular files only: `.csv`, `.tsv`, `.json`, `.pkl`, `.h5`, `.npy`, `.npz`, `.parquet` |
| `structures/` | Protein structure files: `.pdb`, `.pdb.gz`, `.cif`, `.mmcif`, and compressed archives (`.txz`, `.tar.gz`) |
| `notebooks/` | Jupyter notebooks (`.ipynb`) for analysis and figures. Name them `YYYY-MM-DD.<topic>.ipynb` |
| `tests/` | Test inputs and expected outputs for computational tools (e.g. rpxdock runs) |
| `scripts/` | Standalone Python/bash scripts for processing or running jobs |

### What does NOT belong here

- **No large binary files** (> ~50 MB) — use a shared storage location or data archive instead
- **No raw MD trajectories or large simulation outputs** — keep those on the cluster, link to them in the notebook
- **No environment files** — do not commit `conda-lock.yml`, `.env`, or credential files
- **No output that can be regenerated** — if a file can be produced by running a script, it should not be tracked (add it to `.gitignore`)

---

## Directory layout

```
pce26/
├── data/           # Numeric datasets (CSV, NPY, HDF5, …)
├── notebooks/      # Jupyter notebooks — one per analysis session
├── structures/     # PDB and CIF protein structure files
├── tests/          # Test cases for docking / processing pipelines
└── scripts/        # Utility and job-submission scripts
```

---

## Getting started

### 1. Clone the repo

```bash
git clone <repo-url>
cd pce26
```

### 2. Set up your environment

Create a conda environment with the required packages:

```bash
conda create -n pce26 python=3.11
conda activate pce26
pip install numpy pandas matplotlib jupyterlab biopython
```

If a specific tool (e.g. **rpxdock**) is needed, follow its own installation instructions and make sure it is accessible in your `$PATH` before running notebooks.

### 3. Start a notebook

```bash
jupyter lab notebooks/
```

Name every new notebook with today's date and a short topic:

```
notebooks/YYYY-MM-DD.<topic>.ipynb
```

### 4. Working with structures

PDB files live in `structures/`. To load one in Python:

```python
from Bio.PDB import PDBParser
parser = PDBParser(QUIET=True)
structure = parser.get_structure("7MRY", "structures/7MRY.pdb")
```

### 5. Committing changes

Only commit files that belong to the categories listed above. Before committing, check:

```bash
git status          # review what changed
git diff --stat     # make sure no large or unintended files are staged
git add <files>
git commit -m "short description of what changed"
```

---

## Naming conventions

- **Structures** — use the PDB ID and chain: `7MRY.pdb`, `7MRY.A.pdb`
- **Notebooks** — `YYYY-MM-DD.<tool-or-topic>.ipynb` (e.g. `2026-03-06.rpxdock.ipynb`)
- **Data files** — descriptive snake_case: `rpxdock_scores_7mry.csv`

---

## Questions?

Open an issue in this repo or reach out to Chloe
