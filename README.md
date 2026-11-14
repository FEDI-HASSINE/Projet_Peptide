# Projet Peptide

`Projet_Peptide` is a reproducible Python project for collecting, cleaning,
analysing, and preparing peptide sequences for downstream machine-learning
experiments.

The project focuses on the data and validation steps that must be completed
before training a generative model: source tracking, sequence filtering,
duplicate detection, similarity clustering, physicochemical descriptors,
position-weight matrices, data-leakage analysis, and model-ready splits.

## Project objectives

The pipeline is designed to produce a clean and auditable peptide corpus from
raw biological data sources. It provides:

- reproducible downloads and SHA-256 provenance records;
- strict filtering to canonical amino-acid sequences of length 10 to 50;
- exact and near-duplicate analysis using k-mer indexing and Levenshtein
  distance;
- sequence clustering and similarity analysis;
- physicochemical descriptors such as charge, pI, GRAVY, hydrophobic moment,
  Boman index, and molecular mass;
- composition statistics, entropy, KL divergence, and Jensen-Shannon
  divergence;
- position-weight matrices and information profiles aligned from the N- and
  C-termini;
- toxin flagging and biosafety checks;
- family-aware train, validation, and test splits intended to limit data
  leakage;
- configuration templates for LSTM and Transformer experiments.

## Repository structure

```text
.
├── configs/                  # LSTM and Transformer experiment configurations
├── data/
│   ├── raw/                  # Original downloaded files; never edit in place
│   └── processed/            # Filtered and derived datasets
├── Docs/                    # Methodology, reports, and project documentation
├── models/                  # Model artefacts
├── notebooks/               # Exploratory notebooks
├── peplib/                  # Reusable peptide-processing library
├── results/                 # Tables and figures produced by the analyses
├── tests/                   # Automated tests
├── PROVENANCE.md            # Sources, licences, queries, and file hashes
├── requirements.txt         # Pinned Python dependencies
└── splits_manifest.json     # Split-generation manifest
```

## Installation

The project requires Python 3.10 or later.

```bash
python -m venv .venv
```

Activate the virtual environment:

```bash
# Linux or macOS
source .venv/bin/activate

# Windows PowerShell
.venv\Scripts\Activate.ps1
```

Install the pinned dependencies:

```bash
python -m pip install -r requirements.txt
```

## Running the tests

Run the complete test suite from the repository root:

```bash
pytest
```

The tests cover filtering, statistics, descriptors, alignment, clustering,
deduplication, tokenization, data loading, biosafety checks, and generated
artefacts.

## Main workflows

### Download source data

The download modules retrieve the documented UniProt-based datasets and record
their provenance:

```bash
python -m peplib.download_uniprot
python -m peplib.download_t8_t9
```

Additional source preparation modules are available for APD3 and the peptide
challenge dataset:

```bash
python -m peplib.download_apd3
python -m peplib.download_challenge
```

### Filter and prepare the corpus

The filtering pipeline validates the amino-acid alphabet, removes missing
values, enforces the 10–50 residue length window, and writes processed data:

```bash
python -m peplib.filters
python -m peplib.build_corpus_v02
```

### Run quality and safety checks

```bash
python -m peplib.check_toxins
python -m peplib.check_toxins_v02
python -m peplib.check_near_duplicates
python -m peplib.dedup
python -m peplib.cluster
```

Toxin flags are an essential review step. Flagged sequences must be excluded
from any generative corpus unless they have been explicitly reviewed and
approved under the project’s safety policy.

### Generate analyses and figures

```bash
python -m peplib.stats
python -m peplib.pwm
python -m peplib.figures
python -m peplib.check_overlaps
```

The generated tables and figures are written to `results/` and documented
outputs under `data/processed/`.

## Reproducibility

The project uses a fixed seed of `42` for Python, NumPy, and PyTorch workflows.
Raw files are preserved separately from processed outputs, and their source,
query, licence, date, line count, and SHA-256 hash are recorded in
[`PROVENANCE.md`](PROVENANCE.md).

Do not modify files in `data/raw/`. Regenerate processed files from the raw
inputs using the library workflows instead.

## Data provenance and responsible use

The repository contains provenance information for UniProt-based datasets and
documents the applicable usage conditions. Some biological databases contain
toxin-related sequences. The project therefore flags toxin candidates and
requires them to be excluded from generated or redistributed peptide sets.

This repository is intended for research and educational use. It does not
provide medical advice, therapeutic recommendations, or validated biological
activity claims.

## Model configuration

The `configs/` directory contains validated baseline configurations for:

- a two-layer LSTM language model in `configs/lstm_base.yaml`;
- a four-layer Transformer model in `configs/transformer_base.yaml`.

These files define model parameters, training settings, data paths, and
sampling grids. They are configuration templates for downstream experiments;
the core repository currently concentrates on corpus preparation and analysis.

## Documentation

- [`PROVENANCE.md`](PROVENANCE.md): data sources, licences, queries, hashes,
  and biosafety notes.
- [`Docs/Semaine1_peptides_peplib_4etudiants.md`](Docs/Semaine1_peptides_peplib_4etudiants.md):
  project specification and workflow.
- [`Docs/rapport_T25.md`](Docs/rapport_T25.md): analysis of near-duplicate
  sequences at 90% identity.
- [`notebooks/NB01_peptides.ipynb`](notebooks/NB01_peptides.ipynb):
  exploratory peptide analysis.

## Current status

The repository contains the tested `peplib` analysis library, corpus
preparation workflows, provenance records, quality checks, and baseline model
configurations. Training and evaluation of the final generative models are
planned as a subsequent stage after the corpus and leakage controls are
finalised.

## Licence and attribution

Before redistributing datasets or derived sequences, consult the source
licences and the project provenance records. Cite the original data providers
and this repository when using the code or analyses.
