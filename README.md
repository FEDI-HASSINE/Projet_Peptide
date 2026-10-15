# Projet Peptide — Semaine 1 (peplib v0.1)

## Environnement
- Python >= 3.10 (venv)
- `pip install -r requirements.txt`
- Seed : 42 (random, NumPy, PyTorch)

## Arborescence
- `data/raw/` : fichiers bruts (jamais modifiés, empreinte SHA-256 dans PROVENANCE.md)
- `data/processed/` : CSV filtrés (sequence, source, id)
- `peplib/` : bibliothèque Python (constantes, filtres, stats, PWM, descripteurs, dédoublonnage)
- `tests/` : un fichier de tests par module
- `notebooks/NB01_peptides.ipynb` : appelle peplib uniquement

## Commandes
```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
pip install -r requirements.txt
pytest
```
