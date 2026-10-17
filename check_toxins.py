"""
Point d'entrée de biosécurité (T10) : génère ``data/processed/toxin_flags.json``.

La logique de détection vit dans :mod:`peplib.biosafety` et est couverte par
``tests/test_biosafety.py``. Ce script ne fait qu'appeler le module.
"""

from peplib.biosafety import run

RAW_FILE = "data/raw/uniprot_reviewed_20261017.tsv"
OUT_FLAGS = "data/processed/toxin_flags.json"


if __name__ == "__main__":
    toxin_count = run(RAW_FILE, OUT_FLAGS)
    print(f"Nombre de toxines détectées et flaggées : {toxin_count}")
    print(f"Rapport de biosécurité enregistré sous {OUT_FLAGS}")
