"""
Point d'entrée de biosécurité (T10) : génère ``data/processed/toxin_flags.json``.

Usage :
    python -m peplib.check_toxins

La logique de détection vit dans :mod:`peplib.biosafety` et est couverte par
``tests/test_biosafety.py``.
"""

from peplib.biosafety import run

RAW_FILE = "data/raw/uniprot_reviewed_20261017.tsv"
OUT_FLAGS = "data/processed/toxin_flags.json"


def main(raw_file=RAW_FILE, out_file=OUT_FLAGS):
    toxin_count = run(raw_file, out_file)
    print(f"Nombre de toxines détectées et flaggées : {toxin_count}")
    print(f"Rapport de biosécurité enregistré sous {out_file}")


if __name__ == "__main__":
    main()
