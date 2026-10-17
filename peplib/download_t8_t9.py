"""
Point d'entrée T8/T9 : télécharge les corpus AMP et bactériocines.

Usage :
    python -m peplib.download_t8_t9

La logique vit dans :mod:`peplib.download` et est couverte par
``tests/test_download.py``.
"""

from peplib.download import (
    compute_sha256,
    count_lines,
    download_t8,
    download_t9,
    raw_path,
)


def _report(label, filepath, total_sequences):
    lines = count_lines(filepath)
    sha = compute_sha256(filepath)
    print(f"{label} terminé ! Fichier : {filepath}")
    print(f"   Lignes brutes : {lines} ({total_sequences} séquences)")
    print(f"   SHA-256 : {sha}\n")


def main():
    print("1. Téléchargement de T8 (peptides antimicrobiens)...")
    amp_path = raw_path("amp_database")
    _report("T8", amp_path, download_t8(amp_path))

    print("2. Téléchargement de T9 (bactériocines)...")
    bact_path = raw_path("bactibase_bacteriocins")
    _report("T9", bact_path, download_t9(bact_path))


if __name__ == "__main__":
    main()
