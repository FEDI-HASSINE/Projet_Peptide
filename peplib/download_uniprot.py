"""
Point d'entrée T7 : télécharge le corpus général depuis UniProtKB.

Usage :
    python -m peplib.download_uniprot

La logique vit dans :mod:`peplib.download` et est couverte par
``tests/test_download.py``.
"""

from peplib.download import compute_sha256, count_lines, download_t7, raw_path


def main():
    filepath = raw_path("uniprot_reviewed")

    print("Début du téléchargement depuis UniProt (T7)...")
    total_sequences = download_t7(filepath)

    raw_lines = count_lines(filepath)
    sha256 = compute_sha256(filepath)

    print("\n" + "=" * 50)
    print("TÉLÉCHARGEMENT UNIPROT TERMINÉ (T7) !")
    print(f"Fichier sauvegardé : {filepath}")
    print(
        f"Nombre de lignes brutes : {raw_lines} "
        f"(soit {total_sequences} séquences + 1 ligne d'en-tête)"
    )
    print(f"Empreinte SHA-256 : {sha256}")
    print("=" * 50)


if __name__ == "__main__":
    main()
