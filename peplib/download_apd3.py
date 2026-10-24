"""Connecteur T34 pour le fichier APD 2024 (site APD6).

Ce module ne télécharge rien : il lit un FASTA déjà présent dans
``data/raw/`` (immuable) et écrit la vue filtrée dans ``data/processed/``.

Le FASTA APD ne contient que l'identifiant ``APxxxxx`` et la séquence :
les métadonnées (nom, organisme, cible, activité, unité) ne sont pas dans
ce fichier. Aucune colonne n'est fabriquée à vide ; la dérogation est
actée par écrit dans ``PROVENANCE.md`` (bloc T34).

Usage :
    python -m peplib.download_apd3
"""

from pathlib import Path
import re

import pandas as pd

from peplib.filters import filter_dataframe


APD_ID = re.compile(r"^AP\d+$")
SOURCE = "APD3"
DEFAULT_RAW = Path("data/raw/naturalAMPs_APD2024a.fasta")
DEFAULT_OUT = Path("data/processed/amp_apd3.csv")


def load_apd3(path: str | Path) -> pd.DataFrame:
    """Lit les peptides APD sans modifier le fichier brut.

    Un en-tête qui n'est pas un identifiant ``APxxxxx`` est ignoré : c'est
    le cas de la ligne d'annonce du site (``>Your search led to 3306
    peptides``), qui n'est pas une séquence. Les lignes vides séparant les
    enregistrements sont ignorées.
    """
    records = []
    peptide_id = None
    sequence_parts = []

    def append_record():
        if peptide_id is not None:
            records.append({
                "sequence": "".join(sequence_parts),
                "source": SOURCE,
                "id": peptide_id,
            })

    with open(path, encoding="utf-8-sig") as file:
        for line in file:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                append_record()
                header = line[1:].strip()
                peptide_id = header if APD_ID.fullmatch(header) else None
                sequence_parts = []
            elif peptide_id is not None:
                sequence_parts.append(line)

    append_record()

    return pd.DataFrame(records, columns=["sequence", "source", "id"])


def prepare_apd3(path: str | Path):
    """Applique les règles de filtrage communes (alphabet 20 AA, 10-50)."""
    raw = load_apd3(path)
    cleaned, rejections = filter_dataframe(raw)
    return cleaned, rejections


def write_processed(df: pd.DataFrame, path: str | Path = DEFAULT_OUT) -> Path:
    """Écrit la vue filtrée en UTF-8 et retourne le chemin écrit."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, encoding="utf-8")
    return path


def main(raw_path: str | Path = DEFAULT_RAW,
         out_path: str | Path = DEFAULT_OUT) -> Path:
    """Point d'entrée : parse, filtre, écrit, affiche le bilan."""
    raw = load_apd3(raw_path)
    cleaned, rejections = prepare_apd3(raw_path)
    written = write_processed(cleaned, out_path)

    print("Source :", SOURCE)
    print("Entrées brutes :", len(raw))
    print("Entrées conservées :", len(cleaned))
    print("Rejets :", rejections)
    print("Fichier écrit :", written)
    return written


if __name__ == "__main__":
    main()
