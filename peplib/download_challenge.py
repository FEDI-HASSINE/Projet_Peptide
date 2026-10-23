
"""Connecteur T34 pour AMP Challenge 2027."""

from pathlib import Path
import re

import pandas as pd

from peplib.filters import filter_dataframe



def parse_header(header: str) -> dict:
    """Extrait les métadonnées d'un en-tête AMP Challenge."""
    header = header.removeprefix(">").strip()

    parts = header.split()
    if not parts:
        raise ValueError("En-tête FASTA vide")

    metadata = {"id": parts[0]}

    for part in parts[1:]:
        if "=" in part:
            key, value = part.split("=", 1)
            metadata[key.lower()] = value

    return metadata



def load_challenge(path: str | Path) -> pd.DataFrame:
    """Lit un FASTA brut sans le modifier."""
    records = []
    header = None
    sequence_parts = []

    def append_record():
        if header is not None:
            record = parse_header(header)
            record["sequence"] = "".join(sequence_parts)
            record["source"] = "AMPChallenge"
            records.append(record)

    with open(path, encoding="utf-8-sig") as file:
        for line in file:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                append_record()
                header = line
                sequence_parts = []
            elif header is not None:
                sequence_parts.append(line)
            else:
                # Certaines sources ajoutent un commentaire avant le FASTA.
                if not line.startswith((";", "#")):
                    raise ValueError("Texte inattendu avant le premier en-tête FASTA")

    append_record()

    return pd.DataFrame(
        records,
        columns=["sequence", "source", "id", "len", "charge",
                 "disulfide", "activity", "dbs"]
    )


def prepare_challenge(path: str | Path):
    """Normalise et filtre selon les règles de la Semaine 1."""
    raw = load_challenge(path)
    cleaned, rejections = filter_dataframe(raw)
    return cleaned, rejections


if __name__ == "__main__":
    path = Path("data/raw/challenge_antibacterial_20261022.fasta")
    raw = load_challenge(path)
    cleaned, rejections = prepare_challenge(path)

    print("Entrées brutes :", len(raw))
    print("Entrées conservées :", len(cleaned))
    print("Rejets :", rejections)
    print("Colonnes :", raw.columns.tolist())
    print("Premier en-tête interprété :", raw.iloc[0].to_dict())
