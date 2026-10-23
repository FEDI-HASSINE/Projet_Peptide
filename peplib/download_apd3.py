
"""Connecteur T34 pour le fichier APD 2024 (site APD6)."""

from pathlib import Path
import re

import pandas as pd

from peplib.filters import filter_dataframe


APD_ID = re.compile(r"^AP\d+$")


def load_apd3(path: str | Path) -> pd.DataFrame:
    """Lit les peptides APD sans modifier le fichier brut."""
    records = []
    peptide_id = None
    sequence_parts = []

    def append_record():
        if peptide_id is not None:
            records.append({
                "id": peptide_id,
                "sequence": "".join(sequence_parts),
                "source": "APD",
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
    """Applique les règles de filtrage communes."""
    raw = load_apd3(path)
    cleaned, rejections = filter_dataframe(raw)
    return cleaned, rejections


if __name__ == "__main__":
    path = Path("data/raw/naturalAMPs_APD2024a.fasta")

    raw = load_apd3(path)
    cleaned, rejections = prepare_apd3(path)

    print("Entrées brutes :", len(raw))
    print("Entrées conservées :", len(cleaned))
    print("Rejets :", rejections)
    print("Premier peptide :", raw.iloc[0].to_dict())
