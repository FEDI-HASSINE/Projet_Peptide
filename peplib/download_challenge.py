"""Connecteur T34 pour le fichier AMP Challenge.

Ce module ne télécharge rien : il lit un FASTA déjà présent dans
``data/raw/`` (immuable) et écrit la vue filtrée dans ``data/processed/``.

Les en-têtes du Challenge portent les métadonnées
(``len=``, ``charge=``, ``disulfide=``, ``activity=``, ``dbs=``). Le champ
``dbs=`` est conservé en colonne ``dbs_amont`` : c'est une information de
provenance multi-bases (13 bases possibles), il ne pilote ni filtre ni
split. Tout jeton d'en-tête qui n'est pas au format ``clé=valeur`` est
ignoré, et une clé absente produit une valeur ``NaN`` (aucune valeur
inventée).

Usage :
    python -m peplib.download_challenge
"""

from pathlib import Path

import pandas as pd

from peplib.filters import filter_dataframe


SOURCE = "AMPChallenge"
DEFAULT_RAW = Path("data/raw/challenge_antibacterial_20261022.fasta")
DEFAULT_OUT = Path("data/processed/challenge.csv")
COLUMNS = [
    "sequence", "source", "id", "len", "charge",
    "disulfide", "activity", "dbs_amont",
]


def parse_header(header: str) -> dict:
    """Extrait les métadonnées d'un en-tête AMP Challenge.

    Le premier jeton est l'identifiant ; les suivants sont lus au format
    ``clé=valeur``. Un jeton sans ``=`` est ignoré. Un en-tête vide lève
    une ``ValueError``.
    """
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
    """Lit un FASTA brut sans le modifier.

    Les lignes vides sont ignorées. Un commentaire (``;`` ou ``#``) est
    toléré avant le premier en-tête ; tout autre texte à cet endroit lève
    une ``ValueError``.
    """
    records = []
    header = None
    sequence_parts = []

    def append_record():
        if header is not None:
            record = parse_header(header)
            record["sequence"] = "".join(sequence_parts)
            record["source"] = SOURCE
            record["dbs_amont"] = record.pop("dbs", None)
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
                    raise ValueError(
                        "Texte inattendu avant le premier en-tête FASTA"
                    )

    append_record()

    return pd.DataFrame(records, columns=COLUMNS)


def prepare_challenge(path: str | Path):
    """Normalise et filtre selon les règles de la Semaine 1 (10-50, 20 AA)."""
    raw = load_challenge(path)
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
    raw = load_challenge(raw_path)
    cleaned, rejections = prepare_challenge(raw_path)
    written = write_processed(cleaned, out_path)

    print("Source :", SOURCE)
    print("Entrées brutes :", len(raw))
    print("Entrées conservées :", len(cleaned))
    print("Rejets :", rejections)
    print("Colonnes :", cleaned.columns.tolist())
    print("Fichier écrit :", written)
    return written


if __name__ == "__main__":
    main()
