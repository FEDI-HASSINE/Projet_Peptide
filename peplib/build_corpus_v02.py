
"""T37 : fusion et déduplication exacte du corpus v0.2."""

from pathlib import Path
import pandas as pd

from peplib.dedup import exact_dedup_stats


SOURCES = {
    "uniprot": "uniprot.csv",
    "amp": "amp.csv",
    "bactibase": "bactibase.csv",
    "apd3": "amp_apd3.csv",
    "challenge": "challenge.csv",
}


def build_corpus(processed_dir="data/processed"):
    processed_dir = Path(processed_dir)
    frames = []

    for name, filename in SOURCES.items():
        path = processed_dir / filename
        df = pd.read_csv(path, dtype=str)

        required = {"sequence", "source", "id"}
        if not required.issubset(df.columns):
            raise ValueError(f"Colonnes manquantes dans {filename}")

        if df["sequence"].isna().any():
            raise ValueError(f"Séquences manquantes : {filename}")

        df["sequence"] = df["sequence"].str.strip().str.upper()
        df["source"] = df["source"].fillna(name)
        df["id"] = df["id"].fillna("")

        frames.append(df)

    merged = pd.concat(frames, ignore_index=True)
    stats = exact_dedup_stats(merged["sequence"])

    # Chaque séquence est conservée une seule fois.
    # Les sources et identifiants sont regroupés.
    rows = []

    for sequence, group in merged.groupby("sequence", sort=True):
        row = {
            "sequence": sequence,
            "sources": "|".join(sorted(set(group["source"]))),
            "ids": "|".join(sorted(set(
                str(x) for x in group["id"] if x
            ))),
            "n_occurrences": len(group),
            "n_sources": group["source"].nunique(),
        }

        # Conserver les métadonnées disponibles de Challenge.
        for column in ("len", "charge", "disulfide", "activity", "dbs_amont"):
            if column in group.columns:
                values = group[column].dropna().astype(str).unique()
                row[column] = "|".join(sorted(values)) if len(values) else ""

        rows.append(row)

    corpus = pd.DataFrame(rows)
    output = processed_dir / "corpus_v0.2.csv"
    corpus.to_csv(output, index=False)

    if len(corpus) != stats["unique"]:
        raise ValueError("Incohérence du nombre de séquences uniques")

    return corpus, stats


if __name__ == "__main__":
    corpus, stats = build_corpus()
    print("Entrées brutes :", stats["total"])
    print("Séquences uniques :", stats["unique"])
    print("Doublons supprimés :", stats["duplicates"])
    print("Taux de redondance :", round(stats["redundancy_rate"] * 100, 2), "%")
    print("Fichier écrit : data/processed/corpus_v0.2.csv")
