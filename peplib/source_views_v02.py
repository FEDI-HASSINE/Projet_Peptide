
"""Génération des vues par source du corpus v0.2.

Consommation : toute extraction destinée à la génération de séquences doit
joindre `data/processed/toxin_flags_v02.csv` par `sequence` et traiter
`toxicity_status = "unknown"` comme `"flagged"` (mode d'échec sûr, relecture
PR #19, point W4).
"""

import re
from pathlib import Path

import pandas as pd


SOURCE_ALIASES = {
    "uniprot": "UniProt",
    "uniprotkb": "UniProt",
    "amp": "UniProt",
    "bactibase": "UniProt",
    "apd3": "APD3",
    "apd6": "APD3",
    "amp_apd3": "APD3",
    "challenge": "AMPChallenge",
    "ampchallenge": "AMPChallenge",
    "amp_challenge": "AMPChallenge",
}

SOURCES = ("UniProt", "APD3", "AMPChallenge")


def normalize_sources(value):
    """Normalise les sources séparées par + ou |."""
    if pd.isna(value):
        return set()

    normalized = set()

    for raw_name in re.split(r"[+|]", str(value)):
        name = raw_name.strip()
        key = name.lower().replace(" ", "").replace("-", "_")

        if not name:
            continue

        if key not in SOURCE_ALIASES:
            raise ValueError(
                f"Source inconnue : {name!r}. "
                "Vérifier la colonne sources du corpus."
            )

        normalized.add(SOURCE_ALIASES[key])

    return normalized


def build_source_views(
    corpus_path="data/processed/corpus_v0.2.csv",
    output_dir="data/processed",
):
    """Génère une vue CSV pour chaque source du corpus."""

    corpus = pd.read_csv(corpus_path, dtype=str)

    required = {"sequence", "sources"}
    missing = required - set(corpus.columns)

    if missing:
        raise ValueError(
            f"Colonnes manquantes : {sorted(missing)}"
        )

    if corpus["sequence"].isna().any():
        raise ValueError("Séquences manquantes dans le corpus")

    if corpus["sequence"].duplicated().any():
        raise ValueError("Séquences dupliquées dans le corpus")

    source_sets = corpus["sources"].apply(normalize_sources)

    if source_sets.apply(len).eq(0).any():
        raise ValueError(
            "Certaines séquences n'ont aucune source reconnue"
        )

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    results = {}

    for source in SOURCES:
        mask = source_sets.apply(
            lambda values: source in values
        )

        view = corpus.loc[mask].copy()

        filename = (
            output_dir / f"view_{source.lower()}_v02.csv"
        )

        view.to_csv(filename, index=False)
        results[source] = len(view)

    # Vérifier que toutes les séquences sont couvertes.
    total_unique = len(corpus)
    covered = int(source_sets.apply(bool).sum())

    if covered != total_unique:
        raise AssertionError(
            "Certaines séquences ne sont pas couvertes"
        )

    print(f"Corpus total : {total_unique}")
    print(f"Séquences couvertes : {covered}")

    for source, count in results.items():
        print(f"{source}: {count} séquences")

    print(
        f"Total des appartenances : {sum(results.values())}"
    )

    print("Validation : OK")

    return results


if __name__ == "__main__":
    build_source_views()
