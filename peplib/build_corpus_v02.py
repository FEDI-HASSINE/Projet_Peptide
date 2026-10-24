"""T37 : fusion et deduplication exacte du corpus v0.2.

Trois atomes de source, conformement a la specification de la semaine 2 :
UniProt, APD3 et AMPChallenge. Les fichiers `uniprot.csv`, `amp.csv` et
`bactibase.csv` sont trois requetes UniProtKB (cf. PROVENANCE.md) : ils
alimentent tous les trois l'atome `UniProt` et ne comptent donc pas comme
des sources independantes.

Regles :
- une sequence = une ligne ;
- colonne `sources` = atomes presents, separes par `+`, tries ;
- colonne `ids` = identifiants disponibles, separes par `|`, tries ;
- colonne `longueur` = longueur de la sequence, pour toutes les lignes ;
- metadonnees Challenge conservees telles quelles (charge, disulfide,
  activity, dbs_amont) ;
- validation : sequence non vide, alphabet canonique, longueur 10-50.
  Toute entree hors de ces bornes leve une `ValueError` : le corpus ne
  peut pas contenir une sequence invalide par accident des entrees.
"""

from pathlib import Path

import pandas as pd

from peplib.dedup import exact_dedup_stats

#: Fichier d'entree -> atome de source (specification T37).
SOURCES = {
    "uniprot.csv": "UniProt",
    "amp.csv": "UniProt",
    "bactibase.csv": "UniProt",
    "amp_apd3.csv": "APD3",
    "challenge.csv": "AMPChallenge",
}

#: Atomes de source, dans l'ordre d'affichage.
ATOMES = ("UniProt", "APD3", "AMPChallenge")

ALPHABET = r"[ACDEFGHIKLMNPQRSTVWY]+"
LONGUEUR_MIN = 10
LONGUEUR_MAX = 50

METADONNEES = ("charge", "disulfide", "activity", "dbs_amont")


def _normaliser(df, filename):
    """Valide et normalise les sequences d'un fichier d'entree."""
    manquantes = {"sequence", "id"} - set(df.columns)
    if manquantes:
        raise ValueError(
            f"Colonnes manquantes dans {filename} : {sorted(manquantes)}"
        )

    sequences = df["sequence"].fillna("").str.strip().str.upper()

    vides = sequences == ""
    if vides.any():
        raise ValueError(
            f"{filename} : {int(vides.sum())} sequence(s) vide(s)"
        )

    hors_alphabet = ~sequences.str.fullmatch(ALPHABET)
    if hors_alphabet.any():
        exemple = sequences[hors_alphabet].iloc[0]
        raise ValueError(
            f"{filename} : {int(hors_alphabet.sum())} sequence(s) hors "
            f"alphabet canonique (ex. {exemple!r})"
        )

    longueurs = sequences.str.len()
    hors_fenetre = (longueurs < LONGUEUR_MIN) | (longueurs > LONGUEUR_MAX)
    if hors_fenetre.any():
        exemple = sequences[hors_fenetre].iloc[0]
        raise ValueError(
            f"{filename} : {int(hors_fenetre.sum())} sequence(s) hors "
            f"fenetre {LONGUEUR_MIN}-{LONGUEUR_MAX} (ex. longueur "
            f"{len(exemple)})"
        )

    df = df.copy()
    df["sequence"] = sequences
    return df


def _joindre(merged, colonne, separateur):
    """Valeurs distinctes triees par sequence, jointes par `separateur`."""
    from itertools import groupby

    sous = merged[["sequence", colonne]].dropna()
    sous = sous[sous[colonne].astype(str) != ""]
    sous = sous.drop_duplicates().sort_values(["sequence", colonne])

    valeurs = {}
    for sequence, groupe in groupby(
        zip(sous["sequence"], sous[colonne]), key=lambda paire: paire[0]
    ):
        valeurs[sequence] = separateur.join(valeur for _, valeur in groupe)
    return valeurs


def build_corpus(processed_dir="data/processed"):
    """Fusionne les cinq fichiers d'entree en un corpus deduplique."""
    processed_dir = Path(processed_dir)
    frames = []

    for filename, atome in SOURCES.items():
        df = pd.read_csv(processed_dir / filename, dtype=str)
        df = _normaliser(df, filename)
        df["atome"] = atome
        df["id"] = df["id"].fillna("")
        frames.append(df)

    merged = pd.concat(frames, ignore_index=True)
    stats = exact_dedup_stats(merged["sequence"])

    tailles = merged.groupby("sequence", sort=True).size()
    atomes = merged.groupby("sequence", sort=True)["atome"].nunique()
    index = tailles.index

    corpus = pd.DataFrame(index=index)
    corpus["sources"] = index.map(_joindre(merged, "atome", "+"))
    corpus["ids"] = index.map(_joindre(merged, "id", "|"))
    corpus["n_occurrences"] = tailles
    corpus["n_sources"] = atomes
    corpus["longueur"] = index.str.len()

    for colonne in METADONNEES:
        if colonne in merged.columns:
            corpus[colonne] = index.map(_joindre(merged, colonne, "|")).fillna("")

    corpus = corpus.reset_index()

    ordre = ["sequence", "sources", "ids", "n_occurrences", "n_sources",
             "longueur"] + [c for c in METADONNEES if c in corpus.columns]
    corpus = corpus[ordre]

    if len(corpus) != stats["unique"]:
        raise ValueError("Incoherence du nombre de sequences uniques")
    if not corpus["sequence"].is_unique:
        raise ValueError("Sequences dupliquees dans le corpus")

    corpus.to_csv(processed_dir / "corpus_v0.2.csv", index=False)
    return corpus, stats


if __name__ == "__main__":
    corpus, stats = build_corpus()
    print("Entrees brutes :", stats["total"])
    print("Sequences uniques :", stats["unique"])
    print("Doublons supprimes :", stats["duplicates"])
    print("Taux de redondance :", round(stats["redundancy_rate"] * 100, 2), "%")
    print("Fichier ecrit : data/processed/corpus_v0.2.csv")
