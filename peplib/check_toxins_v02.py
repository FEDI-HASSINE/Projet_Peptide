
"""Signalement des toxines du corpus v0.2.

Ce contrôle recherche des annotations textuelles de toxicité.
L'absence de signalement ne constitue PAS une preuve de sécurité.
"""

from pathlib import Path
import pandas as pd

from peplib.biosafety import TOXIN_PATTERN


def check_toxins_v02(
    corpus_path="data/processed/corpus_v0.2.csv",
    uniprot_path="data/raw/uniprot_reviewed_20261017.tsv",
    output_path="data/processed/toxin_flags_v02.csv",
):
    corpus = pd.read_csv(corpus_path, dtype=str)
    uniprot = pd.read_csv(uniprot_path, sep="\t", dtype=str)

    required_corpus = {"sequence", "sources", "ids"}
    required_uniprot = {"Sequence", "Protein names", "Keywords"}

    if not required_corpus.issubset(corpus.columns):
        raise ValueError("Colonnes manquantes dans le corpus")

    if not required_uniprot.issubset(uniprot.columns):
        raise ValueError("Colonnes manquantes dans UniProt")

    if corpus["sequence"].isna().any():
        raise ValueError("Séquences manquantes dans le corpus")

    uniprot = uniprot.copy()
    uniprot["sequence"] = uniprot["Sequence"].fillna("").str.strip().str.upper()

    annotations = (
        uniprot["Protein names"].fillna("")
        + " "
        + uniprot["Keywords"].fillna("")
    )

    uniprot["is_toxin_annotation"] = annotations.str.contains(
        TOXIN_PATTERN, case=False, regex=True, na=False
    )

    # Plusieurs enregistrements UniProt peuvent partager une séquence.
    # Un signalement positif doit être conservé.
    flags = uniprot.groupby("sequence")["is_toxin_annotation"].any()

    result = corpus[["sequence", "sources", "ids"]].copy()
    result["has_uniprot_annotation"] = result["sequence"].isin(flags.index)

    result["is_toxin_annotation"] = (
        result["sequence"].map(flags).fillna(False).astype(bool)
    )

    result["toxicity_status"] = result.apply(
        lambda row: (
            "flagged"
            if row["is_toxin_annotation"]
            else "not_flagged_by_uniprot"
            if row["has_uniprot_annotation"]
            else "unknown"
        ),
        axis=1,
    )

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output, index=False)

    return result


if __name__ == "__main__":
    result = check_toxins_v02()

    print("Total :", len(result))
    print("Annotations UniProt disponibles :",
          int(result["has_uniprot_annotation"].sum()))
    print("Toxines signalees :",
          int(result["is_toxin_annotation"].sum()))
    print("Statuts :")
    print(result["toxicity_status"].value_counts().to_string())
    print("Sortie : data/processed/toxin_flags_v02.csv")
