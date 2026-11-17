
"""T60 — Proximite des peptides generes avec le train.

Reutilise les fonctions de la semaine 2 :
- candidats_val_train : prefiltrage par k-mers
- meilleures_paires : alignement Smith-Waterman
- identite_couverture : convention identite + couverture

Aucun acces au split de test.
"""

import pandas as pd

from peplib.constants import AA_LIST
from peplib.fuite import candidats_val_train, meilleures_paires


AA_SET = set(AA_LIST)


def _validate_sequences(sequences, name):
    """Verifie le format des sequences d'entree."""
    if not sequences:
        raise ValueError(f"{name} ne doit pas etre vide")

    for i, sequence in enumerate(sequences):
        if not isinstance(sequence, str):
            raise ValueError(f"{name}[{i}] n'est pas une chaine")

        if not 1 <= len(sequence) <= 50:
            raise ValueError(
                f"{name}[{i}] doit avoir entre 1 et 50 residus"
            )

        if any(aa not in AA_SET for aa in sequence):
            raise ValueError(
                f"{name}[{i}] contient des acides amines invalides"
            )


def proximity_details(
    generated: list[str],
    train: list[str],
    threshold: float = 0.90,
) -> pd.DataFrame:
    """Calcule la meilleure proximite de chaque peptide avec le train.

    La proximite est mesuree par :
        best_min = max_train(min(identite_SW, couverture_SW))

    Attention :
    le prefiltrage vise la detection des paires au seuil demande.
    Les scores sous ce seuil peuvent etre sous-estimes.
    """

    _validate_sequences(generated, "generated")
    _validate_sequences(train, "train")

    if not 0 < threshold <= 1:
        raise ValueError("threshold doit etre dans ]0, 1]")

    # Recherche de candidats croises genere -> train.
    candidates = candidats_val_train(
        generated,
        train,
        tau=threshold,
        k=3,
    )

    # Pour les peptides plus courts que 3 residus,
    # le filtre 3-mers n'est pas applicable.
    # On compare donc ces paires exhaustivement.
    for i, seq_gen in enumerate(generated):
        if len(seq_gen) < 3:
            candidates.update(
                (i, j) for j in range(len(train))
            )

    for j, seq_train in enumerate(train):
        if len(seq_train) < 3:
            candidates.update(
                (i, j) for i in range(len(generated))
            )

    # Fonction de l'equipe : calcul SW sur les paires retenues.
        details = meilleures_paires(
        generated,
        train,
        candidates,
    ).copy()

    # meilleures_paires fournit deja la colonne sequence.
    # On evite donc de l'ajouter une deuxieme fois.
    if "sequence" not in details.columns:
        details.insert(0, "sequence", generated)

    details["near_train_90"] = (
        details["best_min"] >= 0.90
    )

    details["near_train_threshold"] = (
        details["best_min"] >= threshold
    )

    return details


def proximity_summary(
    generated: list[str],
    train: list[str],
    threshold: float = 0.90,
) -> dict:
    """Resume la proximite des peptides generes avec le train."""

    details = proximity_details(
        generated,
        train,
        threshold=threshold,
    )

    near_count = int(
        details["near_train_threshold"].sum()
    )

    return {
        "n_generated": len(generated),
        "n_train": len(train),
        "threshold": threshold,
        "n_near_train": near_count,
        "fraction_near_train": near_count / len(generated),
        "fraction_near_train_90": float(
            details["near_train_90"].mean()
        ),
        "mean_best_min_prefiltered": float(
            details["best_min"].mean()
        ),
        "max_best_min_prefiltered": float(
            details["best_min"].max()
        ),
    }
