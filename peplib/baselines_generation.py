
"""T61 — Baselines de reference pour la generation de peptides.

Baselines :
- aleatoire uniforme
- unigramme
- bigramme
- validation reelle
- validation melangee

Les modeles n-grammes sont ajustes uniquement sur le train.

Correction B2 (PR #26) :
- Longueur minimale : 10 AA.
- Longueur maximale : 50 AA.
- EOS interdit avant la generation de 10 AA.
"""

import csv
import random
from pathlib import Path

import torch

from peplib.constants import AA_LIST
from peplib.ngram import NgramModel, TARGETS
from peplib.sampling import sample_next_token
from peplib.tokenizer import VOCAB, EOS_ID, ID_TO_TOKEN


def load_sequences(csv_path: str | Path) -> list[str]:
    """Charge les peptides d'un CSV contenant la colonne sequence."""

    with open(csv_path, newline="", encoding="utf-8-sig") as file:
        reader = csv.DictReader(file)

        if "sequence" not in (reader.fieldnames or []):
            raise ValueError("Colonne sequence absente")

        sequences = [
            row["sequence"].strip()
            for row in reader
            if row.get("sequence") and row["sequence"].strip()
        ]

    if not sequences:
        raise ValueError("Aucune sequence disponible")

    return sequences


def random_baseline(
    reference: list[str],
    seed: int = 42,
) -> list[str]:
    """Genere des peptides aleatoires aux longueurs des references."""

    rng = random.Random(seed)

    return [
        "".join(rng.choices(AA_LIST, k=len(seq)))
        for seq in reference
    ]


def validation_baseline(
    validation: list[str],
) -> list[str]:
    """Copie les sequences reelles de validation."""

    return list(validation)


def shuffled_validation_baseline(
    validation: list[str],
    seed: int = 42,
) -> list[str]:
    """Melange les AA de chaque peptide sans changer sa composition."""

    rng = random.Random(seed)
    result = []

    for sequence in validation:
        letters = list(sequence)
        rng.shuffle(letters)
        result.append("".join(letters))

    return result


def generate_ngram(
    model: NgramModel,
    n_sequences: int,
    seed: int = 42,
    max_residues: int = 50,
    temperature: float = 1.0,
    top_k: int | None = None,
    top_p: float = 1.0,
) -> list[str]:
    """Genere des peptides avec un unigramme ou bigramme ajuste.

    Regles de generation :
    - Chaque peptide contient entre 10 et 50 AA.
    - EOS termine la generation apres au moins 10 AA.
    - EOS est interdit lors des 10 premiers tirages.
    - Si EOS n'apparait pas, arret force a max_residues.
    """

    if model.order not in (1, 2):
        raise ValueError(
            "Seuls unigramme et bigramme sont supportes"
        )

    if n_sequences <= 0:
        raise ValueError("n_sequences doit etre positif")

    # B2 : le maximum autorise ne peut pas etre inferieur a 10.
    if not 10 <= max_residues <= 50:
        raise ValueError(
            "max_residues doit etre entre 10 et 50"
        )

    generator = torch.Generator().manual_seed(seed)
    generated = []

    for _ in range(n_sequences):
        residues = []

        for _step in range(max_residues):

            # Unigramme : pas de contexte.
            # Bigramme : dernier AA ou BOS au debut.
            if model.order == 1:
                context = ()
            else:
                context = (
                    residues[-1] if residues else "<BOS>",
                )

            logits = torch.full(
                (len(VOCAB),),
                float("-inf"),
            )

            # Probabilites des 20 AA et de EOS.
            for token in TARGETS:
                logits[VOCAB[token]] = model.log_prob_target(
                    token,
                    context=context,
                )

            # CORRECTION B2 :
            # EOS ne peut pas etre genere avant 10 AA.
            if len(residues) < 10:
                logits[EOS_ID] = float("-inf")

            token_id = sample_next_token(
                logits,
                temperature=temperature,
                top_k=top_k,
                top_p=top_p,
                generator=generator,
            )

            # EOS termine la generation.
            if token_id == EOS_ID:
                break

            # Sinon, ajouter l'acide amine genere.
            residues.append(ID_TO_TOKEN[token_id])

        generated.append("".join(residues))

    return generated


def build_reference_baselines(
    train: list[str],
    validation: list[str],
    seed: int = 42,
    alpha: float = 1.0,
) -> dict[str, list[str]]:
    """Construit les cinq references de T61."""

    if not train or not validation:
        raise ValueError(
            "Train et validation doivent etre non vides"
        )

    # Apprentissage uniquement sur le train.
    unigram = NgramModel(
        order=1,
        alpha=alpha,
    ).fit(train)

    bigram = NgramModel(
        order=2,
        alpha=alpha,
    ).fit(train)

    n = len(validation)

    return {
        "random": random_baseline(
            validation,
            seed=seed,
        ),
        "unigram": generate_ngram(
            unigram,
            n_sequences=n,
            seed=seed,
        ),
        "bigram": generate_ngram(
            bigram,
            n_sequences=n,
            seed=seed,
        ),
        "validation": validation_baseline(
            validation,
        ),
        "shuffled_validation": shuffled_validation_baseline(
            validation,
            seed=seed,
        ),
    }
