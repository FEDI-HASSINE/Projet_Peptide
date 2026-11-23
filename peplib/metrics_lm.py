"""T51 — Métriques des modèles de langage de peptides."""

import math
from collections import Counter

import torch
import torch.nn.functional as F

from peplib.constants import AA_LIST
from peplib.datasets import LABEL_IGNORE_INDEX
from peplib.tokenizer import PAD_ID, BOS_ID


def perplexity_from_nll(total_nll: float, n_tokens: int) -> float:
    """Calcule la perplexité à partir de la NLL totale."""
    if n_tokens <= 0:
        raise ValueError("n_tokens doit être positif")

    if not math.isfinite(total_nll) or total_nll < 0:
        raise ValueError("total_nll doit être fini et >= 0")

    return math.exp(total_nll / n_tokens)


def perplexity_from_probabilities(probabilities: list[float]) -> float:
    """Calcule la perplexité depuis les probabilités des bonnes cibles."""
    if not probabilities:
        raise ValueError("Liste de probabilités vide")

    if any(not math.isfinite(p) or not 0 < p <= 1 for p in probabilities):
        raise ValueError("Probabilités invalides")

    total_nll = -sum(math.log(p) for p in probabilities)

    return perplexity_from_nll(total_nll, len(probabilities))



def perplexity_from_logits(
    logits: torch.Tensor,
    labels: torch.Tensor,
) -> float:
    """
    Calcule la perplexité d'un modèle causal.

    logits : (batch, longueur, vocabulaire)
    labels : (batch, longueur), avec -100 pour PAD.
    """
    if logits.ndim != 3:
        raise ValueError("logits doit avoir 3 dimensions")

    if labels.shape != logits.shape[:2]:
        raise ValueError("Dimensions incompatibles")

    valid = labels != LABEL_IGNORE_INDEX

    if not valid.any().item():
        raise ValueError("Aucun token valide")

    if ((labels[valid] == PAD_ID) | (labels[valid] == BOS_ID)).any().item():
        raise ValueError("PAD et BOS ne doivent pas être des cibles")

    total_nll = F.cross_entropy(
        logits.reshape(-1, logits.shape[-1]),
        labels.reshape(-1),
        ignore_index=LABEL_IGNORE_INDEX,
        reduction="sum",
    )

    n_tokens = int(valid.sum().item())

    return perplexity_from_nll(float(total_nll.item()), n_tokens)


def positional_unigram_perplexity(
    train_sequences: list[str],
    val_sequences: list[str],
    alpha: float = 1.0,
) -> float:
    """Perplexite d'un unigramme positionnel lisse (add-alpha).

    Chaque position i a sa propre distribution sur les 20 acides amines,
    plus un symbole EOS qui cloture la sequence. C'est une reference plus
    exigeante que l'unigramme global (qui ignore la position et la
    longueur) : c'est la barre utilisee pour le gate de la Semaine 3.

    Parameters
    ----------
    train_sequences : sequences d'entrainement (10-50 residus).
    val_sequences   : sequences de validation.
    alpha           : pseudo-compte add-alpha (> 0).

    Returns
    -------
    float : perplexite en nats/token, EOS inclus.

    Notes
    -----
    L'ordre des tokens n'est pas modelise : seule la position compte. Le
    denominateur d'une position inclut les 20 acides amines (chacun avec
    alpha) et EOS (avec alpha), soit 21 * alpha de masse de lissage.
    """
    if not train_sequences:
        raise ValueError("train_sequences ne doit pas etre vide")
    if not val_sequences:
        raise ValueError("val_sequences ne doit pas etre vide")
    if not math.isfinite(alpha) or alpha <= 0:
        raise ValueError("alpha doit etre fini et > 0")

    allowed = set(AA_LIST)

    def _check(sequences: list[str], name: str) -> None:
        for index, sequence in enumerate(sequences):
            if not isinstance(sequence, str) or not sequence:
                raise ValueError(f"{name}[{index}] doit etre une chaine non vide")
            unknown = set(sequence) - allowed
            if unknown:
                raise ValueError(
                    f"{name}[{index}] contient des caracteres hors alphabet : "
                    f"{sorted(unknown)}"
                )

    _check(train_sequences, "train_sequences")
    _check(val_sequences, "val_sequences")

    max_len = max(
        max(len(s) for s in train_sequences),
        max(len(s) for s in val_sequences),
    )

    aa_counts = [Counter() for _ in range(max_len + 1)]
    eos_counts: Counter = Counter()
    for sequence in train_sequences:
        for position, residue in enumerate(sequence):
            aa_counts[position][residue] += 1
        eos_counts[len(sequence)] += 1

    total_nll = 0.0
    n_tokens = 0

    for sequence in val_sequences:
        for position, residue in enumerate(sequence):
            numerator = aa_counts[position][residue] + alpha
            denominator = (
                sum(aa_counts[position].values())
                + len(AA_LIST) * alpha
                + eos_counts[position]
                + alpha
            )
            total_nll -= math.log(numerator / denominator)
            n_tokens += 1

        position = len(sequence)
        numerator = eos_counts[position] + alpha
        denominator = (
            sum(aa_counts[position].values())
            + len(AA_LIST) * alpha
            + eos_counts[position]
            + alpha
        )
        total_nll -= math.log(numerator / denominator)
        n_tokens += 1

    return perplexity_from_nll(total_nll, n_tokens)