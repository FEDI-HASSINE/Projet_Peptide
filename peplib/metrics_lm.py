"""T51 — Métriques des modèles de langage de peptides."""

import math
import torch
import torch.nn.functional as F

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