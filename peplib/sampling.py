
"""T57 — Echantillonnage des peptides generes.

Fonctionnalites :
- temperature
- top-k
- top-p (nucleus sampling)
- interdiction des tokens PAD et BOS
- reproductibilite avec torch.Generator
"""

import math

import torch

from peplib.tokenizer import PAD_ID, BOS_ID, VOCAB


def sampling_probabilities(
    logits: torch.Tensor,
    temperature: float = 1.0,
    top_k: int | None = None,
    top_p: float = 1.0,
) -> torch.Tensor:
    """Transforme les logits en probabilites d'echantillonnage.

    Args:
        logits: vecteur 1D contenant un score par token.
        temperature: controle la diversite (> 0).
        top_k: conserve au maximum k tokens, ou None.
        top_p: seuil de probabilite cumulee dans ]0, 1].

    Returns:
        Un vecteur de probabilites de meme taille que logits.
    """

    if logits.ndim != 1:
        raise ValueError("logits doit etre un vecteur 1D")

    if logits.numel() != len(VOCAB):
        raise ValueError(
            f"Le vocabulaire doit contenir {len(VOCAB)} tokens"
        )

    if not math.isfinite(temperature) or temperature <= 0:
        raise ValueError("temperature doit etre finie et > 0")

    if top_k is not None:
        if isinstance(top_k, bool) or not isinstance(top_k, int):
            raise ValueError("top_k doit etre un entier")
        if not 1 <= top_k <= logits.numel():
            raise ValueError("top_k invalide")

    if not math.isfinite(top_p) or not 0 < top_p <= 1:
        raise ValueError("top_p doit appartenir a ]0, 1]")

    scores = logits.clone().float()

    # Autoriser -inf pour un token deja exclu,
    # mais refuser NaN et +inf.
    if torch.isnan(scores).any() or torch.isposinf(scores).any():
        raise ValueError("logits contient NaN ou +inf")

    # PAD et BOS ne doivent jamais etre generes.
    scores[PAD_ID] = float("-inf")
    scores[BOS_ID] = float("-inf")

    if not torch.isfinite(scores).any():
        raise ValueError("Aucun token autorise disponible")

    # 1. Temperature
    scores = scores / temperature

    # 2. Top-k : selectionner au maximum k indices.
    # Cette version gere les egalites de probabilites.
    if top_k is not None:
        top_indices = torch.topk(
            scores,
            k=top_k,
        ).indices

        filtered_scores = torch.full_like(
            scores,
            float("-inf"),
        )

        filtered_scores[top_indices] = scores[top_indices]
        scores = filtered_scores

    # 3. Top-p : garder le plus petit ensemble de tokens
    # dont la probabilite cumulee atteint au moins top_p.
    if top_p < 1.0:
        sorted_scores, sorted_indices = torch.sort(
            scores,
            descending=True,
        )

        sorted_probs = torch.softmax(
            sorted_scores,
            dim=0,
        )

        cumulative = torch.cumsum(
            sorted_probs,
            dim=0,
        )

        # Conserver le token qui fait franchir le seuil.
        remove = cumulative - sorted_probs >= top_p
        remove[0] = False

        scores[sorted_indices[remove]] = float("-inf")

    probabilities = torch.softmax(scores, dim=0)

    if not torch.isfinite(probabilities).all():
        raise ValueError("Probabilites d'echantillonnage invalides")

    return probabilities


def sample_next_token(
    logits: torch.Tensor,
    temperature: float = 1.0,
    top_k: int | None = None,
    top_p: float = 1.0,
    generator: torch.Generator | None = None,
) -> int:
    """Echantillonne un prochain token, EOS inclus.

    PAD et BOS sont toujours interdits.
    """

    probabilities = sampling_probabilities(
        logits=logits,
        temperature=temperature,
        top_k=top_k,
        top_p=top_p,
    )

    token_id = torch.multinomial(
        probabilities,
        num_samples=1,
        generator=generator,
    )

    return int(token_id.item())
