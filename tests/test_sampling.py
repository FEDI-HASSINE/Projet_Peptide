
"""Tests T57 : temperature, top-k, top-p et baselines n-grammes."""

import pytest
import torch

from peplib.ngram import NgramModel, TARGETS
from peplib.tokenizer import VOCAB, PAD_ID, BOS_ID, EOS_ID

from peplib.sampling import (
    sampling_probabilities,
    sample_next_token,
)


def test_forbidden_tokens():
    logits = torch.zeros(23)

    probs = sampling_probabilities(logits)

    assert probs[PAD_ID].item() == 0
    assert probs[BOS_ID].item() == 0
    assert probs.sum().item() == pytest.approx(1.0)


def test_top_k():
    logits = torch.arange(23, dtype=torch.float32)

    probs = sampling_probabilities(logits, top_k=2)

    assert int((probs > 0).sum()) == 2


def test_top_p():
    logits = torch.full((23,), -10.0)
    logits[3] = 3.0
    logits[4] = 2.0
    logits[5] = 1.0

    probs = sampling_probabilities(
        logits,
        top_p=0.8,
    )

    assert probs.sum().item() == pytest.approx(1.0)
    assert probs[3].item() > 0


def test_reproducibility():
    logits = torch.zeros(23)

    g1 = torch.Generator().manual_seed(42)
    g2 = torch.Generator().manual_seed(42)

    result1 = sample_next_token(
        logits,
        generator=g1,
    )

    result2 = sample_next_token(
        logits,
        generator=g2,
    )

    assert result1 == result2


def test_invalid_temperature():
    with pytest.raises(ValueError):
        sampling_probabilities(
            torch.zeros(23),
            temperature=0,
        )


def test_top_k_one():
    logits = torch.zeros(23)
    logits[8] = 10.0

    probs = sampling_probabilities(
        logits,
        top_k=1,
    )

    assert probs[8].item() == pytest.approx(1.0)
    assert int((probs > 0).sum()) == 1


def test_temperature_effect():
    logits = torch.zeros(23)
    logits[3] = 3.0
    logits[4] = 1.0

    cold = sampling_probabilities(
        logits,
        temperature=0.5,
    )

    hot = sampling_probabilities(
        logits,
        temperature=2.0,
    )

    assert cold[3] > hot[3]


def test_forbidden_tokens_even_if_high_logits():
    logits = torch.zeros(23)
    logits[PAD_ID] = 1000.0
    logits[BOS_ID] = 1000.0

    probs = sampling_probabilities(logits)

    assert probs[PAD_ID].item() == 0
    assert probs[BOS_ID].item() == 0
    assert probs.sum().item() == pytest.approx(1.0)


@pytest.mark.parametrize(
    "order,context",
    [
        (1, ()),
        (2, ("A",)),
    ],
)
def test_sampling_with_ngram(order, context):
    """Integration avec unigramme et bigramme."""

    sequences = [
        "ACACAC",
        "ACAC",
        "AGAC",
        "AACG",
    ]

    model = NgramModel(
        order=order,
        alpha=1.0,
    )
    model.fit(sequences)

    # Convertir les log-probabilites du modele
    # en logits du vocabulaire de 23 tokens.
    logits = torch.full(
        (len(VOCAB),),
        float("-inf"),
    )

    for token in TARGETS:
        token_id = VOCAB[token]

        logits[token_id] = model.log_prob_target(
            token,
            context=context,
        )

    probs = sampling_probabilities(
        logits,
        temperature=1.0,
        top_k=5,
        top_p=0.9,
    )

    assert torch.isfinite(probs).all()
    assert probs.sum().item() == pytest.approx(1.0)

    assert probs[PAD_ID].item() == 0.0
    assert probs[BOS_ID].item() == 0.0

    assert 1 <= int((probs > 0).sum()) <= 5

    # Verifier la reproductibilite du tirage.
    g1 = torch.Generator().manual_seed(42)
    g2 = torch.Generator().manual_seed(42)

    token1 = sample_next_token(
        logits,
        top_k=5,
        top_p=0.9,
        generator=g1,
    )

    token2 = sample_next_token(
        logits,
        top_k=5,
        top_p=0.9,
        generator=g2,
    )

    assert token1 == token2
    assert token1 not in (PAD_ID, BOS_ID)


def test_top_k_with_equal_scores():
    """Top-k respecte k meme si les scores sont egaux."""

    logits = torch.zeros(23)

    probs = sampling_probabilities(
        logits,
        top_k=5,
        top_p=1.0,
    )

    assert int((probs > 0).sum()) == 5
    assert probs.sum().item() == pytest.approx(1.0)

    assert probs[PAD_ID].item() == 0
    assert probs[BOS_ID].item() == 0


def test_top_p_keeps_threshold_crossing_token():
    """Top-p conserve le token qui franchit le seuil."""

    logits = torch.full((23,), float("-inf"))

    logits[3] = 0.0
    logits[4] = 0.0
    logits[5] = 0.0

    probs = sampling_probabilities(
        logits,
        top_p=0.5,
    )

    assert int((probs > 0).sum()) == 2
    assert probs.sum().item() == pytest.approx(1.0)


def test_only_eos_available():
    """EOS reste generable meme si tous les AA sont exclus."""

    logits = torch.full((23,), float("-inf"))
    logits[EOS_ID] = 5.0

    probs = sampling_probabilities(logits)

    assert probs[EOS_ID].item() == pytest.approx(1.0)

    token = sample_next_token(logits)
    assert token == EOS_ID


def test_all_allowed_tokens_masked():
    """Refuser une distribution sans aucun token possible."""

    logits = torch.full((23,), float("-inf"))

    with pytest.raises(ValueError):
        sampling_probabilities(logits)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"temperature": -1.0},
        {"temperature": float("nan")},
        {"temperature": float("inf")},
        {"top_k": 0},
        {"top_k": 24},
        {"top_p": 0.0},
        {"top_p": 1.5},
    ],
)
def test_invalid_sampling_parameters(kwargs):
    with pytest.raises(ValueError):
        sampling_probabilities(
            torch.zeros(23),
            **kwargs,
        )


def test_invalid_logit_dimensions():
    with pytest.raises(ValueError):
        sampling_probabilities(torch.zeros(2, 23))


def test_invalid_vocab_size():
    with pytest.raises(ValueError):
        sampling_probabilities(torch.zeros(21))
