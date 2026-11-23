import math
import pytest
import torch

from peplib.metrics_lm import perplexity_from_logits

from peplib.metrics_lm import (
    perplexity_from_nll,
    perplexity_from_probabilities,
)


def test_uniform_perplexity_21():
    probabilities = [1 / 21] * 100
    assert perplexity_from_probabilities(probabilities) == pytest.approx(21)


def test_perfect_predictions():
    assert perplexity_from_probabilities([1.0, 1.0]) == pytest.approx(1)


def test_perplexity_from_nll():
    nll = 100 * math.log(21)
    assert perplexity_from_nll(nll, 100) == pytest.approx(21)


def test_invalid_probabilities():
    with pytest.raises(ValueError):
        perplexity_from_probabilities([0.0, 0.5])


def test_invalid_token_count():
    with pytest.raises(ValueError):
        perplexity_from_nll(1.0, 0)
        

def test_uniform_logits_perplexity_21():
    # 23 tokens au total, mais PAD et BOS sont interdits en sortie.
    logits = torch.zeros(2, 4, 23)
    logits[:, :, 0] = float("-inf")
    logits[:, :, 1] = float("-inf")

    labels = torch.tensor([
        [3, 4, 2, -100],
        [5, 2, -100, -100],
    ])

    result = perplexity_from_logits(logits, labels)

    assert result == pytest.approx(21.0)


def test_perplexity_ignore_padding():
    logits = torch.zeros(1, 4, 23)
    logits[:, :, 0:2] = float("-inf")

    labels = torch.tensor([[3, 4, 2, -100]])

    result = perplexity_from_logits(logits, labels)

    assert result == pytest.approx(21.0)


def test_perplexity_reject_bos_target():
    logits = torch.zeros(1, 2, 23)
    labels = torch.tensor([[1, 2]])

    with pytest.raises(ValueError):
        perplexity_from_logits(logits, labels)        
# ----------------------------------------------------------------------
# Unigramme positionnel (barre du gate Semaine 3)
# ----------------------------------------------------------------------

def test_positional_unigram_toy_hand_checked():
    """Corpus ["AA"] : p(A|0)=p(A|1)=p(EOS|2)=2/22=1/11 -> perplexite 11."""
    from peplib.metrics_lm import positional_unigram_perplexity

    value = positional_unigram_perplexity(["AA"], ["AA"], alpha=1.0)
    assert value == pytest.approx(11.0)


def test_positional_unigram_is_finite_on_small_corpus():
    from peplib.metrics_lm import positional_unigram_perplexity

    train = ["ACDE", "ACDF", "GGGG"]
    val = ["ACDE", "GGGG"]
    value = positional_unigram_perplexity(train, val, alpha=1.0)
    assert math.isfinite(value)
    assert value > 1.0


def test_positional_unigram_rejects_bad_input():
    from peplib.metrics_lm import positional_unigram_perplexity

    with pytest.raises(ValueError):
        positional_unigram_perplexity([], ["ACDE"])
    with pytest.raises(ValueError):
        positional_unigram_perplexity(["ACDE"], [])
    with pytest.raises(ValueError):
        positional_unigram_perplexity(["ACDE"], ["ACDE"], alpha=0.0)
    with pytest.raises(ValueError):
        positional_unigram_perplexity(["ACDE"], ["ACXDE"])
