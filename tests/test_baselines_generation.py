
"""Tests T61 — Baselines de reference."""

from collections import Counter

import pytest

from peplib.constants import AA_LIST
from peplib.ngram import NgramModel
from peplib.baselines_generation import (
    random_baseline,
    validation_baseline,
    shuffled_validation_baseline,
    generate_ngram,
    build_reference_baselines,
)


TRAIN = [
    "ACDEFG",
    "ACACAC",
    "GGGGGG",
    "WYACDE",
]

VAL = [
    "ACDE",
    "GGGG",
    "ACAC",
    "WY",
]


def test_random_lengths():
    generated = random_baseline(VAL)

    assert [len(s) for s in generated] == [
        len(s) for s in VAL
    ]


def test_random_validity():
    generated = random_baseline(VAL)

    assert all(
        set(s).issubset(set(AA_LIST))
        for s in generated
    )


def test_random_reproducible():
    assert random_baseline(VAL, 42) == random_baseline(VAL, 42)


def test_validation_baseline():
    result = validation_baseline(VAL)

    assert result == VAL
    assert result is not VAL


def test_shuffle_preserves_composition():
    shuffled = shuffled_validation_baseline(VAL)

    for original, modified in zip(VAL, shuffled):
        assert Counter(original) == Counter(modified)
        assert len(original) == len(modified)


def test_shuffle_reproducible():
    result1 = shuffled_validation_baseline(VAL, 42)
    result2 = shuffled_validation_baseline(VAL, 42)

    assert result1 == result2


@pytest.mark.parametrize("order", [1, 2])
def test_generate_ngram(order):
    model = NgramModel(order=order, alpha=1.0).fit(TRAIN)

    generated = generate_ngram(
        model,
        n_sequences=10,
        seed=42,
    )

    assert len(generated) == 10
    assert all(10 <= len(s) <= 50 for s in generated)
    assert all(set(s).issubset(set(AA_LIST)) for s in generated)


@pytest.mark.parametrize("order", [1, 2])
def test_generate_ngram_reproducible(order):
    model = NgramModel(order=order, alpha=1.0).fit(TRAIN)

    first = generate_ngram(
        model,
        n_sequences=10,
        seed=42,
    )

    second = generate_ngram(
        model,
        n_sequences=10,
        seed=42,
    )

    assert first == second


def test_build_five_baselines():
    result = build_reference_baselines(TRAIN, VAL, seed=42)

    assert set(result) == {
        "random",
        "unigram",
        "bigram",
        "validation",
        "shuffled_validation",
    }

    assert all(len(seqs) == len(VAL) for seqs in result.values())


def test_invalid_ngram_order():
    model = NgramModel(order=3, alpha=1.0).fit(TRAIN)

    with pytest.raises(ValueError):
        generate_ngram(model, n_sequences=10)


def test_empty_train():
    with pytest.raises(ValueError):
        build_reference_baselines([], VAL)


@pytest.mark.parametrize("order", [1, 2])
def test_ngram_minimum_length_10(order):
    """B2 : tous les peptides generes ont au moins 10 AA."""

    model = NgramModel(order=order, alpha=1.0).fit(TRAIN)

    generated = generate_ngram(
        model,
        n_sequences=100,
        seed=42,
    )

    assert len(generated) == 100
    assert all(10 <= len(s) <= 50 for s in generated)


@pytest.mark.parametrize("order", [1, 2])
def test_ngram_minimum_max_residues(order):
    """Un maximum de 10 doit produire des peptides de 10 AA."""

    model = NgramModel(order=order, alpha=1.0).fit(TRAIN)

    generated = generate_ngram(
        model,
        n_sequences=20,
        seed=42,
        max_residues=10,
    )

    assert len(generated) == 20
    assert all(len(s) == 10 for s in generated)


def test_ngram_rejects_max_length_below_10():
    model = NgramModel(order=1, alpha=1.0).fit(TRAIN)

    with pytest.raises(ValueError):
        generate_ngram(
            model,
            n_sequences=5,
            max_residues=9,
        )
