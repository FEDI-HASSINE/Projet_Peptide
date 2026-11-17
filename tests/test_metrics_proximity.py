
"""Tests T60 — Proximite des peptides generes avec le train."""

import pytest

from peplib.identite import identite_couverture
from peplib.metrics_proximity import (
    proximity_details,
    proximity_summary,
)


def test_identical_peptide():
    details = proximity_details(
        ["ACDEFG"],
        ["ACDEFG"],
    )

    assert details.iloc[0]["best_min"] == pytest.approx(1.0)
    assert bool(details.iloc[0]["near_train_90"])


def test_substring_of_train():
    """Un fragment inclus dans le train est tres proche localement."""
    details = proximity_details(
        ["ACDE"],
        ["GGACDEGG"],
    )

    assert details.iloc[0]["best_min"] == pytest.approx(1.0)


def test_proximity_summary():
    generated = ["ACDE", "ACDE", "GGGG"]
    train = ["ACDE"]

    result = proximity_summary(generated, train)

    assert result["n_generated"] == 3
    assert result["n_train"] == 1
    assert result["n_near_train"] >= 2
    assert result["fraction_near_train"] >= 2 / 3


def test_short_peptides():
    """Les sequences plus courtes que 3 AA restent comparables."""
    details = proximity_details(
        ["AC"],
        ["AC"],
    )

    assert details.iloc[0]["best_min"] == pytest.approx(1.0)


def test_empty_generated():
    with pytest.raises(ValueError):
        proximity_details([], ["ACDE"])


def test_empty_train():
    with pytest.raises(ValueError):
        proximity_details(["ACDE"], [])


def test_invalid_sequence():
    with pytest.raises(ValueError):
        proximity_details(["ACX"], ["ACDE"])


def test_invalid_threshold():
    with pytest.raises(ValueError):
        proximity_details(
            ["ACDE"],
            ["ACDE"],
            threshold=1.5,
        )


def test_threshold_matches_exhaustive_on_small_corpus():
    """Le prefiltre ne doit pas manquer les paires >= 90 %."""

    generated = [
        "ACDEFG",
        "ACDEFA",
        "GGGGGG",
        "ACDE",
    ]

    train = [
        "ACDEFG",
        "GGACDEGG",
        "WWWWWW",
    ]

    details = proximity_details(
        generated,
        train,
        threshold=0.90,
    )

    for i, seq_gen in enumerate(generated):
        # Comparaison exhaustive : oracle de test.
        best_exhaustive = max(
            min(*identite_couverture(seq_gen, seq_train))
            for seq_train in train
        )

        assert bool(details.iloc[i]["near_train_90"]) == (
            best_exhaustive >= 0.90
        )
