
"""Tests T59 — Qualite statistique des peptides generes."""

import math

import pytest

from peplib.metrics_generation import (
    is_valid_peptide,
    validity_rate,
    uniqueness_rate,
    sequence_length_stats,
    amino_acid_distribution,
    js_divergence,
    ks_length_statistic,
    generation_metrics,
)


def test_valid_peptide():
    assert is_valid_peptide("ACDEFG") is True


def test_invalid_peptide():
    assert is_valid_peptide("ACX") is False
    assert is_valid_peptide("acde") is False
    assert is_valid_peptide("") is False


def test_validity_rate():
    sequences = ["ACD", "GGG", "ACX", "WY"]

    assert validity_rate(sequences) == pytest.approx(0.75)


def test_uniqueness_rate():
    sequences = ["ACD", "ACD", "GGG", "WY"]

    assert uniqueness_rate(sequences) == pytest.approx(0.75)


def test_lengths():
    sequences = ["AC", "ACDE", "GGGGGG"]

    stats = sequence_length_stats(sequences)

    assert stats["mean"] == pytest.approx(4.0)
    assert stats["min"] == 2
    assert stats["max"] == 6
    assert stats["median"] == pytest.approx(4.0)


def test_aa_distribution_sums_to_one():
    distribution = amino_acid_distribution(["ACD", "GGG"])

    assert len(distribution) == 20
    assert distribution.sum() == pytest.approx(1.0)


def test_js_identical():
    sequences = ["ACDE", "GGGG", "WY"]

    assert js_divergence(sequences, sequences) == pytest.approx(0.0)


def test_js_different():
    generated = ["AAAA"]
    reference = ["CCCC"]

    assert js_divergence(generated, reference) == pytest.approx(
        math.log(2)
    )


def test_ks_identical_lengths():
    generated = ["ACD", "GGGG", "WY"]
    reference = ["FFF", "AAAA", "CC"]

    assert ks_length_statistic(generated, reference) == pytest.approx(0.0)


def test_ks_different_lengths():
    generated = ["AA", "CC"]
    reference = ["AAAAAA", "GGGGGG"]

    assert ks_length_statistic(generated, reference) == pytest.approx(1.0)


def test_generation_metrics():
    generated = ["ACDE", "ACDE", "GGGG", "ACX"]
    reference = ["ACDE", "FFFF", "GGGG"]

    metrics = generation_metrics(generated, reference)

    assert metrics["n_generated"] == 4
    assert metrics["n_valid"] == 3
    assert metrics["validity_rate"] == pytest.approx(0.75)
    assert metrics["uniqueness_rate"] == pytest.approx(0.75)
    assert metrics["length_mean"] == pytest.approx(4.0)

    assert 0 <= metrics["js_divergence"] <= math.log(2)
    assert 0 <= metrics["ks_length"] <= 1


def test_empty_input():
    with pytest.raises(ValueError):
        validity_rate([])

    with pytest.raises(ValueError):
        uniqueness_rate([])


def test_no_valid_peptide():
    with pytest.raises(ValueError):
        sequence_length_stats(["XXX", "123"])

    with pytest.raises(ValueError):
        amino_acid_distribution(["XXX"])
