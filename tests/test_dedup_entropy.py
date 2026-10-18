"""T24 - Tests redondance (3-mers, Levenshtein) + entropie de composition.

Ces tests se SAUTENT tant que dedup.py (Sofiene, T23) et stats.py (Yassine, T16)
ne sont pas pushes. Si les noms de fonctions different, adapter UNIQUEMENT les imports.
Noms supposes : levenshtein, identity, build_kmer_index, candidate_pairs (dedup.py)
                composition, composition_entropy (stats.py)
"""
import math

import pytest

dedup = pytest.importorskip("peplib.dedup")
stats = pytest.importorskip("peplib.stats")

levenshtein = dedup.levenshtein
identity = dedup.identity
build_kmer_index = dedup.build_kmer_index
candidate_pairs = dedup.candidate_pairs
composition = stats.composition
composition_entropy = stats.composition_entropy


def _pairs(seqs):
    """Paires candidates sous forme d'ensemble de frozenset({i, j})."""
    return {frozenset(p) for p in candidate_pairs(build_kmer_index(seqs, 3))}


# ---------- Levenshtein / identite ----------
def test_levenshtein_kitten_sitting():
    assert levenshtein("kitten", "sitting") == 3


def test_levenshtein_basic():
    assert levenshtein("", "abc") == 3
    assert levenshtein("abc", "abc") == 0
    assert levenshtein("abc", "abd") == 1


def test_identity_formula():
    # 1 - d / max(|a|,|b|)
    assert identity("AAAAAAAAAA", "AAAAAAAAAA") == pytest.approx(1.0)
    assert identity("AAAAAAAAAA", "AAAAAAAAAC") == pytest.approx(0.9)


# ---------- Index de 3-mers ----------
def test_identical_sequences_are_candidates():
    seqs = ["ACDEFGHIKL", "ACDEFGHIKL"]
    assert frozenset({0, 1}) in _pairs(seqs)


def test_sequences_without_common_3mer_are_not_candidates():
    seqs = ["AAAAAAAAAA", "CCCCCCCCCC"]
    assert frozenset({0, 1}) not in _pairs(seqs)


def test_one_shared_3mer_is_enough():
    seqs = ["AAAAAACDEF", "CDEFGGGGGG"]  # partagent CDE
    assert frozenset({0, 1}) in _pairs(seqs)


# ---------- Entropie de composition H = -sum f log f ----------
def test_composition_sums_to_one():
    f = composition(["ACDEFGHIKL", "MNPQRSTVWY", "AAAAAAAAAA"])
    assert sum(f.values()) == pytest.approx(1.0, abs=1e-12)


def test_entropy_single_amino_acid_is_zero():
    assert composition_entropy(["AAAAAAAAAA"]) == pytest.approx(0.0, abs=1e-12)


def test_entropy_uniform_20_aa_is_log20():
    seq = "ACDEFGHIKLMNPQRSTVWY"
    assert composition_entropy([seq]) == pytest.approx(math.log(20), abs=1e-12)


def test_entropy_two_aa_equal_is_log2():
    assert composition_entropy(["AAAACCCC"]) == pytest.approx(math.log(2), abs=1e-12)