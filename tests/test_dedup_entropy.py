"""T24 - Tests redondance (3-mers, Levenshtein) + entropie de composition.

Imports adaptés aux NOMS REELS de peplib (la version initiale de ce
fichier supposait des noms différents : levenshtein, identity,
composition, composition_entropy — ils n'existent pas).

Base des entropies : NATS (np.log), cohérente avec JS calculée en nats ;
les bases de log ne sont PAS unifiées dans le projet (règle spec).
"""
import math

import pytest

dedup = pytest.importorskip("peplib.dedup")
stats = pytest.importorskip("peplib.stats")

levenshtein = dedup.levenshtein_distance
identity = dedup.sequence_identity
find_candidate_pairs = dedup.find_candidate_pairs
composition = stats.compute_composition
composition_entropy = stats.compute_entropy


def _pairs(seqs):
    """Paires candidates sous forme d'ensemble de frozenset({i, j})."""
    return {frozenset(p) for p in find_candidate_pairs(seqs, k=3, min_shared=1)}


# ---------- Levenshtein / identite ----------
def test_levenshtein_kitten_sitting():
    # Ancre spec ; aussi présente dans tests/test_dedup.py (redondance
    # sans danger entre fichiers de test).
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
    # find_candidate_pairs prend les SEQUENCES (pas l'index) avec
    # min_shared=1 : un seul 3-mer partagé suffit pour devenir candidat.
    seqs = ["AAAAAACDEF", "CDEFGGGGGG"]  # partagent CDE
    assert frozenset({0, 1}) in _pairs(seqs)


# ---------- Entropie de composition H = -sum f log f ----------
def test_composition_sums_to_one():
    f = composition(["ACDEFGHIKL", "MNPQRSTVWY", "AAAAAAAAAA"])
    assert sum(f.values()) == pytest.approx(1.0, abs=1e-12)


def test_entropy_single_amino_acid_is_zero():
    # compute_entropy prend le dict de composition.
    assert composition_entropy(composition(["AAAAAAAAAA"])) == pytest.approx(0.0, abs=1e-12)


def test_entropy_uniform_20_aa_is_LN20():
    # H uniforme sur 20 AA = ln(20) nats (PAS log2 : base nats).
    seq = "ACDEFGHIKLMNPQRSTVWY"
    assert composition_entropy(composition([seq])) == pytest.approx(math.log(20), abs=1e-12)


def test_entropy_two_aa_equal_is_LN2():
    # 50/50 sur 2 lettres = ln(2) nats.
    assert composition_entropy(composition(["AAAACCCC"])) == pytest.approx(math.log(2), abs=1e-12)
