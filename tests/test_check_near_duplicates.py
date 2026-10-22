"""Tests du module de vérification indépendante des quasi-doublons (T25).

Anchored to the Spec §6 : Levenshtein("kitten", "sitting") == 3 ;
frontière d'identité exactement 0.90 = quasi-doublon compté ;
0.80 hors seuil. Le bandé est croisé avec l'implémentation de référence
peplib.dedup.levenshtein_distance sur 100 paires aléatoires (seed 42).
"""

import random

import pytest

from peplib.check_near_duplicates import levenshtein_within, near_duplicate_stats
from peplib.dedup import levenshtein_distance


def test_levenshtein_within_ancres():
    """Ancre obligatoire du spec : distance(kitten, sitting) = 3."""
    assert levenshtein_distance("kitten", "sitting") == 3
    assert levenshtein_within("kitten", "sitting", 3) is True
    assert levenshtein_within("kitten", "sitting", 2) is False


def test_bande_dmax_zero():
    """Un budget 0 ne tolère qu'une séquence strictement identique."""
    assert levenshtein_within("ACTGTGAGTC", "ACTGTGAGTC", 0) is True
    assert levenshtein_within("ACTGTGAGTC", "ACTGTGAGTA", 0) is False
    assert levenshtein_within("ACTGTGAGTC", "ACTGTGAGTCA", 0) is False


def test_frontiere_identite_090():
    "Cas limite du spec : 1 édition sur une longueur 10 -> 0,90 pile, compté."
    a = "ACDEFGHIKL"
    b = "ACDEFGHIKM"  # 1 substitution
    st = near_duplicate_stats([a, b])
    assert st["unique_sequences"] == 2
    assert st["near_duplicate_pairs"] == 1
    assert st["sequences_in_near_dups"] == 2
    assert st["near_dup_rate"] == pytest.approx(1.0)


def test_sous_seuil_exclu():
    "2 éditions sur une longueur 10 -> 0,80 < 0,90, exclu."
    st = near_duplicate_stats(["ACDEFGHIKL", "ACDEAGHIKM"])  # 2 substitutions
    assert st["near_duplicate_pairs"] == 0
    assert st["sequences_in_near_dups"] == 0
    assert st["near_dup_rate"] == 0.0


def test_doublon_exact_reduction():
    """Les doublons exacts sont comptés une seule fois AVANT le quasi-test."""
    st = near_duplicate_stats(["ACDEFGHIKL", "ACDEFGHIKL", "ACDEFGHIKM"])
    assert st["unique_sequences"] == 2
    assert st["near_duplicate_pairs"] == 1


def test_liste_vide():
    st = near_duplicate_stats([])
    assert st["unique_sequences"] == 0
    assert st["near_duplicate_pairs"] == 0
    # "if n else 0.0" protege la division par zero
    assert st["near_dup_rate"] == 0.0


def test_bande_equivalent_leven_complet():
    """Propriété : le bandé ne renvoie jamais True si la distance > d_max
    ni False si la distance <= d_max ; équivalence stricte avec le Levenshtein
    de référence peplib.dedup sur 100 paires aléatoires (seed 42)."""
    rng = random.Random(42)
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    for _ in range(100):
        n = rng.randint(5, 15)
        a = "".join(rng.choice(alphabet) for _ in range(n))
        b_list = list(a)
        for p in range(len(b_list)):  # chaque position : garder / muter
            r = rng.random()
            if r < 0.15:
                b_list[p] = rng.choice(alphabet)  # substitution
            elif r < 0.22:
                b_list[p] = ""                   # délétion
        b = "".join(b_list)
        if rng.random() < 0.3:                   # insertion possible
            pos = rng.randint(0, len(b))
            b = b[:pos] + rng.choice(alphabet) + b[pos:]
        d_reel = levenshtein_distance(a, b)
        budget = rng.randint(0, 8)
        assert levenshtein_within(a, b, budget) == (d_reel <= budget)
