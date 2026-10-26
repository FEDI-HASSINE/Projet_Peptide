"""T36 - Needleman-Wunsch / Smith-Waterman (Gotoh) : ancres de la spec et calculs à la main.

Convention de gap : L -> G_OPEN + (L-1)*G_EXTEND, avec 11 et 1.
"""
import random

import pytest

from peplib.align import (
    alignment_identity,
    needleman_wunsch,
    score_alignment,
    smith_waterman,
)
from peplib.constants import AA_LIST, blosum62_score

NW = needleman_wunsch
SW = smith_waterman


# ---------- NW : identités (ancres de la spec) ----------
def test_nw_identity_AAAA():
    r = NW("AAAA", "AAAA")
    assert r.score == 16 and r.aligned_x == r.aligned_y == "AAAA"


def test_nw_identity_GGGG():
    assert NW("GGGG", "GGGG").score == 24


# ---------- NW : gap affine (ancres de la spec) ----------
def test_nw_one_gap_costs_open():
    # 4 x 4 = 16, un gap de longueur 1 coûte 11  ->  5
    r = NW("AAAA", "AAAAA")
    assert r.score == 5
    assert r.aligned_x.count("-") == 1 and len(r.aligned_x) == 5
    assert r.aligned_x.replace("-", "") == "AAAA"


def test_nw_gap_of_three_costs_open_plus_two_extends():
    # 16 - (11 + 2*1) = 3 : trois '-' d'un coup coûtent 13, pas 33
    r = NW("AAAA", "AAAAAAA")
    assert r.score == 3
    assert r.aligned_x.count("-") == 3


# ---------- NW : calculs à la main ----------
def test_nw_substitution_by_hand():
    # A/A 4 + C/E -4 + D/D 6 = 6, sans gap
    r = NW("ACD", "AED")
    assert r.score == 6 and r.aligned_x == "ACD" and r.aligned_y == "AED"


def test_nw_internal_gap_by_hand():
    # W/W 11 + gap(11) + W/W 11 = 11  (autres choix : -3)
    r = NW("WKW", "WW")
    assert r.score == 11
    assert r.aligned_x == "WKW" and r.aligned_y == "W-W"


def test_nw_leading_gap_by_hand():
    # W/A vaut -3 : le W de tête est forcément décalé, gap de 1 en tête
    # 4 x A/A = 16, moins 11 -> 5
    r = NW("AAAA", "WAAAA")
    assert r.score == 5
    assert (r.aligned_x, r.aligned_y) == ("-AAAA", "WAAAA")


def test_nw_trailing_gap_by_hand():
    # même calcul, gap en fin : 16 - 11 = 5
    r = NW("AAAA", "AAAAW")
    assert r.score == 5
    assert (r.aligned_x, r.aligned_y) == ("AAAA-", "AAAAW")


def test_nw_two_separated_gaps_by_hand():
    # 6 x W/W = 66, deux gaps de longueur 1 séparés par des appariements
    # 66 - 2*11 = 44 (et non 66 - 22 refacturé ou 66 - 2*33)
    r = NW("WWKWWKWW", "WWWWWW")
    assert r.score == 44
    assert r.aligned_y.count("-") == 2          # x est le plus long : les gaps sont dans y
    assert r.aligned_y == "WW-WW-WW"


def test_nw_empty_sequences():
    assert NW("", "").score == 0
    assert NW("", "AAA").score == -(11 + 2)   # un gap de longueur 3
    assert NW("AA", "").score == -(11 + 1)


# ---------- NW : validation de l'alphabet ----------
def test_nw_rejects_non_canonical_letters():
    with pytest.raises(ValueError):
        NW("AAXA", "AAAA")


@pytest.mark.parametrize("bad", ["X", "B", "Z", "J", "U", "O", "a", "1", "-"])
def test_nw_rejects_bad_letter_in_either_sequence(bad):
    with pytest.raises(ValueError):
        NW("AAA" + bad, "AAAA")
    with pytest.raises(ValueError):
        NW("AAAA", "AAA" + bad)


@pytest.mark.parametrize("bad", ["X", "B", "Z", "J", "U", "O", "a", "1", "-"])
def test_check_and_blosum62_score_agree(bad):
    # une seule règle d'alphabet : _check et blosum62_score refusent les mêmes lettres
    with pytest.raises(ValueError):
        NW("AA" + bad, "AAA")
    with pytest.raises(ValueError):
        blosum62_score(bad, "A")


# ---------- SW ----------
def test_sw_identity_AAAA():
    assert SW("AAAA", "AAAA").score == 16


def test_sw_local_part_is_contiguous_substring():
    r = SW("WWWW", "KKWWWWKK")
    assert r.score == 44                       # 4 x 11
    assert r.aligned_x == r.aligned_y == "WWWW"
    assert "WWWW"[r.x_start:r.x_end] == "WWWW"
    assert "KKWWWWKK"[r.y_start:r.y_end] == "WWWW"
    assert (r.y_start, r.y_end) == (2, 6)


def test_sw_nothing_in_common_gives_zero():
    r = SW("WWWW", "DDDD")                    # W/D = -4
    assert r.score == 0 and r.aligned_x == "" and r.aligned_y == ""


def test_sw_empty_sequence_gives_zero():
    assert SW("", "AAA").score == 0
    assert SW("AAA", "").score == 0


# ---------- score_alignment : score recalculé ----------
def test_score_alignment_adjacent_opposite_gaps_refacture_open():
    # motif impossible en sortie de DP (saisie utilisateur seulement) :
    # un '-' dans x puis un '-' dans y comptent deux ouvertures -> -22
    assert score_alignment("A-", "-C") == -22


def test_score_alignment_rejects_unequal_lengths():
    # sans garde, zip tronque et rend 19 pour trois colonnes comparées
    with pytest.raises(ValueError):
        score_alignment("ACD", "ACDE")


def test_score_alignment_rejects_non_canonical_residue():
    # sans garde, l'indexation directe de BLOSUM62 laisserait fuir un KeyError
    with pytest.raises(ValueError):
        score_alignment("AX", "AA")


# ---------- Identité par alignement ----------
def test_alignment_identity_counts_gaps_in_length():
    assert alignment_identity("ACD-", "ACDE") == pytest.approx(0.75)
    assert alignment_identity("AAAA", "AAAA") == 1.0
    assert alignment_identity("", "") == 0.0


def test_alignment_identity_rejects_unequal_lengths():
    # sans garde, "ACD" / "ACDE" rend 1.0 : trois colonnes comparées, division par 3
    with pytest.raises(ValueError):
        alignment_identity("ACD", "ACDE")


# ---------- Propriétés sur des paires aléatoires (seed 42) ----------
def _random_pairs(n=60):
    rnd = random.Random(42)
    out = []
    for _ in range(n):
        x = "".join(rnd.choice(AA_LIST) for _ in range(rnd.randint(1, 30)))
        y = "".join(rnd.choice(AA_LIST) for _ in range(rnd.randint(1, 30)))
        out.append((x, y))
    return out


@pytest.mark.parametrize("x, y", _random_pairs())
def test_properties_on_random_pairs(x, y):
    g, g_rev = NW(x, y), NW(y, x)
    loc = SW(x, y)
    # symétrie du score
    assert g.score == g_rev.score
    # le score annoncé = score recalculé sur l'alignement rendu
    assert score_alignment(g.aligned_x, g.aligned_y) == g.score
    assert score_alignment(loc.aligned_x, loc.aligned_y) == loc.score
    # retirer les '-' redonne les séquences
    assert g.aligned_x.replace("-", "") == x and g.aligned_y.replace("-", "") == y
    # local : >= 0, >= global, et morceaux contigus
    assert loc.score >= 0 and loc.score >= g.score
    assert loc.aligned_x.replace("-", "") == x[loc.x_start:loc.x_end]
    assert loc.aligned_y.replace("-", "") == y[loc.y_start:loc.y_end]


# ---------- T7 : identite_couverture et identite_globale ----------
from peplib.align import identite_couverture, identite_globale


# --- ancres calculées à la main ---

def test_identite_couverture_perfect_no_gaps():
    # Aucun gap, identique : couverture = 1.0
    assert identite_couverture("AAAA", "AAAA") == pytest.approx(1.0)


def test_identite_globale_perfect_no_gaps():
    assert identite_globale("AAAA", "AAAA") == pytest.approx(1.0)


def test_identite_couverture_gap_in_longer():
    # aligned_x = "AAAA-", aligned_y = "AAAAW"
    # x_sans_gap = 4, y_sans_gap = 5 -> min = 4
    # identiques : A/A A/A A/A A/A (4 match, gap non compté)
    # identite_couverture = 4/4 = 1.0
    r = NW("AAAA", "AAAAW")
    assert r.aligned_x == "AAAA-" and r.aligned_y == "AAAAW"
    assert identite_couverture(r.aligned_x, r.aligned_y) == pytest.approx(1.0)


def test_identite_globale_gap_in_longer():
    # même alignement : max(4, 5) = 5, 4 matches -> 4/5 = 0.8
    r = NW("AAAA", "AAAAW")
    assert identite_globale(r.aligned_x, r.aligned_y) == pytest.approx(0.8)


def test_identite_couverture_internal_substitution():
    # NW("ACD", "AED") -> "ACD" / "AED", score = 6, aucun gap
    # identiques : A/A (1), D/D (1) ; non identique C/E
    # x_len = y_len = 3, min = 3 -> 2/3
    r = NW("ACD", "AED")
    assert r.aligned_x == "ACD" and r.aligned_y == "AED"
    assert identite_couverture(r.aligned_x, r.aligned_y) == pytest.approx(2 / 3)


def test_identite_globale_internal_substitution():
    # max(3, 3) = 3 -> 2/3
    r = NW("ACD", "AED")
    assert identite_globale(r.aligned_x, r.aligned_y) == pytest.approx(2 / 3)


def test_identite_couverture_empty_strings():
    assert identite_couverture("", "") == 0.0


def test_identite_globale_empty_strings():
    assert identite_globale("", "") == 0.0


def test_identite_couverture_rejects_unequal_lengths():
    with pytest.raises(ValueError):
        identite_couverture("ACD", "ACDE")


def test_identite_globale_rejects_unequal_lengths():
    with pytest.raises(ValueError):
        identite_globale("ACD", "ACDE")


def test_identite_all_gaps_in_x():
    # aligned_x = "---", aligned_y = "AAA" : x_sans_gap = 0, denom = min(0,3) = 0
    assert identite_couverture("---", "AAA") == 0.0
    assert identite_globale("---", "AAA") == 0.0


# --- invariants sur paires aléatoires ---

@pytest.mark.parametrize("x, y", _random_pairs())
def test_identite_couverture_globale_invariants(x, y):
    g = NW(x, y)
    cov = identite_couverture(g.aligned_x, g.aligned_y)
    glob = identite_globale(g.aligned_x, g.aligned_y)
    ai = alignment_identity(g.aligned_x, g.aligned_y)

    # bornes
    assert 0.0 <= glob <= cov <= 1.0
    # couverture >= alignment_identity
    assert cov >= ai - 1e-12
    # symmetrie : NW(y, x) donne le meme score d'identite
    g_rev = NW(y, x)
    cov_rev = identite_couverture(g_rev.aligned_x, g_rev.aligned_y)
    glob_rev = identite_globale(g_rev.aligned_x, g_rev.aligned_y)
    assert cov == pytest.approx(cov_rev, abs=1e-12)
    assert glob == pytest.approx(glob_rev, abs=1e-12)


@pytest.mark.parametrize("x, y", _random_pairs())
def test_identite_couverture_is_one_for_identical_sequences(x, y):
    # Une sequence avec elle-meme -> couverture 1.0, globale 1.0
    g = NW(x, x)
    assert identite_couverture(g.aligned_x, g.aligned_y) == pytest.approx(1.0)
    assert identite_globale(g.aligned_x, g.aligned_y) == pytest.approx(1.0)
