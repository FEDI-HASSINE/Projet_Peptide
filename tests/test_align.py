"""T36 - Needleman-Wunsch / Smith-Waterman (Gotoh) : ancres de la spec et calculs a la main.
Convention de gap : L -> G_OPEN + (L-1)*G_EXTEND, avec 11 et 1."""
import random

import pytest

from peplib.align import (
    alignment_identity,
    needleman_wunsch,
    score_alignment,
    smith_waterman,
)
from peplib.constants import AA_LIST

NW = needleman_wunsch
SW = smith_waterman


# ---------- NW : identites (ancres de la spec) ----------
def test_nw_identity_AAAA():
    r = NW("AAAA", "AAAA")
    assert r.score == 16 and r.aligned_x == r.aligned_y == "AAAA"


def test_nw_identity_GGGG():
    assert NW("GGGG", "GGGG").score == 24


# ---------- NW : gap affine (ancres de la spec) ----------
def test_nw_one_gap_costs_open():
    # 4 x 4 = 16, un gap de longueur 1 coute 11  ->  5
    r = NW("AAAA", "AAAAA")
    assert r.score == 5
    assert r.aligned_x.count("-") == 1 and len(r.aligned_x) == 5
    assert r.aligned_x.replace("-", "") == "AAAA"


def test_nw_gap_of_three_costs_open_plus_two_extends():
    # 16 - (11 + 2*1) = 3 : trois '-' d'un coup coutent 13, pas 33
    r = NW("AAAA", "AAAAAAA")
    assert r.score == 3
    assert r.aligned_x.count("-") == 3


# ---------- NW : calculs a la main ----------
def test_nw_substitution_by_hand():
    # A/A 4 + C/E -4 + D/D 6 = 6, sans gap
    r = NW("ACD", "AED")
    assert r.score == 6 and r.aligned_x == "ACD" and r.aligned_y == "AED"


def test_nw_internal_gap_by_hand():
    # W/W 11 + gap(11) + W/W 11 = 11  (autres choix : -3)
    r = NW("WKW", "WW")
    assert r.score == 11
    assert r.aligned_x == "WKW" and r.aligned_y == "W-W"


def test_nw_empty_sequences():
    assert NW("", "").score == 0
    assert NW("", "AAA").score == -(11 + 2)   # un gap de longueur 3
    assert NW("AA", "").score == -(11 + 1)


def test_nw_rejects_non_canonical_letters():
    with pytest.raises(ValueError):
        NW("AAXA", "AAAA")


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


# ---------- Identite par alignement ----------
def test_alignment_identity_counts_gaps_in_length():
    assert alignment_identity("ACD-", "ACDE") == pytest.approx(0.75)
    assert alignment_identity("AAAA", "AAAA") == 1.0
    assert alignment_identity("", "") == 0.0


# ---------- Proprietes sur des paires aleatoires (seed 42) ----------
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
    # symetrie du score
    assert g.score == g_rev.score
    # le score annonce = score recalcule sur l'alignement rendu
    assert score_alignment(g.aligned_x, g.aligned_y) == g.score
    assert score_alignment(loc.aligned_x, loc.aligned_y) == loc.score
    # retirer les '-' redonne les sequences
    assert g.aligned_x.replace("-", "") == x and g.aligned_y.replace("-", "") == y
    # local : >= 0, >= global, et morceaux contigus
    assert loc.score >= 0 and loc.score >= g.score
    assert loc.aligned_x.replace("-", "") == x[loc.x_start:loc.x_end]
    assert loc.aligned_y.replace("-", "") == y[loc.y_start:loc.y_end]