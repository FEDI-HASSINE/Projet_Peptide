"""T19 - Tests des descripteurs sur des peptides calculables a la main (spec section 6)."""
import math

import pytest

from peplib.constants import EISENBERG, N_TERMINUS_PKA, C_TERMINUS_PKA
from peplib.descriptors import (
    calculate_charge,
    calculate_gravy,
    calculate_hydrophobic_moment,
    calculate_molecular_mass,
    calculate_pi,
)


# ---------- GRAVY ----------
@pytest.mark.parametrize("seq, expected", [("A", 1.8), ("GG", -0.4), ("AG", 0.7)])
def test_gravy_hand_values(seq, expected):
    assert calculate_gravy(seq) == pytest.approx(expected, abs=1e-12)


# ---------- Masse monoisotopique ----------
def test_mass_glycine():
    assert calculate_molecular_mass("G") == pytest.approx(75.0320, abs=1e-4)


def test_mass_aa():
    assert calculate_molecular_mass("AA") == pytest.approx(160.0848, abs=1e-4)


# ---------- Charge ----------
def test_charge_sign_at_ph7():
    assert calculate_charge("DE", 7.0) < 0
    assert calculate_charge("KK", 7.0) > 0


def test_charge_glycine_ph1_by_hand():
    # Formule ecrite a la main, independante du code de descriptors.py
    expected = 1 / (1 + 10 ** (1 - N_TERMINUS_PKA)) - 1 / (1 + 10 ** (C_TERMINUS_PKA - 1))
    assert calculate_charge("G", 1.0) == pytest.approx(expected, abs=1e-3)


# ---------- pI ----------
@pytest.mark.parametrize("seq", ["G", "KKD", "DEDE", "KRKR", "ACDEFGHIKLMNPQRSTVWY"])
def test_charge_at_pi_is_zero(seq):
    assert calculate_charge(seq, calculate_pi(seq)) == pytest.approx(0.0, abs=1e-3)


@pytest.mark.parametrize("seq", ["G", "KKD", "ACDEFGHIKLMNPQRSTVWY"])
def test_charge_strictly_decreasing_in_ph(seq):
    phs = [i * 0.5 for i in range(0, 29)]  # 0 -> 14
    charges = [calculate_charge(seq, p) for p in phs]
    assert all(a > b for a, b in zip(charges, charges[1:]))


# ---------- Moment hydrophobe ----------
@pytest.mark.parametrize("aa", ["A", "R", "G", "L"])
def test_hydrophobic_moment_single_residue(aa):
    assert calculate_hydrophobic_moment(aa) == pytest.approx(abs(EISENBERG[aa]))


def test_hydrophobic_moment_k0_is_abs_sum_over_n():
    # angle = 0 : tous les vecteurs sont alignes -> |somme H| (divisee par N dans notre code)
    seq = "ARL"
    expected = abs(sum(EISENBERG[a] for a in seq)) / len(seq)
    assert calculate_hydrophobic_moment(seq, angle=0.0) == pytest.approx(expected)


def test_hydrophobic_moment_is_never_negative_or_nan():
    for seq in ["AKLLKKLLKL", "GGGGGG", "DEDEDE"]:
        mu = calculate_hydrophobic_moment(seq)
        assert mu >= 0 and math.isfinite(mu)