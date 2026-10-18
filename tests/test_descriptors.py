import pytest

from peplib.descriptors import calculate_charge
from peplib.descriptors import calculate_pi
from peplib.descriptors import calculate_gravy
from peplib.descriptors import calculate_hydrophobic_moment
from peplib.descriptors import calculate_boman_index
from peplib.descriptors import calculate_molecular_mass

def test_charge_lysine_at_ph7():
    assert calculate_charge("K", 7) == pytest.approx(
        0.9976959794, abs=1e-8
    )


def test_charge_kkd_at_ph7():
    assert calculate_charge("KKD", 7) == pytest.approx(
        0.9982013572, abs=1e-8
    )


def test_charge_glycine_at_ph1():
    assert calculate_charge("G", 1) == pytest.approx(
        0.9562891502, abs=1e-8
    )


def test_sequence_is_case_insensitive():
    assert calculate_charge("kkd", 7) == pytest.approx(
        calculate_charge("KKD", 7)
    )


def test_invalid_sequence():
    with pytest.raises(ValueError):
        calculate_charge("KXZ", 7)


def test_empty_sequence():
    with pytest.raises(ValueError):
        calculate_charge("", 7)

def test_pi_is_between_zero_and_fourteen():
    pi = calculate_pi("KKD")
    assert 0 <= pi <= 14


def test_charge_at_pi_is_close_to_zero():
    pi = calculate_pi("KKD")
    assert abs(calculate_charge("KKD", pi)) < 1e-5


def test_pi_is_case_insensitive():
    assert calculate_pi("kkd") == pytest.approx(
        calculate_pi("KKD")
    )


def test_pi_rejects_empty_sequence():
    with pytest.raises(ValueError):
        calculate_pi("")

def test_gravy_gk():
    assert calculate_gravy("GK") == pytest.approx(-2.15)


def test_gravy_is_case_insensitive():
    assert calculate_gravy("gk") == pytest.approx(
        calculate_gravy("GK")
    )


def test_hydrophobic_moment_single_amino_acid():
    assert calculate_hydrophobic_moment("A") == pytest.approx(0.62)


def test_hydrophobic_moment_is_case_insensitive():
    assert calculate_hydrophobic_moment("gk") == pytest.approx(
        calculate_hydrophobic_moment("GK")
    )


def test_hydrophobic_moment_rejects_empty_sequence():
    with pytest.raises(ValueError):
        calculate_hydrophobic_moment("")


def test_boman_single_amino_acid():
    assert calculate_boman_index("R") == pytest.approx(14.92)


def test_boman_average():
    # A = -1.81 et R = 14.92
    assert calculate_boman_index("AR") == pytest.approx(6.555)


def test_boman_is_case_insensitive():
    assert calculate_boman_index("ar") == pytest.approx(
        calculate_boman_index("AR")
    )


def test_boman_rejects_empty_sequence():
    with pytest.raises(ValueError):
        calculate_boman_index("")




def test_molecular_mass_alanine():
    assert calculate_molecular_mass("A") == pytest.approx(89.04767)


def test_molecular_mass_two_residues():
    # Glycine + alanine + une molécule d'eau
    expected = 57.02146 + 71.03711 + 18.01056
    assert calculate_molecular_mass("GA") == pytest.approx(expected)


def test_molecular_mass_is_case_insensitive():
    assert calculate_molecular_mass("ga") == pytest.approx(
        calculate_molecular_mass("GA")
    )


def test_molecular_mass_rejects_empty_sequence():
    with pytest.raises(ValueError):
        calculate_molecular_mass("")
