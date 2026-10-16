import pytest

from peplib.constants import (
    AMINO_ACIDS,
    MONOISOTOPIC_RESIDUE_MASS,
    SIDE_CHAIN_PKA,
    N_TERMINUS_PKA,
    C_TERMINUS_PKA,
    KYTE_DOOLITTLE,
    EISENBERG,
)


def test_twenty_standard_amino_acids():
    assert len(AMINO_ACIDS) == 20
    assert len(set(AMINO_ACIDS)) == 20


@pytest.mark.parametrize(
    "table",
    [
        MONOISOTOPIC_RESIDUE_MASS,
        KYTE_DOOLITTLE,
        EISENBERG,
    ],
)
def test_complete_amino_acid_tables(table):
    assert set(table) == set(AMINO_ACIDS)


def test_side_chain_pka_contains_ionizable_residues():
    assert set(SIDE_CHAIN_PKA) == set("CDEHKRY")


def test_terminal_pka_values():
    assert N_TERMINUS_PKA == 9.69
    assert C_TERMINUS_PKA == 2.34


def test_known_masses():
    assert MONOISOTOPIC_RESIDUE_MASS["A"] == pytest.approx(71.03711)
    assert MONOISOTOPIC_RESIDUE_MASS["G"] == pytest.approx(57.02146)


def test_known_hydrophobicity_values():
    assert KYTE_DOOLITTLE["A"] == pytest.approx(1.8)
    assert KYTE_DOOLITTLE["G"] == pytest.approx(-0.4)
    assert EISENBERG["A"] == pytest.approx(0.62)


def test_required_side_chain_pka_values():
    assert SIDE_CHAIN_PKA["D"] == pytest.approx(3.9)
    assert SIDE_CHAIN_PKA["E"] == pytest.approx(4.07)
    assert SIDE_CHAIN_PKA["C"] == pytest.approx(8.18)
    assert SIDE_CHAIN_PKA["Y"] == pytest.approx(10.46)
    assert SIDE_CHAIN_PKA["H"] == pytest.approx(6.04)
    assert SIDE_CHAIN_PKA["K"] == pytest.approx(10.54)
    assert SIDE_CHAIN_PKA["R"] == pytest.approx(12.48)
