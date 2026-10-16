from peplib.constants import AA_LIST, MASS, KYTE_DOOLITTLE, EISENBERG


def test_twenty_amino_acids():
    assert len(AA_LIST) == 20
    assert set(AA_LIST) == set("ACDEFGHIKLMNPQRSTVWY")


def test_masses_positive_and_complete():
    assert len(MASS) == 20
    assert all(m > 0 for m in MASS.values())


def test_scales_have_20_values():
    assert len(KYTE_DOOLITTLE) == 20
    assert len(EISENBERG) == 20


def test_known_values():
    assert MASS["G"] == 57.02146
    assert MASS["A"] == 71.03711
    assert KYTE_DOOLITTLE["A"] == 1.8
    assert KYTE_DOOLITTLE["G"] == -0.4