"""T33 - BLOSUM62 : saisie, symetrie, diagonale, ancres de la source citee."""
import pytest

from peplib.constants import AA_LIST, BLOSUM62, blosum62_score


def test_400_values():
    assert set(BLOSUM62) == set(AA_LIST)
    assert all(set(row) == set(AA_LIST) for row in BLOSUM62.values())
    assert sum(len(row) for row in BLOSUM62.values()) == 400


def test_values_are_integers():
    assert all(isinstance(v, int) for row in BLOSUM62.values() for v in row.values())


def test_symmetric():
    for a in AA_LIST:
        for b in AA_LIST:
            assert BLOSUM62[a][b] == BLOSUM62[b][a], (a, b)


def test_diagonal_positive_and_dominates_row():
    for a in AA_LIST:
        assert BLOSUM62[a][a] > 0
        assert BLOSUM62[a][a] == max(BLOSUM62[a].values()), a


@pytest.mark.parametrize("a, b, expected", [
    ("A", "A", 4), ("W", "W", 11), ("C", "C", 9), ("G", "G", 6),
    ("L", "I", 2), ("I", "V", 3), ("D", "E", 2), ("K", "R", 2),
    ("F", "Y", 3), ("H", "Y", 2),
    ("A", "W", -3), ("W", "G", -2), ("D", "L", -4), ("C", "E", -4),
])
def test_anchor_values(a, b, expected):
    assert blosum62_score(a, b) == expected
    assert blosum62_score(b, a) == expected


def test_conservative_cheaper_than_radical():
    assert blosum62_score("L", "I") > blosum62_score("A", "W")


def test_matrix_sum_is_known_checksum():
    # somme des 400 valeurs (controle de saisie global, verifiee contre le fichier NCBI)
    assert sum(v for row in BLOSUM62.values() for v in row.values()) == -426


@pytest.mark.parametrize("residu", ["X", "B", "Z", "J", "U", "O", "*", "-"])
def test_ambiguous_residues_rejected(residu):
    # contrat docstring : toute valeur hors des 20 AA canoniques leve
    # une ValueError explicite (matrice sans ambigus, rejet T12/T13)
    with pytest.raises(ValueError):
        blosum62_score(residu, "A")
    with pytest.raises(ValueError):
        blosum62_score("A", residu)
