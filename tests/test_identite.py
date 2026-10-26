"""Tests des identités de couverture et globale (A5, fragment provisoire)."""

import pytest

from peplib.identite import identite_couverture, identite_globale

# Fragment exact de 20 résidus d'un peptide de 50 : ancre du §5 du doc Ajustements
FRAGMENT = "KTIDWREEKGVVDAGYSQLC"          # 20 résidus canoniques
# Préfixe et suffixe de 15 résidus canoniques quoi qu'il arrive
PEPTIDE_50 = "ADPQSIRNLDEVWAK" + FRAGMENT + "MHTFLKGQESRVDCP"   # 15 + 20 + 15

assert len(PEPTIDE_50) == 50


def test_veille_fragment_exact_d_un_peptide_de_50():
    """Ancre §5 : couverture 1,0 pour un fragment exact, alors que l'identité
    globale vaut 20/50 = 0,4 (sous le seuil 0,90 ; la fuite passe le global)."""
    identite, couverture = identite_couverture(PEPTIDE_50, FRAGMENT)

    assert identite == 1.0
    assert couverture == 1.0
    # Le global ne détecte pas la fuite : c'est l'objet de toute la mise en garde.
    # 20 colonnes identiques / 50 colonnes globales (l'alignement global aligne
    # les 30 résidus complémentaires du 50-mer à des gaps).
    assert identite_globale(PEPTIDE_50, FRAGMENT) == pytest.approx(0.4)


def test_identite_globale_s_alignes_parfaitement():
    assert identite_globale(FRAGMENT, FRAGMENT) == 1.0


def test_identite_couverture_s_alignes_parfaitement():
    identite, couverture = identite_couverture(FRAGMENT, FRAGMENT)
    assert identite == 1.0
    assert couverture == 1.0


def test_alignement_local_vide_rend_zeros():
    """Paire dont toute colonne est négative : SW rend le plancher 0 (alignement
    local vide), le couple est (0,0 ; 0,0) sans division par zéro."""
    # BLOSUM62('G','C') = -3 ; le couple le plus court ne peut s'aligner
    identite, couverture = identite_couverture("G", "C")
    assert (identite, couverture) == (0.0, 0.0)


def test_sequence_vide_leve():
    with pytest.raises(ValueError):
        identite_couverture("", "KTIDWREEKGVVDAGYSQLC")
    with pytest.raises(ValueError):
        identite_couverture("KTIDWREEKGVVDAGYSQLC", "")


def test_couple_borne_dans_0_1():
    """Couple (identité, couverture) borné sur une paire partiellement homologue."""
    identite, couverture = identite_couverture(FRAGMENT, FRAGMENT + "AAAA")
    assert 0.0 <= identite <= 1.0
    assert 0.0 <= couverture <= 1.0
    # le fragment de FRAGMENT est aligné pur : couverture 1,0, identite >= 0,8
    assert couverture == 1.0
    assert identite >= 0.8
