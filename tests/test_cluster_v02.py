"""Tests du clustering v0.2 au critère de couverture (A5)."""

import math
from pathlib import Path

import pytest

from peplib.cluster_v02 import (
    candidats_par_couverture,
    borne_intacts,
    cluster_v02,
    tableau_t13,
)
from peplib.identite import identite_couverture

FRAGMENT = "KTIDWREEKGVVDAGYSQLC"             # ancre identitaire (20 résidus)
PEPTIDE_50 = "ADPQSIRNLDEVWAK" + FRAGMENT + "MHTFLKGQESRVDCP"


def test_borne_intacts_valeurs_exactes():
    """Valeurs issues de la preuve A6 : intacts >= 0,63·m − 3, arrondi bas."""
    assert borne_intacts(10) == 3       # ⌊3,3⌋
    assert borne_intacts(11) == 3       # ⌊3,93⌋
    assert borne_intacts(20) == 9       # ⌊9,6⌋
    assert borne_intacts(50) == 28      # ⌊28,5⌋
    with pytest.raises(ValueError):
        borne_intacts(10, tau=1.5)


def test_le_fragment_est_candidat():
    """Un fragment exact partage tous les 3-mers du fragment : candidat交融."""
    seqs = [PEPTIDE_50, FRAGMENT]
    stats = {}
    candidats = candidats_par_couverture(seqs, tau=0.90, stats=stats)
    assert (0, 1) in candidats
    # borne invoquée : intactes >= ⌊0,63·20 − 3⌋ = 9 ; le fragment en fournit 18
    assert stats["candidats_couverture"] >= 1


def test_paires_independantes_sont_coupees():
    """Deux séquences sans aucun 3-mer commun ne survivent PAS au préfiltre."""
    seqs = ["AAAAAAAAAAAA", "CGCGCGCGCGCG"]
    candidats = candidats_par_couverture(seqs, tau=0.90)
    assert candidats == set()


def test_fragment_rejoint_le_peptide_et_pas_les_autres():
    """Lien par couverture : le fragment et son peptide long forment une famille ;
    une séquence indépendante reste seule."""
    independent = "MMMMMMMMMMMMMMMMMMMMMM"     # homopes sans aucun 3-mer fragmentaire...
    # attention : vérifions qu'il n'a aucun 3-mer commun avec FRAGMENT
    assert not ({FRAGMENT[i : i + 3] for i in range(len(FRAGMENT) - 2)}
                & {independent[i : i + 3] for i in range(len(independent) - 2)})
    seqs = [PEPTIDE_50, FRAGMENT, independent]
    etiquettes = cluster_v02(seqs, tau=0.90)
    fragmente = etiquettes[0] == etiquettes[1]
    assert fragmente
    assert etiquettes[2] != etiquettes[0]


def test_identiques_sont_dans_la_meme_famille():
    a = FRAGMENT
    b = FRAGMENT + "AAAA"
    seqs = [a, b, "APWQ" * 4 + "T"]   # autre séquence débranchée
    etiquettes = cluster_v02(seqs)
    assert etiquettes[0] == etiquettes[1]
    assert etiquettes[2] != etiquettes[0]


def test_cluster_determinisme():
    seqs = [PEPTIDE_50, FRAGMENT, "RPNQWVDEEGILKHFSTACM", "TTTTTTTTTT"]
    e1 = cluster_v02(seqs)
    e2 = cluster_v02(seqs)
    assert e1 == e2


def test_tableau_t13(tmp_path):
    sortie = Path(tmp_path) / "t13.csv"
    # 5 etiquettes : 1 famille de 3, 1 de 2 (pas de singletons ici)
    df = tableau_t13([0, 0, 0, 1, 1], sortie=sortie)
    assert df["nb_familles"][0] == 2
    assert df["taille_max"][0] == 3
    assert df["singletons"][0] == 0
    assert df["fraction_plus_grande"][0] == 0.6
    # lecture : le fichier existe à la sortie donnée
    assert sortie.exists()


def test_criterer_couverture_seuil_0_9_exactement():
    """Un alignement avec exactement 90 % d'identité et de couverture est assimilé.

    18 colonnes identiques sur 20 exactement (18/20 = 0,9) : la comparaison
    flottante est exacte pour W ≤ 200 (cf. docstring du module), vérifier ici.
    """
    a = "ACDEFGHIKLNPQRSTVWYA"      # 20
    # b : 18 résidus alignés exactement + 2 substitutions au score BLOSUM62
    # positif (Y-F = +3, A-S = +1) pour que SW aligne les 20 colonnes.
    b = a[:18] + "FS"               # colonne identité 18/20 = 0,9, couverture 20/20
    identite, couverture = identite_couverture(a, b)
    assert couverture == 1.0
    assert identite == pytest.approx(0.9)
    etiquettes = cluster_v02([a, b])
    assert etiquettes[0] == etiquettes[1]
