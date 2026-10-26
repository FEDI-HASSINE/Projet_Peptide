"""Identités d'alignement pour le critère de couverture (A5, AJUSTEMENTS §3/§5).

MODULE DE DÉROGATION : ``identite_couverture`` et ``identite_globale`` sont les
fonctions d'identité prévues par l'ajustement A5, attribuées à Youssef. À
l'instant du fetch, elles n'existaient pas dans ``main`` et leur absence
bloquait toute la chaîne A5/A6 ; ce fragment maison provisoire est donc écrit
par Sofiene sous la règle alias A3 : si Youssef fournit ses fonctions dans
``align.py``, c'est cette API qui est mise en alias (à conserver jusqu'à la fin
de semaine), pas le contenu.

DÉFINITION IMPOSÉE (reprendue verbatim de la définition de mission) :
  - alignement local par Smith-Waterman, gaps affines (o = 11, e = 1),
    substitutions BLOSUM62 (conventions DECISIONS, reprise telles quelles) ;
  - pour l'alignement local de score max :
      identité  = colonnes identiques / colonnes de l'alignement local
      couverture = largeur du fragment / min(len(a), len(b))
    où la largeur du fragment est le plus grand nombre de résidus réels
    (non gaps) alignés de l'une des deux séquences ;
  - ``identite_globale`` : identité colonne par colonne sur l'alignement
    global NW, mêmes paramètres (mesure complémentaire, cohérence §5).

Ancre du §5 : un fragment exact de 20 résidus d'un peptide de 50 a une
identité globale de 0,4 (le test s'échoue sous le seuil) mais une couverture
de 1,0 : la fuite n'est détectable que par le couple (identité, couverture).

Le couple (identité, couverture) >= (0,90 ; 0,90) est le critère de liaison
des clusters v0.2 (``peplib/cluster_v02.py``).
"""

from peplib.align import (
    G_EXTEND,
    G_OPEN,
    alignment_identity,
    needleman_wunsch,
    smith_waterman,
)


def identite_couverture(
    a: str, b: str, gap_open: int = G_OPEN, gap_extend: int = G_EXTEND
) -> tuple:
    """Retourne (identité de l'alignement local max, couverture du fragment).

    L'alignement local optimal est calculé une seule fois ; les deux nombres
    portent sur le même alignement. Si aucune colonne positive n'existe
    (meilleur alignement local vide, plancher SW = 0), le couple est (0,0 ; 0,0).
    Une séquence vide lève ValueError (division par zéro sinon).
    """
    if not a or not b:
        raise ValueError("identite_couverture exige des séquences non vides")

    al = smith_waterman(a, b, gap_open, gap_extend)
    colonnes = len(al.aligned_x)
    if colonnes == 0:                       # alignement local vide (plancher 0)
        return 0.0, 0.0

    identite = alignment_identity(al.aligned_x, al.aligned_y)
    couvert_a = al.x_end - al.x_start
    couvert_b = al.y_end - al.y_start
    couverture = max(couvert_a, couvert_b) / min(len(a), len(b))
    return identite, couverture


def identite_globale(
    a: str, b: str, gap_open: int = G_OPEN, gap_extend: int = G_EXTEND
) -> float:
    """Identité colonne par colonne sur l'alignement global NW.

    Cette mesure est le complément de la couverture, pas son remplacement :
    identique à 1 - Levenshtein/max(len) en esprit, mais avec substitutions
    pondérées par BLOSUM62 et gaps affines (docstring d'align.py).
    """
    if not a or not b:                       # alignement vide : identité 0
        return 0.0
    al = needleman_wunsch(a, b, gap_open, gap_extend)
    return alignment_identity(al.aligned_x, al.aligned_y)
