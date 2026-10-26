"""
Alignement de peptides : Needleman-Wunsch (global) et Smith-Waterman (local),
avec pénalités de gap affines (Gotoh). Écrit à la main, sans outil externe (T36).

CONVENTION DE GAP (à ne pas changer sans mettre à jour les tests) :
    un gap de longueur L coûte  G_OPEN + (L - 1) * G_EXTEND
    -> le 1er '-' coûte G_OPEN, chaque '-' suivant coûte G_EXTEND.
Valeurs par défaut G_OPEN = 11, G_EXTEND = 1 : valeurs proposées par la
spécification de la semaine 2 (T36).
Attention : deux conventions existent (open + (L-1) x extend, la nôtre, contre
open + L x extend). Toute divergence d'ancre avec un autre outil vient de cette
convention, pas du code.

Score d'un alignement = somme des scores BLOSUM62 des paires alignées
                        - somme des pénalités de gap.

Trois matrices (Gotoh) :
    H[i][j] : meilleur score d'un alignement des i premières lettres de x
              et des j premières de y
    E[i][j] : même chose, mais finissant par un gap dans x (lettre de y seule)
    F[i][j] : même chose, mais finissant par un gap dans y (lettre de x seule)
Local (Smith-Waterman) : même moteur avec max(0, ...) et départ du meilleur H.

Empreinte mémoire : les trois matrices sont allouées en entier, en O(n*m). C'est
négligeable pour la fenêtre de 10-50 résidus (3 x 51 x 51 entiers), mais T35
appellera la fonction sur des milliers de paires : garder un appel par paire
plutôt que de conserver une matrice globale.

Identité par alignement = positions identiques / longueur de l'alignement (gaps inclus).
Elle diffère de l'identité Levenshtein (T23) : celle-ci vaut 1 - distance / max(len),
avec un coût unitaire par édition ; ici les substitutions sont pondérées par BLOSUM62
et les gaps sont affines.
"""
from typing import NamedTuple

from peplib.constants import BLOSUM62, blosum62_score

G_OPEN = 11
G_EXTEND = 1
# Sentinelle "moins l'infini" entière. Le score d'un alignement est borné par
# 50 * 11 = 550 (fenêtre maximale de 50 résidus) et par -(50 * 11 + 50 * 4) = -750 :
# la sentinelle n'est jamais atteinte, elle sert seulement à interdire un état au départ.
_NEG = -10 ** 9


class Alignment(NamedTuple):
    score: int
    aligned_x: str
    aligned_y: str
    x_start: int   # indices dans x : x[x_start:x_end] est la partie alignée
    x_end: int
    y_start: int
    y_end: int


def _check(seq: str, name: str) -> None:
    """Refuse tout résidu hors alphabet canonique, en nommant l'argument fautif.

    Même règle que `blosum62_score` (les clés de BLOSUM62, seule source de
    vérité) : la duplication se limite à l'agrégation des lettres fautives, pour
    signaler en une fois tout ce qui cloche dans une séquence.
    """
    bad = sorted({c for c in seq if c not in BLOSUM62})
    if bad:
        raise ValueError(f"{name} contient des lettres non canoniques : {bad}")


def _align(x: str, y: str, gap_open: int, gap_extend: int, local: bool) -> Alignment:
    _check(x, "x")
    _check(y, "y")
    n, m = len(x), len(y)
    H = [[0] * (m + 1) for _ in range(n + 1)]
    E = [[_NEG] * (m + 1) for _ in range(n + 1)]
    F = [[_NEG] * (m + 1) for _ in range(n + 1)]

    if not local:
        for j in range(1, m + 1):
            E[0][j] = -(gap_open + (j - 1) * gap_extend)
            H[0][j] = E[0][j]
        for i in range(1, n + 1):
            F[i][0] = -(gap_open + (i - 1) * gap_extend)
            H[i][0] = F[i][0]

    best, best_i, best_j = 0, 0, 0
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            E[i][j] = max(H[i][j - 1] - gap_open, E[i][j - 1] - gap_extend)
            F[i][j] = max(H[i - 1][j] - gap_open, F[i - 1][j] - gap_extend)
            h = max(H[i - 1][j - 1] + BLOSUM62[x[i - 1]][y[j - 1]], E[i][j], F[i][j])
            if local:
                h = max(0, h)
            H[i][j] = h
            if local and h > best:
                best, best_i, best_j = h, i, j

    if local:
        i, j, score = best_i, best_j, best
    else:
        i, j, score = n, m, H[n][m]
    end_i, end_j = i, j

    ax, ay, state = [], [], "H"
    while True:
        if local and state == "H" and H[i][j] == 0:
            break
        if not local and i == 0 and j == 0:
            break
        if state == "H":
            if i > 0 and j > 0 and H[i][j] == H[i - 1][j - 1] + BLOSUM62[x[i - 1]][y[j - 1]]:
                ax.append(x[i - 1]); ay.append(y[j - 1]); i -= 1; j -= 1
            elif H[i][j] == E[i][j]:
                state = "E"
            else:
                state = "F"
        elif state == "E":          # gap dans x : on consomme y[j-1]
            ax.append("-"); ay.append(y[j - 1])
            state = "H" if E[i][j] == H[i][j - 1] - gap_open else "E"
            j -= 1
        else:                       # gap dans y : on consomme x[i-1]
            ax.append(x[i - 1]); ay.append("-")
            state = "H" if F[i][j] == H[i - 1][j] - gap_open else "F"
            i -= 1

    return Alignment(score, "".join(reversed(ax)), "".join(reversed(ay)),
                     i, end_i, j, end_j)


def needleman_wunsch(
    x: str, y: str, gap_open: int = G_OPEN, gap_extend: int = G_EXTEND
) -> Alignment:
    """Alignement global de x et y (de bout en bout) : score + alignement."""
    return _align(x, y, gap_open, gap_extend, local=False)


def smith_waterman(
    x: str, y: str, gap_open: int = G_OPEN, gap_extend: int = G_EXTEND
) -> Alignment:
    """Meilleur alignement local de x et y : score >= 0 + sous-chaînes alignées."""
    return _align(x, y, gap_open, gap_extend, local=True)


def score_alignment(
    aligned_x: str, aligned_y: str, gap_open: int = G_OPEN, gap_extend: int = G_EXTEND
) -> int:
    """Recalcule le score d'un alignement déjà construit (contrôle indépendant).

    Contrat : les deux chaînes alignées ont la même longueur et ne contiennent
    que des résidus canoniques ou le gap '-'. Toute autre entrée lève ValueError :
    sans cette garde, zip s'arrête à la chaîne la plus courte et le score porte
    silencieusement sur un préfixe.

    Deux gaps adjacents de sens opposés (un '-' dans x suivi d'un '-' dans y)
    comptent deux ouvertures distinctes : le DP ne produit jamais ce motif, il
    n'apparaît que dans une saisie utilisateur, et cette règle est volontaire.
    """
    if len(aligned_x) != len(aligned_y):
        raise ValueError(
            f"alignements de longueurs différentes : {len(aligned_x)} != {len(aligned_y)}"
        )
    _check(aligned_x.replace("-", ""), "aligned_x")
    _check(aligned_y.replace("-", ""), "aligned_y")
    total, in_gap = 0, None
    for a, b in zip(aligned_x, aligned_y):
        if a == "-" or b == "-":
            kind = "x" if a == "-" else "y"
            total -= gap_extend if in_gap == kind else gap_open
            in_gap = kind
        else:
            total += blosum62_score(a, b)
            in_gap = None
    return total


def alignment_identity(aligned_x: str, aligned_y: str) -> float:
    """Positions identiques / longueur de l'alignement (gaps inclus).

    Contrat : deux chaînes alignées de même longueur ; des longueurs différentes
    lèvent ValueError (sinon zip compare un préfixe et rend une identité fausse).
    Convention : deux chaînes vides valent 0.0, un ratio 0/0 n'ayant pas de sens.
    """
    if len(aligned_x) != len(aligned_y):
        raise ValueError(
            f"alignements de longueurs différentes : {len(aligned_x)} != {len(aligned_y)}"
        )
    if not aligned_x:
        return 0.0
    same = sum(1 for a, b in zip(aligned_x, aligned_y) if a == b and a != "-")
    return same / len(aligned_x)


def identite_couverture(aligned_x: str, aligned_y: str) -> float:
    """Positions identiques / longueur de la séquence la plus courte (sans gaps).

    Convention dite « coverage » : normalisée par min(len_x_sans_gap, len_y_sans_gap).
    C'est le critère retenu pour les clusters v0.2 (seuil 0,90, A5).

    Deux chaînes vides (après retrait des gaps) -> 0.0.
    Longueurs d'alignement différentes -> ValueError.

    Relation avec alignment_identity : identite_couverture >= alignment_identity,
    car le denominateur est au plus egal a la longueur de l'alignement.
    """
    if len(aligned_x) != len(aligned_y):
        raise ValueError(
            f"alignements de longueurs differentes : {len(aligned_x)} != {len(aligned_y)}"
        )
    x_len = len(aligned_x.replace("-", ""))
    y_len = len(aligned_y.replace("-", ""))
    denom = min(x_len, y_len)
    if denom == 0:
        return 0.0
    same = sum(
        1 for a, b in zip(aligned_x, aligned_y)
        if a == b and a != "-"
    )
    return same / denom


def identite_globale(aligned_x: str, aligned_y: str) -> float:
    """Positions identiques / longueur de la séquence la plus longue (sans gaps).

    Convention dite « global » : normalisee par max(len_x_sans_gap, len_y_sans_gap).
    Utilisee comme second critere en A5/T14 (seuil 0,70 identite globale,
    parametrages (o, e) in {(11,1),(10,1),(8,2)}).

    Deux chaînes vides (apres retrait des gaps) -> 0.0.
    Longueurs d'alignement differentes -> ValueError.

    Relation : identite_globale <= identite_couverture <= 1.0.
    """
    if len(aligned_x) != len(aligned_y):
        raise ValueError(
            f"alignements de longueurs differentes : {len(aligned_x)} != {len(aligned_y)}"
        )
    x_len = len(aligned_x.replace("-", ""))
    y_len = len(aligned_y.replace("-", ""))
    denom = max(x_len, y_len)
    if denom == 0:
        return 0.0
    same = sum(
        1 for a, b in zip(aligned_x, aligned_y)
        if a == b and a != "-"
    )
    return same / denom