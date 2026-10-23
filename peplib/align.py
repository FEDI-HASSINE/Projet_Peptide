"""
Alignement de peptides : Needleman-Wunsch (global) et Smith-Waterman (local),
avec penalites de gap affines (Gotoh). Ecrit a la main, sans outil externe (T36).

CONVENTION DE GAP (a ne pas changer sans mettre a jour les tests) :
    un gap de longueur L coute  G_OPEN + (L - 1) * G_EXTEND
    -> le 1er '-' coute G_OPEN, chaque '-' suivant coute G_EXTEND.
Valeurs par defaut G_OPEN = 11, G_EXTEND = 1 : ce sont les valeurs par defaut de
BLASTP pour BLOSUM62. Attention : BLAST compte un gap de longueur L a 11 + L*1,
soit 1 point de plus par gap que notre convention (la notre = BLAST avec open 10).
Toute divergence d'ancre avec un autre outil vient de cette convention, pas du code.

Score d'un alignement = somme des scores BLOSUM62 des paires alignees
                        - somme des penalites de gap.

Trois matrices (Gotoh) :
    H[i][j] : meilleur score d'un alignement des i premieres lettres de x
              et des j premieres de y
    E[i][j] : meme chose, mais finissant par un gap dans x (lettre de y seule)
    F[i][j] : meme chose, mais finissant par un gap dans y (lettre de x seule)
Local (Smith-Waterman) : meme moteur avec max(0, ...) et depart du meilleur H.

Identite par alignement = positions identiques / longueur de l'alignement (gaps inclus).
Elle differe de l'identite Levenshtein (T23) : celle-ci vaut 1 - distance / max(len),
avec un cout unitaire par edition ; ici les substitutions sont ponderees par BLOSUM62
et les gaps sont affines.
"""
from typing import NamedTuple

from peplib.constants import BLOSUM62

G_OPEN = 11
G_EXTEND = 1
_NEG = -10 ** 9  # "moins l'infini" entier


class Alignment(NamedTuple):
    score: int
    aligned_x: str
    aligned_y: str
    x_start: int   # indices dans x : x[x_start:x_end] est la partie alignee
    x_end: int
    y_start: int
    y_end: int


def _check(seq, name):
    bad = sorted({c for c in seq if c not in BLOSUM62})
    if bad:
        raise ValueError(f"{name} contient des lettres non canoniques : {bad}")


def _align(x, y, gap_open, gap_extend, local):
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


def needleman_wunsch(x, y, gap_open=G_OPEN, gap_extend=G_EXTEND):
    """Alignement global de x et y (de bout en bout) : score + alignement."""
    return _align(x, y, gap_open, gap_extend, local=False)


def smith_waterman(x, y, gap_open=G_OPEN, gap_extend=G_EXTEND):
    """Meilleur alignement local de x et y : score >= 0 + sous-chaines alignees."""
    return _align(x, y, gap_open, gap_extend, local=True)


def score_alignment(aligned_x, aligned_y, gap_open=G_OPEN, gap_extend=G_EXTEND):
    """Re-calcule le score d'un alignement deja construit (controle independant)."""
    total, in_gap = 0, None
    for a, b in zip(aligned_x, aligned_y):
        if a == "-" or b == "-":
            kind = "x" if a == "-" else "y"
            total -= gap_extend if in_gap == kind else gap_open
            in_gap = kind
        else:
            total += BLOSUM62[a][b]
            in_gap = None
    return total


def alignment_identity(aligned_x, aligned_y):
    """Positions identiques / longueur de l'alignement (gaps inclus)."""
    if not aligned_x:
        return 0.0
    same = sum(1 for a, b in zip(aligned_x, aligned_y) if a == b and a != "-")
    return same / len(aligned_x)