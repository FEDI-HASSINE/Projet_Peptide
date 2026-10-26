"""Vérification exhaustive des alignements NW et SW (A9 / T8, AJUSTEMENTS §3).

Principe : pour des séquences courtes (longueur <= 4, alphabet réduit {A, C, G}),
l'alignement optimal peut être obtenu en énumérant TOUS les chemins d'édition
possibles, sans programmation dynamique. Le score le plus haut parmi les chemins
est l'optimum exact ; toute divergence avec ``needleman_wunsch`` ou
``smith_waterman`` est un contre-exemple.

Rôle de l'agent : vérificateur. Le module ``peplib/align.py`` (Youssef, PR #13)
n'est pas modifié ; une divergence échoue ici et fournit les paires
contre-exemples pour la PR.

Convention de gap cohérente avec celle du module : un gap de longueur L coûte
G_OPEN + (L - 1) * G_EXTEND ; cette convention est appliquée via le score
indépendant ``score_alignment`` (re-nomme ici ``score``).
"""

import random

from peplib.align import needleman_wunsch, smith_waterman
from peplib.align import score_alignment
from peplib.seed import set_seed

ALPHABET_REDUIT = "ACG"
N_PAIRS = 500


def _chemins_edition(x, y):
    """Énumère les (aligné_x, aligné_y) de tous les chemins d'édition globaux.

    Les gaps de sens opposés adjacents ('-' sur une ligne puis l'autre) ne sont
    jamais générés : le DP n'en produit pas et ``score_alignment`` les coterait
    comme deux ouvertures distinctes, ce qui fausserait la comparaison.
    """
    chemins = []

    def recursion(i, j, ax, ay, dernier):
        if i == len(x) and j == len(y):
            chemins.append((ax, ay))
            return
        if i < len(x) and j < len(y):
            recursion(i + 1, j + 1, ax + x[i], ay + y[j], None)
        if j < len(y) and dernier != "y":       # gap dans x : consomme y
            recursion(i, j + 1, ax + "-", ay + y[j], "x")
        if i < len(x) and dernier != "x":       # gap dans y : consomme x
            recursion(i + 1, j, ax + x[i], ay + "-", "y")

    recursion(0, 0, "", "", None)
    return chemins


def _score_max_global(x, y):
    """Optimum exact au sens global : meilleur chemin mis au score indépendant."""
    return max(score_alignment(ax, ay) for ax, ay in _chemins_edition(x, y))


def _paires_aleatoires(rng, n):
    """n paires de séquences de longueur 1..4 sur l'alphabet réduit."""
    paires = []
    for _ in range(n):
        a = "".join(rng.choices(ALPHABET_REDUIT, k=rng.randint(1, 4)))
        b = "".join(rng.choices(ALPHABET_REDUIT, k=rng.randint(1, 4)))
        paires.append((a, b))
    return paires


def _sous_chaines(s):
    return [s[i:j] for i in range(len(s)) for j in range(i + 1, len(s) + 1)]


def test_verif_exhaustive_scores_globaux():
    """NW au sens global == énumération exhaustive des chemins (>= paires)."""
    set_seed(42)
    rng = random
    for a, b in _paires_aleatoires(rng, N_PAIRS):
        optimum = _score_max_global(a, b)
        al = needleman_wunsch(a, b)
        assert al.score == optimum, (
            f"contre-exemple global : ({a!r}, {b!r}) "
            f"NW={al.score} vs énumération={optimum}"
        )
        # Cohérence interne : l'alignement rendu reproduit bien son score.
        assert score_alignment(al.aligned_x, al.aligned_y) == al.score, (
            f"score incohérent avec la construction : ({a!r}, {b!r})"
        )


def test_verif_exhaustive_scores_locaux():
    """SW == max exact sur toutes les paires de sous-chaînes de a et de b.

    Un alignement local optimal est un alignement global de deux sous-chaînes ;
    réciproquement tout alignement global de sous-chaînes est un local, donc le
    maximum énuméré ici est l'optimum SW exact. Le plancher SW = 0 (alignement
    local vide autorisé, docstring du module) est inclu dans l'énumération :
    il survient quand toute colonne possible est négative (ex. ('G','C') -> -3).
    """
    set_seed(42)
    rng = random
    for a, b in _paires_aleatoires(rng, N_PAIRS):
        optimum_local = max(
            0,
            max(
                _score_max_global(sa, sb)
                for sa in _sous_chaines(a)
                for sb in _sous_chaines(b)
            ),
        )
        al = smith_waterman(a, b)
        assert al.score == optimum_local, (
            f"contre-exemple local : ({a!r}, {b!r}) "
            f"SW={al.score} vs énumération={optimum_local}"
        )
