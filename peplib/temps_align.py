"""Mesure du temps d'alignement et extrapolation « tous-contre-tous » (A9/T9).

Le budget de temps est une exigence de l'ajustement A9 : temps moyen
d'alignement pour des longueurs 10 / 30 / 50, puis extrapolation du coût
d'un alignement tous-contre-tous du corpus v0.2 (42 575 séquences).

Les séquences de mesure sont aléatoires non homologues (worst case
représentatif d'une énumération tous-contre-tous). La graine est dérivée
par la graine 42 et l'ordre d'itération (longueur croissante, NW puis SW)
reste figé pour que la mesure soit reproduisible.
"""

import random
import time
from pathlib import Path

import pandas as pd

from peplib.align import needleman_wunsch, smith_waterman

ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
LONGUEURS = (10, 30, 50)      # fenêtre du corpus : 10 à 50
N_PAIRS_MESURE = 100          # >= 100 paires par ligne (exigé A9)
N_CORPUS_V02 = 42_575         # tailles de l'extrapolation tous-contre-tous
OUT_PATH = "results/week2/temps_align.csv"


def run(
    longueurs=LONGUEURS,
    n_pairs=N_PAIRS_MESURE,
    n_corpus=N_CORPUS_V02,
    out_path=OUT_PATH,
    seed=42,
):
    """Mesure le temps moyen par alignement et l'extrapolation associée.

    Retourne un DataFrame avec une par (longueur, méthode) :
    longueur, methode, ms_par_alignement, paires_tous_contre_tous,
    extrapolation_heures.
    """
    rng = random.Random(seed)     # graine dérivée = 42, ordre d'itération figée
    lignes = []
    paires_tous = n_corpus * (n_corpus - 1) / 2

    for longueur in longueurs:
        seqs = [
            "".join(rng.choices(ALPHABET, k=longueur))
            for _ in range(n_pairs + 1)
        ]
        for methode, fn in (("NW", needleman_wunsch), ("SW", smith_waterman)):
            debut = time.perf_counter()
            for i in range(n_pairs):
                fn(seqs[i], seqs[i + 1])
            total_s = time.perf_counter() - debut
            ms = total_s * 1000 / n_pairs
            heures = ms / 1000 * paires_tous / 3600
            lignes.append(
                {
                    "longueur": longueur,
                    "methode": methode,
                    "ms_par_alignement": round(ms, 4),
                    "paires_tous_contre_tous": int(paires_tous),
                    "extrapolation_heures": round(heures, 2),
                }
            )

    df = pd.DataFrame(lignes)
    sortie = Path(out_path)
    sortie.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(sortie, index=False, lineterminator="\n")
    return df
