"""Variantes du clustering v0.2 à rapporter (pas à verrouiller) — ajustement A5.

Tableau demandé en revue :
  - critère couverture 0,90 aux conventions de gap (11,1), (10,1) et (8,2) ;
  - critère identité globale NW (identique ≥ 0,90, même três conventions) ;
  - critère couverture 0,70 (o = 11, e = 1).
Le préfiltre multiset est calculé UNE SEULE FOIS (la borne A6 dépend de tau
seulement) et réutilisé pour chaque configuration ; le coût de chaque config
se limite aux alignements SW/NW sur les candidats (1 061 432 paires du run
principal).
"""

import pandas as pd

from peplib.align import needleman_wunsch
from peplib.cluster import UnionFind
from peplib.cluster_v02 import (
    TAU,
    cluster_v02,
    candidats_par_couverture,
    tableau_t13,
)
from peplib.identite import identite_globale

OUT_PATH = "results/week2/variantes_clusters_v02_t13.csv"
CONVENTIONS_GAP = ((11, 1), (10, 1), (8, 2))


def etiquettes_depuis_paires(sequences, paires) -> list:
    """Composantes connexes sur un ensemble de paires (i, j) déjà filtées."""
    unions = UnionFind(len(sequences))
    for i, j in sorted(paires):
        unions.union(i, j)
    racines = [unions.find(i) for i in range(len(sequences))]
    numerotation = {r: n for n, r in enumerate(sorted(set(racines)))}
    return [numerotation[r] for r in racines]


def run(input_path="data/processed/corpus_v0.2.csv", out_path=OUT_PATH):
    """Calcule et écrit le tableau des variantes T13 (console + CSV).

    Configurations exécutées : couverture 0,90 pour (o,e) ∈ {(11,1),(10,1),
    (8,2)} + identité globale NW 0,90 (o=11, e=1). La variante couverture
    0,70 SUR LE CORPUS ENTIER est différée : à tau = 0,70 la borne K(m) =
    ⌊0,49·m − 3⌋ tombe à 1 pour m = 10, ce qui garde en candidat toute
    paire partageant une seule occurrence de 3-mer (dizaines de millions
    de SW) — en revanche le seuil 0,70 EST mesuré exactement par la courbe
    de fuite (``peplib/fuite.py``), sur val->train.
    """
    df = pd.read_csv(input_path, dtype=str)
    sequences = sorted(df["sequence"].unique().tolist())
    print(f"Séquences uniques : {len(sequences)}", flush=True)

    stats = {}
    candidats = candidats_par_couverture(sequences, tau=TAU, stats=stats)
    print(f"Candidats réutilisés : {len(candidats)}", flush=True)

    lignes = []

    for o, e in CONVENTIONS_GAP:
        print(f"Variante couverture 0,90 gaps ({o},{e})...", flush=True)
        etiquettes = cluster_v02(
            sequences, tau=TAU, gap_open=o, gap_extend=e, candidats=candidats
        )
        t = tableau_t13(etiquettes, sortie=None)
        lignes.append(
            {
                "critere": "couverture",
                "tau": 0.90,
                "gap": f"{o},{e}",
                "nb_familles": int(t["nb_familles"][0]),
                "taille_max": int(t["taille_max"][0]),
                "singletons": int(t["singletons"][0]),
                "fraction_plus_grande": float(t["fraction_plus_grande"][0]),
            }
        )
        print(t.to_string(index=False), flush=True)

    # Identité globale NW ≥ 0,90 (mesure complémentaire §5) — sans couverture.
    print("Variante identité globale 0,90...", flush=True)
    retenues = set()
    for i, j in sorted(candidats):
        if identite_globale(sequences[i], sequences[j]) >= 0.90:
            retenues.add((i, j))
    etiquettes = etiquettes_depuis_paires(sequences, retenues)
    t = tableau_t13(etiquettes, sortie=None)
    lignes.append(
        {
            "critere": "identite_globale",
            "tau": 0.90,
            "gap": f"{CONVENTIONS_GAP[0][0]},{CONVENTIONS_GAP[0][1]}",
            "nb_familles": int(t["nb_familles"][0]),
            "taille_max": int(t["taille_max"][0]),
            "singletons": int(t["singletons"][0]),
            "fraction_plus_grande": float(t["fraction_plus_grande"][0]),
        }
    )
    print(t.to_string(index=False), flush=True)

    resultat = pd.DataFrame(lignes)
    resultat.to_csv(out_path, index=False, lineterminator="\n")
    print(f"CSV : {out_path}", flush=True)
    return resultat
