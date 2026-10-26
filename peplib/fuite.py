"""Courbe de fuite (ajustement A6, livrable D3) : val contre train.

Pour chaque séquence de val, la MEILLEURE identité de couverture atteignable
contre train est mesurée parmi les candidats du préfiltre (borne sonore A6,
tau = 0,70 pour couvrir à la fois F(0,90) et F(0,70)). Pour chaque paire
candidate, Smith-Waterman (o = 11, e = 1, BLOSUM62) fournit le couple
(identité, couverture) et l'on retient par séquence val le couple maximisant
la mesure ``min(identité, couverture)`` (le critère de fuite exige que DEUX
valeurs soient hautes ; min est l'agrégation monotone la plus simple —
choix documenté).

Fractions rapportées :
    F(θ) = proportion des séquences de val possédant au moins une paire
    candidate avec (identité ≥ θ ET couverture ≥ θ).

Le split par familles doit donner F(0,90) = 0 (aucune famille à cheval) ;
le split aléatoire de contrôle, attendu fuyard, donne la valeur de référence.

Sorties : ``results/week2/fuite.csv`` (par séquence val : split, meilleure
identity, meilleure couverture, meilleur min) et
``results/week2/courbe_fuite.png`` (deux histogrammes + repères 0,70/0,90).
"""

from pathlib import Path

import pandas as pd

from peplib.cluster_v02 import candidats_par_couverture
from peplib.identite import identite_couverture

TAU_CANDIDATS = 0.70      # borne de préfiltre la plus basse mesure demandée
SEUILS = (0.70, 0.90)


def meilleures_couvertures(
    seq_val: list, seq_train: list, tau_candidats: float = TAU_CANDIDATS
) -> pd.DataFrame:
    """Mesure par sequence de val : meilleure paire contre train.

    Retourne un DataFrame indexé par séquence de val (ordre d'entrée) avec
    colonnes : best_identite, best_couverture, best_min (le max sur les
    paires candidates de min(identité, couverture) ; 0,0 sans candidat).
    """
    toutes = seq_val + seq_train
    candidats = candidats_par_couverture(toutes, tau=tau_candidats)

    best = {
        i: {"best_identite": 0.0, "best_couverture": 0.0, "best_min": 0.0}
        for i in range(len(seq_val))
    }
    n_val = len(seq_val)
    evaluees = 0
    total = len(candidats)
    for i, j in candidats:
        i_val = i < n_val
        j_val = j < n_val
        if i_val and j_val:
            continue          # paire interne à val : hors courbe
        if not i_val and not j_val:
            continue          # paire interne à train : hors courbe
        if i_val:
            iv, it = i, j - n_val
        else:
            iv, it = j, i - n_val   # les deux offsets tombent sur train
        identite, couverture = identite_couverture(
            seq_val[iv], seq_train[it]
        )
        minimum = min(identite, couverture)
        b = best[iv]
        if minimum > b["best_min"]:
            b["best_min"] = minimum
            b["best_identite"] = identite
            b["best_couverture"] = couverture
        elif minimum == b["best_min"] and minimum > 0:
            # franche égalité de min : on garde les max purs (identité/couv)
            if identite > b["best_identite"]:
                b["best_identite"] = identite
            if couverture > b["best_couverture"]:
                b["best_couverture"] = couverture
        evaluees += 1
        if evaluees % 20000 == 0:
            print(f"Paires évaluées : {evaluees}/{total}", flush=True)

    retour = pd.DataFrame(
        [best[i] for i in range(n_val)],
        columns=("best_identite", "best_couverture", "best_min"),
    )
    retour.insert(0, "sequence", seq_val)
    return retour.reset_index(drop=True)


def fractions_F(stats_par_val: pd.DataFrame, seuils=SEUILS) -> dict:
    """F(θ) : proportion de val dont la meilleure paire atteint (θ ; θ)."""
    lignes = {}
    for seuil in seuils:
        masque = (stats_par_val["best_identite"] >= seuil) & (
            stats_par_val["best_couverture"] >= seuil
        )
        lignes[f"F_{round(100 * seuil):02d}"] = float(masque.mean())
    return lignes


def compare_splits(
    val_clusters: list, train_clusters: list,
    val_controle: list, train_controle: list,
    sortie_csv="results/week2/fuite.csv",
    figure="results/week2/courbe_fuite.png",
) -> pd.DataFrame:
    """Fuite mesurée pour les deux splits + CSV + la figure.

    Chaque argument est une liste de séquences (val/train de chaque split).
    """
    lignes = []
    courbes = {}
    for nom, val, train in (
        ("clusters", val_clusters, train_clusters),
        ("controle", val_controle, train_controle),
    ):
        stats = meilleures_couvertures(val, train)
        lignes_split = stats.copy()
        lignes_split.insert(0, "split", nom)
        lignes.append(lignes_split)
        courbes[nom] = stats["best_min"].tolist()

    detail = pd.concat(lignes, ignore_index=True)
    Path(sortie_csv).parent.mkdir(parents=True, exist_ok=True)
    detail.to_csv(sortie_csv, index=False, lineterminator="\n")

    # F(θ) par split : les fractions exploitable en tableau.
    fractions = []
    for nom, stats in (
        ("clusters", detail[detail["split"] == "clusters"]),
        ("controle", detail[detail["split"] == "controle"]),
    ):
        f = fractions_F(stats)
        for nom_fraction, valeur in f.items():
            fractions.append(
                {"split": nom, "fraction": nom_fraction, "valeur": valeur}
            )
    fractions = pd.DataFrame(fractions)

    import matplotlib.pyplot as plt

    Path(figure).parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(9, 5))
    for nom, valeurs in courbes.items():
        plt.hist(
            valeurs, bins=40, alpha=0.6, edgecolor="black",
            label=f"split {nom}", density=True,
        )
    for seuil in (0.70, 0.90):
        plt.axvline(seuil, linestyle="--", linewidth=1,
                    label=f"seuil {seuil:.2f}")
    plt.xlabel("Meilleure couverture identique min(identité, couverture) pour val vs train")
    plt.ylabel("Densité")
    plt.title("Courbe de fuite : val contre train (deux splits)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(figure, dpi=150)
    plt.close()

    return fractions, detail
