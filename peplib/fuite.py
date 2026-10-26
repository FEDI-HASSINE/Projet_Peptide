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


def candidats_val_train(
    seq_val: list,
    seq_train: list,
    tau: float = TAU_CANDIDATS,
    k: int = 3,
    stats: dict | None = None,
) -> set:
    """Candidats val -> train uniquement (borne sonore A6, tau en paramètre).

    Même borne multiset que ``peplib.cluster_v02.candidats_par_couverture``
    mais restreinte aux paires croisées : les groupes de 3-mers sont
    construits séparément pour val et train, et seules les paires (val,
    train) sont accumulées (coût Σ_g |train_g| × |val_g| au lieu de
    Σ_g C(|g|, 2)). Les paires internes val/val et train/train ne sont
    jamais générées. Retour : ensemble de (iv, it), indices respectifs.

    Avec tau = 0,70 : ``kept(0,70) ⊇ kept(0,90)`` (la borne K(m) croît
    avec tau) — un seul passage couvre donc exactement F(0,90) ET F(0,70).
    """
    from peplib.cluster_v02 import borne_intacts

    occurrences_val = []
    groupes_val = {}
    for idx, s in enumerate(seq_val):
        compte = {}
        for pos in range(len(s) - k + 1):
            km = s[pos : pos + k]
            compte[km] = compte.get(km, 0) + 1
            membres = groupes_val.get(km)
            if membres is None:
                groupes_val[km] = [idx]
            elif membres[-1] != idx:
                membres.append(idx)
        occurrences_val.append(compte)

    occurrences_train = []
    groupes_train = {}
    for idx, s in enumerate(seq_train):
        compte = {}
        for pos in range(len(s) - k + 1):
            km = s[pos : pos + k]
            compte[km] = compte.get(km, 0) + 1
            membres = groupes_train.get(km)
            if membres is None:
                groupes_train[km] = [idx]
            elif membres[-1] != idx:
                membres.append(idx)
        occurrences_train.append(compte)

    seuils = {}
    partiel = {}
    retenues = set()
    for km, membres_val in groupes_val.items():
        membres_train = groupes_train.get(km)
        if membres_train is None:
            continue
        for iv in membres_val:
            occ_v = occurrences_val[iv][km]
            long_v = len(seq_val[iv])
            for it in membres_train:
                paire = (iv, it)
                if paire in retenues:
                    continue
                u = partiel.get(paire, 0) + min(
                    occ_v, occurrences_train[it][km]
                )
                m = min(long_v, len(seq_train[it]))
                seuil = seuils.get(m)
                if seuil is None:
                    seuil = borne_intacts(m, tau, k)
                    seuils[m] = seuil
                if u >= seuil:
                    retenues.add(paire)
                    partiel.pop(paire, None)
                else:
                    partiel[paire] = u

    if stats is not None:
        stats["candidats_val_train"] = len(retenues)
    print(f"Candidats val-train (tau = {tau}) : {len(retenues)}", flush=True)
    return retenues


def meilleures_paires(
    seq_val: list, seq_train: list, paires: set
) -> pd.DataFrame:
    """Cœur de mesure : SW sur chaque paire croisée (iv, it) ; max par val.

    Retourne un DataFrame (ordre de val) : best_identite, best_couverture,
    best_min (le max sur les paires de min(identité, couverture) ; 0,0 si
    aucune paire candidate).
    """
    best = {
        i: {"best_identite": 0.0, "best_couverture": 0.0, "best_min": 0.0}
        for i in range(len(seq_val))
    }
    total = len(paires)
    evaluees = 0
    for iv, it in paires:
        identite, couverture = identite_couverture(
            seq_val[iv], seq_train[it]
        )
        minimum = min(identite, couverture)
        b = best[iv]
        if minimum > b["best_min"]:
            b["best_min"] = minimum
            b["best_identite"] = identite
            b["best_couverture"] = couverture
        evaluees += 1
        if evaluees % 20000 == 0:
            print(f"Paires évaluées : {evaluees}/{total}", flush=True)

    retour = pd.DataFrame(
        [best[i] for i in range(len(seq_val))],
        columns=("best_identite", "best_couverture", "best_min"),
    )
    retour.insert(0, "sequence", seq_val)
    return retour.reset_index(drop=True)


def meilleures_couvertures(
    seq_val: list, seq_train: list, tau_candidats: float = TAU_CANDIDATS
) -> pd.DataFrame:
    """Mesure par séquence de val : meilleure paire contre train.

    Candidats calculés au sens GLOBAL (``candidats_par_couverture`` sur
    val + train), puis convertis en paires croisées — les paires internes
    val/val et train/train sont hors courbe. Utilisé par les tests et les
    petits volumes.
    """
    toutes = seq_val + seq_train
    candidats = candidats_par_couverture(toutes, tau=tau_candidats)
    n_val = len(seq_val)
    croisees = set()
    for i, j in candidats:
        if (i < n_val) == (j < n_val):
            continue              # paire interne val ou train : hors courbe
        if i < n_val:
            croisees.add((i, j - n_val))
        else:
            croisees.add((j, i - n_val))
    return meilleures_paires(seq_val, seq_train, croisees)


def fractions_F(stats_par_val: pd.DataFrame, seuils=SEUILS) -> dict:
    """F(θ) : proportion de val dont la meilleure paire atteint (θ ; θ).

    La mesure porte sur LE couple (identité, couverture) de la meilleure
    paire (celle qui maximise min(identité, couverture)) : ``best_min``
    est le min DE CETTE paire, jamais un mélange de valeurs prises sur
    deux paires differentes (invariant verrouillé par les tests).
    """
    lignes = {}
    for seuil in seuils:
        masque = stats_par_val["best_min"] >= seuil
        lignes[f"F_{round(100 * seuil):02d}"] = float(masque.mean())
    return lignes


def _figure(courbes: dict, figure: str) -> None:
    """Histogrammes de best_min par split + repères 0,70 / 0,90."""
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
    plt.xlabel("Meilleure valeur de min(identité, couverture) pour val vs train")
    plt.ylabel("Densité")
    plt.title("Courbe de fuite : val contre train (deux splits)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(figure, dpi=150)
    plt.close()


def _ecritures(
    detail: pd.DataFrame, courbes: dict, sortie_csv: str, figure: str
) -> pd.DataFrame:
    """CSV de détail + tableau des fractions F(θ) + figure."""
    Path(sortie_csv).parent.mkdir(parents=True, exist_ok=True)
    detail.to_csv(sortie_csv, index=False, lineterminator="\n")

    fractions = []
    for nom in courbes:
        f = fractions_F(detail[detail["split"] == nom])
        for nom_fraction, valeur in f.items():
            fractions.append(
                {"split": nom, "fraction": nom_fraction, "valeur": valeur}
            )
    fractions = pd.DataFrame(fractions)
    chemin_fractions = Path(sortie_csv).with_name("fuite_fractions.csv")
    fractions.to_csv(chemin_fractions, index=False, lineterminator="\n")

    _figure(courbes, figure)
    return fractions


def compare_splits(
    val_clusters: list, train_clusters: list,
    val_controle: list, train_controle: list,
    sortie_csv="results/week2/fuite.csv",
    figure="results/week2/courbe_fuite.png",
) -> pd.DataFrame:
    """Fuite mesurée pour les deux splits + CSV + la figure.

    Chaque argument est une liste de séquences (val/train de chaque split).
    Petits volumes : le préfiltre est calculé au sens global ici ; pour
    les artefacts réels, ``run()`` utilise ``candidats_val_train``.
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
    fractions = _ecritures(detail, courbes, sortie_csv, figure)
    return fractions, detail


def run(
    clusters_dir="data/splits",
    controle_dir="data/splits_control",
    sortie_csv="results/week2/fuite.csv",
    figure="results/week2/courbe_fuite.png",
    tau_candidats=TAU_CANDIDATS,
):
    """Courbe de fuite sur les artefacts réels des deux splits.

    Pour chaque split : candidats val->train restreints (coût croisé),
    SW de chaque candidat, meilleur couple par séquence de val, puis
    F(0,70) et F(0,90) exactes (un seul passage, cf. docstring de
    ``candidats_val_train``). Écrit ``fuite.csv`` (détail par séquence),
    ``fuite_fractions.csv`` (F par split) et la figure.
    """
    detail_lignes = []
    courbes = {}
    for nom, repertoire in (
        ("clusters", clusters_dir),
        ("controle", controle_dir),
    ):
        val = pd.read_csv(
            f"{repertoire}/val_v02.csv", dtype=str
        )["sequence"].tolist()
        train = pd.read_csv(
            f"{repertoire}/train_v02.csv", dtype=str
        )["sequence"].tolist()
        st = {}
        candidats = candidats_val_train(
            val, train, tau=tau_candidats, stats=st
        )
        stats = meilleures_paires(val, train, candidats)
        print(f"Split {nom} :", fractions_F(stats), flush=True)
        lignes_split = stats.copy()
        lignes_split.insert(0, "split", nom)
        detail_lignes.append(lignes_split)
        courbes[nom] = stats["best_min"].tolist()

    detail = pd.concat(detail_lignes, ignore_index=True)
    fractions = _ecritures(detail, courbes, sortie_csv, figure)
    return fractions, detail
