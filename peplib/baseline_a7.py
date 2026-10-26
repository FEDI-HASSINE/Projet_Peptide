"""Evaluation finale des baselines n-gram (A7) — stub de livraison.

Ce module sera rempli par l'evaluation sur les VRAIS splits des que la
branche de mission Sofiene (A5/A6 : splits par familles + split aleatoire
de controle) sera fusionnee dans `main` : les CSV de splits ne sont pas
encore presentes sur la base de la PR #20. Conformement au mecanisme de
repli prevu par le plan de l'encadrant, le livrable est DIFFERE avec la
raison ecrite dans le corps de la PR : `run()` leve FileNotFoundError
tant que les splits manquent — AUCUN nombre n'est invente.

Tableau final attendu (Cas 1, apres fusion des splits) :
  - lignes = {ordres 0, 1, 2, 3} x {split familles, split aleatoire} ;
  - colonnes : split, ordre, alpha (grille DEFAULT_ALPHAS = 0.01/0.1/1
    choisie sur val), entropie_nats_token (21 cibles : 20 AA + EOS),
    perplexite ;
  - sorties : `results/week2/baseline_ngram_A7.csv` + figure
    `results/week2/ecart_baselines.png` (entropie versus ordre, 2 courbes),
    avec entree PROVENANCE (commande, date, seed 42, SHA-256).
"""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd

from peplib.ngram import NgramModel, tune_alpha

OUT_PATH = "results/week2/baseline_ngram_A7.csv"
ORDERS = (0, 1, 2, 3)


def _lister_csv(dossier: Path) -> tuple:
    """Retourne (train_csv, val_csv) en refusant tout fichier manquant."""
    train = dossier / "train_v02.csv"
    val = dossier / "val_v02.csv"
    manquants = [str(p) for p in (train, val) if not p.exists()]
    if manquants:
        raise FileNotFoundError(
            "Splits A6 manquants (branche Sofiene a fusionner d'abord) : "
            + ", ".join(manquants)
        )
    return train, val


def _sequences(csv_path: Path) -> list[str]:
    return (
        pd.read_csv(csv_path, dtype=str)["sequence"]
        .dropna()
        .str.strip()
        .astype(str)
        .tolist()
    )


def run(
    splits_dir: str | Path = "data/splits",
    control_dir: str | Path = "data/splits_control",
    out_csv: str | Path = OUT_PATH,
    graines: int = 42,
) -> pd.DataFrame:
    """Evalue les 4 ordres sur les deux splits (familles + aleatoire).

    Leve FileNotFoundError tant que les splits A6 sont absents. Entropie
    croisee en nats/token sur 21 cibles (20 AA + EOS), lissage add-alpha
    avec alpha choisi sur val dans la grille imposee (seed 42).
    """
    decale = Path(splits_dir)
    control = Path(control_dir)
    lignes = []
    for split, dossier in (("familles", decale), ("aleatoire", control)):
        train_csv, val_csv = _lister_csv(dossier)
        train = _sequences(train_csv)
        val = _sequences(val_csv)
        for ordre in ORDERS:
            best_alpha, _ = tune_alpha(
                order=ordre, train_seqs=train, val_seqs=val
            )
            modele = NgramModel(order=ordre, alpha=best_alpha).fit(train)
            ce = modele.cross_entropy(val)
            lignes.append(
                {
                    "split": split,
                    "ordre": ordre,
                    "alpha": best_alpha,
                    "entropie_nats_token": ce,
                    "perplexite": math.exp(ce),
                }
            )

    resultat = pd.DataFrame(lignes)
    out_csv = Path(out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    resultat.to_csv(out_csv, index=False, lineterminator="\n")
    print(f"CSV : {out_csv}", flush=True)
    return resultat
