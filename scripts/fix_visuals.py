#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Regenere la figure des courbes de perte depuis results/week3/runs.csv.

Aucune valeur n'est codee en dur : la figure ne contient que ce que le CSV
contient. Si le CSV est absent, le script echoue au lieu de produire une
figure plausible.

Usage :
    python scripts/fix_visuals.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent
RUNS_CSV = REPO_ROOT / "results" / "week3" / "runs.csv"
OUT_PNG = REPO_ROOT / "results" / "week3" / "visuals" / "2_courbes_perte.png"

REQUIRED_COLUMNS = ("model", "seed", "epoch", "train_loss", "val_loss", "val_pp")
COLORS = {"lstm": "#4c72b0", "transformer": "#dd8452"}


def load_runs(csv_path: str | Path = RUNS_CSV) -> pd.DataFrame:
    """Charge runs.csv et verifie les colonnes necessaires."""
    csv_path = Path(csv_path)
    if not csv_path.exists():
        raise FileNotFoundError(
            f"{csv_path} introuvable : lancer NB03 (entrainement) avant la figure."
        )
    frame = pd.read_csv(csv_path)
    missing = [c for c in REQUIRED_COLUMNS if c not in frame.columns]
    if missing:
        raise ValueError(f"Colonnes manquantes dans {csv_path} : {missing}")
    if frame.empty:
        raise ValueError(f"{csv_path} est vide")
    return frame


def best_per_run(runs: pd.DataFrame) -> pd.DataFrame:
    """Meilleure epoch (val_loss minimale) de chaque (model, seed)."""
    return runs.loc[runs.groupby(["model", "seed"])["val_loss"].idxmin()]


def build_loss_figure(
    runs: pd.DataFrame,
    out_png: str | Path = OUT_PNG,
) -> pd.DataFrame:
    """Trace train/val par seed pour chaque modele. Retourne les meilleurs runs."""
    best = best_per_run(runs).sort_values(["model", "seed"]).reset_index(drop=True)
    models = sorted(runs["model"].unique())
    out_png = Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)

    figure, axes = plt.subplots(1, len(models), figsize=(6.0 * len(models), 4.4))
    if len(models) == 1:
        axes = [axes]

    for axis, model in zip(axes, models):
        subset = runs[runs["model"] == model]
        color = COLORS.get(model, "#555555")

        for seed in sorted(subset["seed"].unique()):
            per_seed = subset[subset["seed"] == seed].sort_values("epoch")
            axis.plot(
                per_seed["epoch"],
                per_seed["val_loss"],
                color=color,
                alpha=0.45,
                linewidth=1.2,
                label=f"val loss (seed {seed})",
            )

        mean_curve = subset.groupby("epoch")["val_loss"].mean().reset_index()
        axis.plot(
            mean_curve["epoch"],
            mean_curve["val_loss"],
            color=color,
            linewidth=2.4,
            marker="o",
            markersize=3,
            label="val loss (moyenne des seeds)",
        )

        best_model = best[best["model"] == model].iloc[0]
        axis.axvline(
            best_model["epoch"],
            color="black",
            linestyle="--",
            linewidth=1.0,
            label=f"meilleure epoch (seed {int(best_model['seed'])})",
        )
        axis.plot(
            best_model["epoch"],
            best_model["val_loss"],
            marker="*",
            color="black",
            markersize=12,
            zorder=5,
        )
        axis.annotate(
            f"val_loss={best_model['val_loss']:.3f}\nval_pp={best_model['val_pp']:.2f}",
            xy=(best_model["epoch"], best_model["val_loss"]),
            xytext=(8, 10),
            textcoords="offset points",
            fontsize=9,
        )

        axis.set_title(f"{model.upper()} -- perte de validation par seed", fontsize=12, fontweight="bold")
        axis.set_xlabel("Epoch")
        axis.set_ylabel("Cross-entropy (nats/token)")
        axis.grid(True, alpha=0.3)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.legend(fontsize=8)

    figure.suptitle(
        "Courbes de perte reconstruites depuis results/week3/runs.csv",
        fontsize=12,
    )
    figure.tight_layout(rect=(0, 0, 1, 0.96))
    figure.savefig(out_png, dpi=200)
    plt.close(figure)

    return best


def main() -> int:
    runs = load_runs()
    best = build_loss_figure(runs)
    print("[visuals] figure :", OUT_PNG.relative_to(REPO_ROOT).as_posix())
    print("[visuals] meilleure epoch par run :")
    print(best[["model", "seed", "epoch", "val_loss", "val_pp"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
