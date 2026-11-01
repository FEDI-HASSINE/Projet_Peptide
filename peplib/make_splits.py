"""
peplib/make_splits.py — Découpe corpus_v0.2.csv en train/val/test (T52 pré-requis).

Ratios : 80 / 10 / 10, seed=42, shuffle reproductible.
Idempotent : ne recrée rien si les trois fichiers existent déjà.

Usage :
    python -m peplib.make_splits                          # par défaut
    python -m peplib.make_splits --corpus data/processed/corpus_v0.2.csv
"""
from __future__ import annotations

import argparse
import random
from pathlib import Path

import pandas as pd


TRAIN_FRAC = 0.80
VAL_FRAC   = 0.10
TEST_FRAC  = 0.10
SEED       = 42


def make_splits(
    corpus_path: str | Path = "data/processed/corpus_v0.2.csv",
    splits_dir:  str | Path = "data/processed/splits",
    seed: int = SEED,
    force: bool = False,
) -> dict[str, Path]:
    """Découpe le corpus en train/val/test et écrit les CSV.

    Parameters
    ----------
    corpus_path : chemin du CSV source (doit contenir une colonne 'sequence').
    splits_dir  : répertoire de destination.
    seed        : graine aléatoire pour la reproductibilité.
    force       : si True, recrée les splits même s'ils existent déjà.

    Returns
    -------
    dict avec les clés 'train', 'val', 'test' et les Path correspondants.
    """
    corpus_path = Path(corpus_path)
    splits_dir  = Path(splits_dir)

    train_path = splits_dir / "train.csv"
    val_path   = splits_dir / "val.csv"
    test_path  = splits_dir / "test.csv"

    if not force and train_path.exists() and val_path.exists() and test_path.exists():
        print(f"[make_splits] Splits déjà présents dans {splits_dir} — rien à faire.")
        return {"train": train_path, "val": val_path, "test": test_path}

    if not corpus_path.exists():
        raise FileNotFoundError(f"Corpus introuvable : {corpus_path}")

    df = pd.read_csv(corpus_path, dtype=str)
    if "sequence" not in df.columns:
        raise ValueError(f"Colonne 'sequence' absente de {corpus_path}")

    # Garder uniquement les lignes valides (non nulles, longueur 10-50)
    df = df.dropna(subset=["sequence"]).copy()
    df["sequence"] = df["sequence"].str.strip().str.upper()
    mask = df["sequence"].str.match(r'^[ACDEFGHIKLMNPQRSTVWY]{10,50}$')
    df = df[mask].reset_index(drop=True)

    indices = list(range(len(df)))
    rng = random.Random(seed)
    rng.shuffle(indices)

    n = len(indices)
    n_train = int(n * TRAIN_FRAC)
    n_val   = int(n * VAL_FRAC)

    train_idx = indices[:n_train]
    val_idx   = indices[n_train : n_train + n_val]
    test_idx  = indices[n_train + n_val :]

    splits_dir.mkdir(parents=True, exist_ok=True)

    df_train = df.iloc[train_idx][["sequence"]].reset_index(drop=True)
    df_val   = df.iloc[val_idx][["sequence"]].reset_index(drop=True)
    df_test  = df.iloc[test_idx][["sequence"]].reset_index(drop=True)

    df_train.to_csv(train_path, index=False)
    df_val.to_csv(val_path,   index=False)
    df_test.to_csv(test_path, index=False)

    print(
        f"[make_splits] Splits créés dans {splits_dir}\n"
        f"  train : {len(df_train):>6} séquences\n"
        f"  val   : {len(df_val):>6} séquences\n"
        f"  test  : {len(df_test):>6} séquences\n"
        f"  total : {len(df):>6} séquences (seed={seed})"
    )
    return {"train": train_path, "val": val_path, "test": test_path}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Crée les splits train/val/test.")
    parser.add_argument("--corpus",     default="data/processed/corpus_v0.2.csv")
    parser.add_argument("--splits-dir", default="data/processed/splits")
    parser.add_argument("--seed",       type=int, default=SEED)
    parser.add_argument("--force",      action="store_true")
    args = parser.parse_args()
    make_splits(args.corpus, args.splits_dir, args.seed, args.force)
