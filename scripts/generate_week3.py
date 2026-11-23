#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Generation et evaluation sur un checkpoint entraine (Semaine 3, P0-2).

Ce script repond a une question precise : le modele entraine genere-t-il des
sequences nouvelles, ou recopie-t-il le train ?

Il ne filtre pas en silence. Pour chaque checkpoint il enregistre :

- le nombre de tentatives et le taux de validite brut (rejets par cause) ;
- le taux d'unicite et les statistiques de longueur ;
- le nombre de copies exactes du train et de copies en prefixe ;
- la proximite Smith-Waterman au train (0.5 / 0.7 / 0.9) sur un
  sous-echantillon documente, avec intervalle de confiance de Wilson ;
- les proprietes physico-chimiques moyennes (proxys, pas des mesures) ;
- les metadonnees de tracabilite (commit, hash de config, hash du
  checkpoint, device, versions).

Sorties :
    results/week3/generation/generation_<modele>_seed<graine>_<date>.json
    results/week3/generation_summary.csv   (une ligne par run, dedupliquee)

Usage :
    python scripts/generate_week3.py \
        --config configs/lstm_base.yaml \
        --checkpoint checkpoints/lstm/seed_42/best.pt \
        --seed 42 --n 1000
"""
from __future__ import annotations

import argparse
import json
import math
import random
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd
import torch
import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from peplib.descriptors import (  # noqa: E402
    calculate_boman_index,
    calculate_charge,
    calculate_gravy,
    calculate_hydrophobic_moment,
    calculate_pi,
)
from peplib.filters import REJECTION_CAUSES, classify_sequence  # noqa: E402
from peplib.generate import generate_with_stats  # noqa: E402
from peplib.metrics_generation import (  # noqa: E402
    js_divergence,
    sequence_length_stats,
    uniqueness_rate,
)
from peplib.metrics_proximity import proximity_details  # noqa: E402
from peplib.tokenizer import MAX_LEN  # noqa: E402
from peplib.traceability import (  # noqa: E402
    config_sha256,
    environment_info,
    git_commit,
    git_is_dirty,
    sha256_file,
    utc_timestamp,
)
from peplib.train import build_model_from_config  # noqa: E402

SUMMARY_KEY = (
    "model",
    "seed",
    "checkpoint_sha256",
    "config_sha256",
    "n_requested",
    "temperature",
    "top_k",
    "top_p",
)


# ----------------------------------------------------------------------
# Entrees
# ----------------------------------------------------------------------

def load_config(path: str | Path) -> dict:
    with open(path, encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def load_train_sequences(csv_path: str | Path) -> tuple[list[str], int]:
    """Charge les sequences d'entrainement valides (alphabet, 10-50 aa).

    Returns
    -------
    (sequences, n_ecartees) : les sequences retenues et le nombre de lignes
    du CSV qui ne passent pas les filtres (jamaisignorees en silence).
    """
    frame = pd.read_csv(csv_path, dtype=str)
    if "sequence" not in frame.columns:
        raise ValueError(f"Colonne 'sequence' absente de {csv_path}")
    raw = frame["sequence"].dropna().astype(str).str.strip().str.upper().tolist()
    kept = [s for s in raw if classify_sequence(s) is None]
    return kept, len(raw) - len(kept)


def load_checkpoint(model: torch.nn.Module, path: str | Path, device: torch.device) -> dict:
    checkpoint = torch.load(path, map_location=device, weights_only=False)
    if "model_state" not in checkpoint:
        raise ValueError(f"Checkpoint sans 'model_state' : {path}")
    model.load_state_dict(checkpoint["model_state"])
    return checkpoint


# ----------------------------------------------------------------------
# Analyses
# ----------------------------------------------------------------------

def copy_analysis(generated: list[str], train: list[str], prefix_len: int = 10) -> dict:
    """Copies exactes et copies en prefixe d'une sequence du train.

    Une copie en prefixe est une sequence generee G (longueur >= prefix_len)
    telle qu'il existe T dans le train avec ``T.startswith(G)`` et ``T != G``.
    Les deux cas sont des memorisations ; le prefixe est le cas ou le modele
    s'arrete avant la fin d'une sequence apprise.
    """
    if not generated:
        return {
            "n_exact_train": 0,
            "taux_copie_exacte": 0.0,
            "n_prefix_train": 0,
            "taux_copie_prefixe": 0.0,
            "exemples_copie_exacte": [],
            "exemples_copie_prefixe": [],
        }

    train_set = set(train)
    by_prefix: dict[str, list[str]] = defaultdict(list)
    for sequence in train:
        by_prefix[sequence[:prefix_len]].append(sequence)

    exact = [s for s in generated if s in train_set]
    prefixes = []
    for sequence in generated:
        if len(sequence) < prefix_len or sequence in train_set:
            continue
        candidates = by_prefix.get(sequence[:prefix_len])
        if candidates and any(t.startswith(sequence) for t in candidates):
            prefixes.append(sequence)

    return {
        "n_exact_train": len(exact),
        "taux_copie_exacte": len(exact) / len(generated),
        "n_prefix_train": len(prefixes),
        "taux_copie_prefixe": len(prefixes) / len(generated),
        "exemples_copie_exacte": exact[:5],
        "exemples_copie_prefixe": prefixes[:5],
    }


def wilson_interval(successes: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    """Intervalle de confiance 95 % d'une proportion (score de Wilson)."""
    if n <= 0:
        return (0.0, 0.0)
    p_hat = successes / n
    denominator = 1.0 + z * z / n
    centre = (p_hat + z * z / (2 * n)) / denominator
    marge = z * math.sqrt(p_hat * (1 - p_hat) / n + z * z / (4 * n * n)) / denominator
    return (max(0.0, centre - marge), min(1.0, centre + marge))


def proximity_analysis(
    generated: list[str],
    train: list[str],
    seed: int,
    n_prox: int = 200,
    thresholds: tuple[float, ...] = (0.5, 0.7, 0.9),
) -> dict:
    """Proximite au train sur un sous-echantillon deterministe.

    Le prefiltrage par 3-mers est fait au seuil le plus bas (0.5), donc un
    seul passage couvre les trois seuils. Les sequences sans paire candidate
    ont best_min = 0 ; les fractions sont donc des bornes basses
    (cf. docstring de peplib/metrics_proximity.py).
    """
    if not generated or n_prox == 0:
        return {"n_prox": 0, "n_train": len(train), "raison": "aucune sequence"}

    if n_prox < 0 or n_prox >= len(generated):
        sample = list(generated)
    else:
        sample = random.Random(seed).sample(generated, n_prox)

    details = proximity_details(sample, train, threshold=min(thresholds))
    result: dict = {
        "n_prox": len(sample),
        "n_train": len(train),
        "seuil_prefiltre": min(thresholds),
        "prox_moyenne_best_min": float(details["best_min"].mean()),
        "prox_max_best_min": float(details["best_min"].max()),
    }

    for threshold in thresholds:
        count = int((details["best_min"] >= threshold).sum())
        low, high = wilson_interval(count, len(sample))
        tag = f"{threshold:.1f}"
        result[f"n_proche_{tag}"] = count
        result[f"fraction_proche_{tag}"] = count / len(sample)
        result[f"fraction_proche_{tag}_ic95_bas"] = low
        result[f"fraction_proche_{tag}_ic95_haut"] = high

    return result


def property_summary(sequences: list[str]) -> dict:
    """Moyennes des descripteurs (proxys physico-chimiques)."""
    if not sequences:
        return {}
    values = {
        "charge_moyenne": [],
        "gravy_moyen": [],
        "moment_hydrophobe_moyen": [],
        "pi_moyen": [],
        "boman_moyen": [],
    }
    for sequence in sequences:
        values["charge_moyenne"].append(calculate_charge(sequence))
        values["gravy_moyen"].append(calculate_gravy(sequence))
        values["moment_hydrophobe_moyen"].append(calculate_hydrophobic_moment(sequence))
        values["pi_moyen"].append(calculate_pi(sequence))
        values["boman_moyen"].append(calculate_boman_index(sequence))

    return {
        name: round(sum(items) / len(items), 4)
        for name, items in values.items()
    }


def evaluate_checkpoint(
    config_path: str | Path,
    checkpoint_path: str | Path,
    seed: int = 42,
    n: int = 1000,
    temperature: float = 1.0,
    top_k: int = 0,
    top_p: float = 1.0,
    max_len: int = MAX_LEN,
    n_prox: int = 200,
    max_attempts: int | None = None,
    train_csv: str | Path = "data/splits/train_v02.csv",
    out_dir: str | Path = "results/week3",
    device: str | torch.device | None = None,
) -> dict:
    """Genere, mesure et ecrit le rapport d'un checkpoint. Retourne le rapport."""
    config_path = Path(config_path)
    checkpoint_path = Path(checkpoint_path)
    out_dir = Path(out_dir)

    if not checkpoint_path.exists():
        raise FileNotFoundError(f"Checkpoint introuvable : {checkpoint_path}")

    config = load_config(config_path)
    model_name = config.get("model", {}).get("type", "unknown")
    device = torch.device(device) if device is not None else torch.device(
        "cuda" if torch.cuda.is_available() else "cpu"
    )

    from peplib.seed import set_seed

    set_seed(seed)

    model = build_model_from_config(config).to(device)
    checkpoint = load_checkpoint(model, checkpoint_path, device)

    train_sequences, n_train_ecartees = load_train_sequences(train_csv)

    stats = generate_with_stats(
        model=model,
        n=n,
        temperature=temperature,
        top_k=top_k if top_k > 0 else None,
        top_p=top_p,
        max_len=max_len,
        seed=seed,
        device=device,
        max_attempts=max_attempts,
    )
    sequences = stats["sequences"]

    diversity = {
        "taux_unicite": uniqueness_rate(sequences) if sequences else 0.0,
        "n_uniques": len(set(sequences)),
        "longueurs": sequence_length_stats(sequences) if sequences else {},
        "js_vs_train": js_divergence(sequences, train_sequences) if sequences else None,
    }

    copies = copy_analysis(sequences, train_sequences)
    proximity = proximity_analysis(sequences, train_sequences, seed=seed, n_prox=n_prox)
    properties = property_summary(sequences)
    env = environment_info()

    report = {
        "metadata": {
            "timestamp_utc": utc_timestamp(),
            "git_sha": git_commit(REPO_ROOT),
            "git_dirty": git_is_dirty(REPO_ROOT),
            "config_path": config_path.as_posix(),
            "config_sha256": config_sha256(config),
            "checkpoint": checkpoint_path.as_posix(),
            "checkpoint_sha256": sha256_file(checkpoint_path),
            "epoch": int(checkpoint.get("epoch", -1)) + 1,
            "val_loss": float(checkpoint.get("val_loss", float("nan"))),
            "val_pp": float(checkpoint.get("val_pp", float("nan"))),
            "checkpoint_seed": checkpoint.get("seed"),
            "device": device.type,
            **env,
        },
        "sampling": {
            "temperature": temperature,
            "top_k": top_k,
            "top_p": top_p,
            "max_len": max_len,
            "n_requested": n,
        },
        "raw": {
            "n_attempts": stats["n_attempts"],
            "n_raw": stats["n_raw"],
            "n_valid": stats["n_valid"],
            "n_invalid": stats["n_invalid"],
            "taux_validite_brut": stats["n_valid"] / stats["n_attempts"]
            if stats["n_attempts"]
            else 0.0,
            "rejections": stats["rejections"],
            "budget_epuise": stats["exhausted"],
            "max_attempts": stats["max_attempts"],
            "n_train_sequences": len(train_sequences),
            "n_train_ecartees": n_train_ecartees,
        },
        "diversite": diversity,
        "copies": copies,
        "proximite": proximity,
        "proprietes": properties,
        "sequences": sequences,
    }

    # --- Ecriture JSON ---
    stamp = report["metadata"]["timestamp_utc"].replace(":", "").replace("-", "")
    report_dir = out_dir / "generation"
    report_dir.mkdir(parents=True, exist_ok=True)
    json_path = report_dir / f"generation_{model_name}_seed{seed}_{stamp}.json"
    with open(json_path, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, ensure_ascii=False)

    # --- Ligne de synthese ---
    row = {
        "timestamp_utc": report["metadata"]["timestamp_utc"],
        "git_sha": report["metadata"]["git_sha"],
        "git_dirty": report["metadata"]["git_dirty"],
        "model": model_name,
        "seed": seed,
        "device": report["metadata"]["device"],
        "gpu_name": report["metadata"]["gpu_name"] or "",
        "torch_version": report["metadata"]["torch_version"],
        "config_sha256": report["metadata"]["config_sha256"],
        "checkpoint": report["metadata"]["checkpoint"],
        "checkpoint_sha256": report["metadata"]["checkpoint_sha256"],
        "epoch": report["metadata"]["epoch"],
        "val_loss": report["metadata"]["val_loss"],
        "val_pp": report["metadata"]["val_pp"],
        "temperature": temperature,
        "top_k": top_k,
        "top_p": top_p,
        "n_requested": n,
        "n_attempts": stats["n_attempts"],
        "n_valid": stats["n_valid"],
        "n_invalid": stats["n_invalid"],
        "taux_validite_brut": report["raw"]["taux_validite_brut"],
        "budget_epuise": stats["exhausted"],
        **{f"rejet_{cause}": stats["rejections"].get(cause, 0) for cause in REJECTION_CAUSES},
        "n_exact_train": copies["n_exact_train"],
        "taux_copie_exacte": copies["taux_copie_exacte"],
        "n_prefix_train": copies["n_prefix_train"],
        "taux_copie_prefixe": copies["taux_copie_prefixe"],
        "n_prox": proximity.get("n_prox", 0),
        "prox_moyenne_best_min": proximity.get("prox_moyenne_best_min"),
        "prox_max_best_min": proximity.get("prox_max_best_min"),
    }
    for threshold in ("0.5", "0.7", "0.9"):
        for suffix in ("", "_ic95_bas", "_ic95_haut"):
            row[f"fraction_proche_{threshold}{suffix}"] = proximity.get(
                f"fraction_proche_{threshold}{suffix}"
            )
    row.update(
        {
            "taux_unicite": diversity["taux_unicite"],
            "longueur_moyenne": diversity["longueurs"].get("mean"),
            "longueur_min": diversity["longueurs"].get("min"),
            "longueur_max": diversity["longueurs"].get("max"),
            "js_vs_train": diversity["js_vs_train"],
            **properties,
            "json_path": json_path.relative_to(REPO_ROOT).as_posix()
            if json_path.is_relative_to(REPO_ROOT)
            else json_path.as_posix(),
        }
    )

    summary_csv = out_dir / "generation_summary.csv"
    write_summary_row(summary_csv, row)

    report["json_path"] = row["json_path"]
    return report


def write_summary_row(csv_path: str | Path, row: dict) -> None:
    """Ajoute une ligne en remplacant une eventuelle ligne de meme cle.

    La cle (modele, graine, checkpoint, config, n, temperature, top_k, top_p)
    rend le fichier idempotent : relancer la meme evaluation met a jour la
    ligne au lieu de la dupliquer.
    """
    csv_path = Path(csv_path)
    csv_path.parent.mkdir(parents=True, exist_ok=True)

    if csv_path.exists():
        existing = pd.read_csv(csv_path)
        mask = pd.Series(True, index=existing.index)
        for field in SUMMARY_KEY:
            if field not in existing.columns:
                mask &= False
            elif pd.isna(row.get(field)):
                mask &= existing[field].isna()
            else:
                mask &= existing[field] == row[field]
        existing = existing.loc[~mask]
        combined = pd.concat([existing, pd.DataFrame([row])], ignore_index=True)
    else:
        combined = pd.DataFrame([row])

    combined.to_csv(csv_path, index=False)


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generation + evaluation sur checkpoint (Semaine 3)."
    )
    parser.add_argument("--config", required=True, help="YAML de configuration")
    parser.add_argument("--checkpoint", required=True, help="best.pt du run")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--n", type=int, default=1000, help="sequences valides visees")
    parser.add_argument("--temperature", type=float, default=1.0)
    parser.add_argument("--top-k", type=int, default=0, help="0 = desactive")
    parser.add_argument("--top-p", type=float, default=1.0)
    parser.add_argument("--max-len", type=int, default=MAX_LEN)
    parser.add_argument(
        "--n-prox",
        type=int,
        default=200,
        help="sous-echantillon pour Smith-Waterman (0 = sauter, <0 = tout)",
    )
    parser.add_argument(
        "--max-attempts",
        type=int,
        default=None,
        help="budget de tentatives (defaut : max(10*n, 5000))",
    )
    parser.add_argument("--train-csv", default="data/splits/train_v02.csv")
    parser.add_argument("--out-dir", default="results/week3")
    parser.add_argument("--device", default=None, help="cpu / cuda (defaut : auto)")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    print("[generation] checkpoint :", args.checkpoint)
    report = evaluate_checkpoint(
        config_path=args.config,
        checkpoint_path=args.checkpoint,
        seed=args.seed,
        n=args.n,
        temperature=args.temperature,
        top_k=args.top_k,
        top_p=args.top_p,
        max_len=args.max_len,
        n_prox=args.n_prox,
        max_attempts=args.max_attempts,
        train_csv=args.train_csv,
        out_dir=args.out_dir,
        device=args.device,
    )

    raw = report["raw"]
    copies = report["copies"]
    proximity = report["proximite"]

    print(f"[generation] tentatives={raw['n_attempts']} valides={raw['n_valid']} "
          f"taux_validite_brut={raw['taux_validite_brut']:.3f}")
    print(f"[generation] rejets={raw['rejections']}")
    print(f"[generation] copies exactes={copies['n_exact_train']} "
          f"prefixes={copies['n_prefix_train']}")
    if proximity.get("n_prox"):
        for tag in ("0.5", "0.7", "0.9"):
            print(f"[generation] proche>={tag} : "
                  f"{proximity[f'fraction_proche_{tag}']:.3f} "
                  f"(n={proximity['n_prox']})")
    print("[generation] rapport :", report["json_path"])
    if raw["budget_epuise"]:
        print("[generation] ATTENTION : budget de tentatives epuise avant n valides.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
