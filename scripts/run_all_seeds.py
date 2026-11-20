"""
scripts/run_all_seeds.py -- Lance les 6 runs semaine 3 (T53).

Ordre : LSTM 42->43->44, puis Transformer 42->43->44.
Produit results/week3/runs_summary.csv avec moyenne +/- ecart-type par modele.

Usage :
    .venv\\Scripts\\python scripts\\run_all_seeds.py
    .venv\\Scripts\\python scripts\\run_all_seeds.py --only lstm
    .venv\\Scripts\\python scripts\\run_all_seeds.py --only transformer
    .venv\\Scripts\\python scripts\\run_all_seeds.py --resume
"""
from __future__ import annotations

import argparse
import csv
import math
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from peplib.train import train_model

CONFIGS = {
    "lstm":        "configs/lstm_base.yaml",
    "transformer": "configs/transformer_base.yaml",
}
SEEDS = [42, 43, 44]
RESULTS_DIR = Path("results/week3")


def _std(values: list[float]) -> float:
    """Ecart-type non biaise (n-1)."""
    if len(values) < 2:
        return 0.0
    mu = sum(values) / len(values)
    return math.sqrt(sum((x - mu) ** 2 for x in values) / (len(values) - 1))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--only",   choices=["lstm", "transformer"], default=None)
    args = parser.parse_args()

    models_to_run = [args.only] if args.only else ["lstm", "transformer"]
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    all_results: list[dict] = []
    total_t0 = time.perf_counter()

    for model_name in models_to_run:
        for seed in SEEDS:
            print(f"\n{'='*60}")
            print(f"  {model_name.upper()} | seed={seed}")
            print(f"{'='*60}")
            t0 = time.perf_counter()
            try:
                result = train_model(CONFIGS[model_name], seed=seed, resume=args.resume)
                result["status"] = "ok"
            except Exception as exc:
                print(f"[run_all_seeds] ERREUR {model_name} seed={seed}: {exc}")
                result = {
                    "model": model_name, "seed": seed,
                    "best_val_loss": float("nan"),
                    "best_val_pp":   float("nan"),
                    "best_epoch":    0, "status": f"error: {exc}",
                }
            result["total_time_s"] = round(time.perf_counter() - t0, 1)
            all_results.append(result)

    # --- Tableau de resultats ---
    print(f"\n{'='*60}")
    print("  RESULTATS DES RUNS")
    print(f"{'='*60}")
    for r in all_results:
        print(
            f"  {r['model'].upper():12s} seed={r['seed']} | "
            f"val_pp={r.get('best_val_pp', float('nan')):.2f} | "
            f"best_epoch={r.get('best_epoch', '?')} | "
            f"{r.get('total_time_s','?')}s | {r.get('status','?')}"
        )

    # --- Aggregation moyenne +/- ecart-type par modele ---
    print(f"\n{'='*60}")
    print("  AGREGATION (moyenne +/- ecart-type)")
    print(f"{'='*60}")

    summary_rows: list[dict] = []
    for model_name in models_to_run:
        pps = [
            r["best_val_pp"] for r in all_results
            if r["model"] == model_name and math.isfinite(r.get("best_val_pp", float("nan")))
        ]
        if pps:
            mu  = sum(pps) / len(pps)
            sig = _std(pps)
            n   = len(pps)
            print(f"  {model_name.upper():12s} val_pp = {mu:.3f} +/- {sig:.3f} (n={n})")
            summary_rows.append({
                "model":       model_name,
                "n_seeds":     n,
                "val_pp_mean": round(mu, 4),
                "val_pp_std":  round(sig, 4),
                "seeds":       str([r["seed"] for r in all_results if r["model"] == model_name]),
            })

    # CSV detail
    if all_results:
        detail_path = RESULTS_DIR / "runs_detail.csv"
        with open(detail_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(all_results[0].keys()))
            writer.writeheader()
            writer.writerows(all_results)
        print(f"\n[run_all_seeds] Detail : {detail_path}")

    # CSV resume
    if summary_rows:
        summary_path = RESULTS_DIR / "runs_summary.csv"
        with open(summary_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
            writer.writeheader()
            writer.writerows(summary_rows)
        print(f"[run_all_seeds] Resume : {summary_path}")

    total_elapsed = time.perf_counter() - total_t0
    print(f"[run_all_seeds] Duree totale : {total_elapsed/60:.1f} min")


if __name__ == "__main__":
    main()
