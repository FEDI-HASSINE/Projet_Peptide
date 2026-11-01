"""
scripts/run_all_seeds.py — Lance les 6 runs de la Semaine 3 (T53).

Ordre : LSTM seed 42 → LSTM 43 → LSTM 44 → Transformer 42 → Transformer 43 → Transformer 44

Reprend depuis best.pt si --resume est passé.
Écrit un résumé dans results/week3/runs_summary.csv à la fin.

Usage :
    .venv\\Scripts\\python scripts\\run_all_seeds.py
    .venv\\Scripts\\python scripts\\run_all_seeds.py --resume
    .venv\\Scripts\\python scripts\\run_all_seeds.py --only lstm
    .venv\\Scripts\\python scripts\\run_all_seeds.py --only transformer
"""
from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path

# Ajouter la racine du projet au PYTHONPATH
sys.path.insert(0, str(Path(__file__).parent.parent))

from peplib.train import train_model

CONFIGS = {
    "lstm":        "configs/lstm_base.yaml",
    "transformer": "configs/transformer_base.yaml",
}
SEEDS = [42, 43, 44]
RESULTS_DIR = Path("results/week3")


def main() -> None:
    parser = argparse.ArgumentParser(description="Lance les 6 runs Semaine 3.")
    parser.add_argument("--resume", action="store_true", help="Reprendre depuis best.pt si disponible")
    parser.add_argument("--only",   choices=["lstm", "transformer"], default=None,
                        help="Entraîner un seul type de modèle")
    args = parser.parse_args()

    models_to_run = [args.only] if args.only else ["lstm", "transformer"]

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    summary_path = RESULTS_DIR / "runs_summary.csv"

    all_results: list[dict] = []
    total_t0 = time.perf_counter()

    for model_name in models_to_run:
        config_path = CONFIGS[model_name]
        for seed in SEEDS:
            print(f"\n{'='*60}")
            print(f"  {model_name.upper()} | seed={seed}")
            print(f"{'='*60}")
            t0 = time.perf_counter()
            try:
                result = train_model(config_path, seed=seed, resume=args.resume)
                result["status"] = "ok"
            except Exception as exc:
                print(f"[run_all_seeds] ERREUR {model_name} seed={seed}: {exc}")
                result = {
                    "model": model_name, "seed": seed,
                    "best_val_loss": float("nan"), "best_val_pp": float("nan"),
                    "best_epoch": 0, "status": f"error: {exc}",
                }
            result["total_time_s"] = round(time.perf_counter() - t0, 1)
            all_results.append(result)

    # --- Résumé ---
    print(f"\n{'='*60}")
    print("  RÉSUMÉ DES 6 RUNS")
    print(f"{'='*60}")
    for r in all_results:
        print(
            f"  {r['model'].upper():12s} seed={r['seed']} | "
            f"val_pp={r.get('best_val_pp', float('nan')):.2f} | "
            f"best_epoch={r.get('best_epoch', '?')} | "
            f"status={r.get('status', '?')} | "
            f"{r.get('total_time_s', '?')}s"
        )

    # Écriture du CSV résumé
    if all_results:
        fieldnames = list(all_results[0].keys())
        with open(summary_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(all_results)
        print(f"\n[run_all_seeds] Résumé écrit dans {summary_path}")

    total_elapsed = time.perf_counter() - total_t0
    print(f"[run_all_seeds] Durée totale : {total_elapsed/60:.1f} min")


if __name__ == "__main__":
    main()
