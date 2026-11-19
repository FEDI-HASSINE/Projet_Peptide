
"""T61 — Comparaison des baselines sur validation v0.2.

Execution :
    python -m peplib.evaluate_baselines_s3

Utilise uniquement train_v02.csv et val_v02.csv.
Ne lit jamais test_v02.csv.
"""

from pathlib import Path

import pandas as pd

from peplib.baselines_generation import (
    load_sequences,
    build_reference_baselines,
)
from peplib.metrics_generation import generation_metrics


TRAIN_PATH = Path("data/splits/train_v02.csv")
VAL_PATH = Path("data/splits/val_v02.csv")
OUTPUT_DIR = Path("results/week3")


def main():
    print("=== T61 : Evaluation des baselines ===")

    # 1. Charger les sequences
    train = load_sequences(TRAIN_PATH)
    validation = load_sequences(VAL_PATH)

    print(f"Train : {len(train)} peptides")
    print(f"Validation : {len(validation)} peptides")

    # 2. Construire les cinq baselines
    baselines = build_reference_baselines(
        train=train,
        validation=validation,
        seed=42,
        alpha=1.0,
    )

    # 3. Creer le dossier de resultats
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # 4. Calculer les metriques
    results_rows = []

    for name, sequences in baselines.items():
        print(f"Evaluation : {name}")

        metrics = generation_metrics(
            generated=sequences,
            reference=validation,
        )

        results_rows.append({
            "baseline": name,
            **metrics,
        })

    # 5. Construire le tableau
    results = pd.DataFrame(results_rows)

    # 6. Sauvegarder les resultats
    output_csv = OUTPUT_DIR / "baselines_metrics.csv"
    results.to_csv(output_csv, index=False)

    # 7. Afficher les metriques essentielles
    print("\n=== TABLEAU COMPARATIF ===")

    columns = [
    "baseline",
    "n_generated",
    "validity_rate",
    "validity_rate_10_50",
    "uniqueness_rate",
    "length_min",
    "length_mean",
    "length_max",
    "js_divergence",
    "ks_length",
    ]

    print(results[columns].to_string(index=False))

    print(f"\nResultats sauvegardes : {output_csv}")


if __name__ == "__main__":
    main()
