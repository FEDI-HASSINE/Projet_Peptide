"""Contrôle T72 sur des identifiants jouets. CSV sur stdout."""
import csv
import sys

from scipy.stats import binom
import torch

from peplib.diffusion_data import mask_batch


def measurements():
    steps, trials = 50, 100_000
    generator = torch.Generator().manual_seed(42)
    x0 = torch.zeros((2000, 50), dtype=torch.long)
    for t in range(1, steps + 1):
        p = t / steps
        _, mask, _ = mask_batch(x0, t, generator, T=steps)
        hits = int(mask.sum())
        # Intervalle prédictif binomial bilatéral, Bonferroni sur les 50 pas.
        low, high = binom.interval(1 - 0.05 / steps, trials, p)
        yield dict(t=t, trials=trials, hits=hits, frequency=hits / trials,
                   expected=p, low=int(low), high=int(high), passed=int(low) <= hits <= int(high))


if __name__ == "__main__":
    rows = list(measurements())
    writer = csv.DictWriter(sys.stdout, fieldnames=rows[0], lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    sys.exit(int(not all(row["passed"] for row in rows)))
