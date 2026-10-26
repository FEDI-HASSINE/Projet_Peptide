"""A9 — Vérification exhaustive des alignements sur le corpus v0.2.

Sofiene vérifie le code ; Youssef mesure (rôles inversés, 1 fois).

Ce script :
1. Charge corpus_v0.2.csv (ou corpus_filtered_all.csv en repli).
2. Tire N paires de séquences uniques avec seed 42.
3. Applique NW, SW, identite_couverture, identite_globale sur chaque paire.
4. Vérifie les invariants suivants pour chaque paire :
   - score_alignment == score NW annoncé (recalcul indépendant)
   - score_alignment == score SW annoncé
   - score SW >= score NW  (local >= global)
   - score SW >= 0
   - 0 <= identite_couverture <= 1
   - 0 <= identite_globale <= 1
   - identite_globale <= identite_couverture
   - symmetrie NW : score(x,y) == score(y,x)
   - symmetrie des identites : cov(x,y) == cov(y,x), glob(x,y) == glob(y,x)
5. Sauvegarde un rapport CSV dans results/align_invariants_a9.csv.
6. Échoue si le nombre de violations > 0.
"""
import csv
import random
from pathlib import Path

import pytest

from peplib.align import (
    alignment_identity,
    identite_couverture,
    identite_globale,
    needleman_wunsch,
    score_alignment,
    smith_waterman,
)

ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = ROOT / "results"
N_PAIRS = 500   # paires testées (augmentable sans modifier les ancres)
SEED = 42


def _load_sequences() -> list[str]:
    """Charge les séquences uniques du corpus, version v0.2 ou repli."""
    for name in ("corpus_v0.2.csv", "corpus_filtered_all.csv"):
        path = ROOT / "data" / "processed" / name
        if path.exists():
            import pandas as pd
            seqs = (
                pd.read_csv(path)["sequence"]
                .dropna()
                .astype(str)
                .unique()
                .tolist()
            )
            return sorted(seqs)
    return []


def _draw_pairs(sequences: list[str], n: int, seed: int) -> list[tuple[str, str]]:
    """Tire n paires (i, j) sans remise, avec i != j."""
    rnd = random.Random(seed)
    pairs = []
    seqs = sequences[:]
    while len(pairs) < n and len(seqs) >= 2:
        i, j = rnd.sample(range(len(seqs)), 2)
        pairs.append((seqs[i], seqs[j]))
    return pairs


# ---------- test principal ----------

def test_exhaustive_alignments_a9():
    sequences = _load_sequences()
    if len(sequences) < 2:
        pytest.skip("Corpus absent ou trop petit pour A9")

    pairs = _draw_pairs(sequences, N_PAIRS, SEED)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    report_path = RESULTS_DIR / "align_invariants_a9.csv"

    fieldnames = [
        "x", "y",
        "nw_score", "sw_score",
        "cov_xy", "glob_xy",
        "cov_yx", "glob_yx",
        "violations",
    ]

    violations_total = 0

    with open(report_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()

        for x, y in pairs:
            viols = []

            # --- NW ---
            nw_xy = needleman_wunsch(x, y)
            nw_yx = needleman_wunsch(y, x)

            # score recalculé = score annoncé
            recalc_nw = score_alignment(nw_xy.aligned_x, nw_xy.aligned_y)
            if recalc_nw != nw_xy.score:
                viols.append(f"nw_score_recalc:{recalc_nw}!={nw_xy.score}")

            # symétrie NW
            if nw_yx.score != nw_xy.score:
                viols.append(f"nw_sym:{nw_yx.score}!={nw_xy.score}")

            # retrait des gaps redonne les séquences originales
            if nw_xy.aligned_x.replace("-", "") != x:
                viols.append("nw_x_strip")
            if nw_xy.aligned_y.replace("-", "") != y:
                viols.append("nw_y_strip")

            # --- SW ---
            sw_xy = smith_waterman(x, y)

            # score recalculé
            recalc_sw = score_alignment(sw_xy.aligned_x, sw_xy.aligned_y)
            if recalc_sw != sw_xy.score:
                viols.append(f"sw_score_recalc:{recalc_sw}!={sw_xy.score}")

            # local >= 0
            if sw_xy.score < 0:
                viols.append(f"sw_neg:{sw_xy.score}")

            # local >= global
            if sw_xy.score < nw_xy.score:
                viols.append(f"sw<nw:{sw_xy.score}<{nw_xy.score}")

            # morceaux contigus
            if sw_xy.aligned_x.replace("-", "") != x[sw_xy.x_start:sw_xy.x_end]:
                viols.append("sw_x_contiguous")
            if sw_xy.aligned_y.replace("-", "") != y[sw_xy.y_start:sw_xy.y_end]:
                viols.append("sw_y_contiguous")

            # --- identités ---
            cov_xy = identite_couverture(nw_xy.aligned_x, nw_xy.aligned_y)
            glob_xy = identite_globale(nw_xy.aligned_x, nw_xy.aligned_y)

            nw_yx_aligned = nw_yx
            cov_yx = identite_couverture(nw_yx_aligned.aligned_x, nw_yx_aligned.aligned_y)
            glob_yx = identite_globale(nw_yx_aligned.aligned_x, nw_yx_aligned.aligned_y)

            # bornes
            if not (0.0 <= cov_xy <= 1.0):
                viols.append(f"cov_bounds:{cov_xy}")
            if not (0.0 <= glob_xy <= 1.0):
                viols.append(f"glob_bounds:{glob_xy}")

            # globale <= couverture
            if glob_xy > cov_xy + 1e-12:
                viols.append(f"glob>cov:{glob_xy}>{cov_xy}")

            # couverture >= alignment_identity
            ai = alignment_identity(nw_xy.aligned_x, nw_xy.aligned_y)
            if cov_xy < ai - 1e-12:
                viols.append(f"cov<ai:{cov_xy}<{ai}")

            # symétrie des identités (NW(x,y) vs NW(y,x))
            if abs(cov_xy - cov_yx) > 1e-12:
                viols.append(f"cov_sym:|{cov_xy}-{cov_yx}|>1e-12")
            if abs(glob_xy - glob_yx) > 1e-12:
                viols.append(f"glob_sym:|{glob_xy}-{glob_yx}|>1e-12")

            violations_total += len(viols)

            writer.writerow({
                "x": x[:30],
                "y": y[:30],
                "nw_score": nw_xy.score,
                "sw_score": sw_xy.score,
                "cov_xy": round(cov_xy, 6),
                "glob_xy": round(glob_xy, 6),
                "cov_yx": round(cov_yx, 6),
                "glob_yx": round(glob_yx, 6),
                "violations": "|".join(viols) if viols else "OK",
            })

    assert violations_total == 0, (
        f"{violations_total} violation(s) détectée(s) sur {len(pairs)} paires. "
        f"Voir {report_path}"
    )
