"""T25 - Taux de quasi-doublons a 90 % d'identite (rapport chiffre).

Usage (depuis la racine du repo, venv actif) :
    python scripts/t25_near_duplicates.py
Entree  : data/processed/*.csv  (colonnes : sequence, source, id)
Sortie  : results/t25_near_duplicates.csv  + tableau affiche dans le terminal

Definitions
  - exact      : 1 - (sequences uniques / sequences conservees)           (T22)
  - quasi 90 % : part des sequences UNIQUES ayant au moins un voisin
                 d'identite >= 0.90 (identite = 1 - Lev / max(|a|,|b|)),
                 a l'interieur de la source ("intra") ou dans tout le corpus ("globale").
Depend de peplib.dedup (Sofiene, T23) : build_kmer_index, candidate_pairs, identity.
"""
import glob
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from peplib.dedup import build_kmer_index, candidate_pairs, identity

THRESHOLD = 0.90


def main():
    files = sorted(glob.glob("data/processed/*.csv"))
    if not files:
        sys.exit("Aucun CSV dans data/processed/ : attendre T11-T14 (Yassine).")
    df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    df = df.dropna(subset=["sequence"])

    # --- exact (T22) ---
    stats = {}
    for src, g in df.groupby("source"):
        stats[src] = {"conservees": len(g), "uniques": g["sequence"].nunique()}

    # --- quasi-doublons : sequences uniques (premiere source rencontree) ---
    uniq = df.drop_duplicates("sequence", keep="first").reset_index(drop=True)
    seqs = uniq["sequence"].tolist()
    srcs = uniq["source"].tolist()

    near_any = [False] * len(seqs)    # voisin >= 0.90 n'importe ou
    near_intra = [False] * len(seqs)  # voisin >= 0.90 dans la meme source
    n_pairs = 0
    for i, j in candidate_pairs(build_kmer_index(seqs, 3)):
        a, b = seqs[i], seqs[j]
        if abs(len(a) - len(b)) > (1 - THRESHOLD) * max(len(a), len(b)):
            continue  # identite >= 0.9 impossible : on evite un Levenshtein inutile
        if identity(a, b) >= THRESHOLD:
            n_pairs += 1
            near_any[i] = near_any[j] = True
            if srcs[i] == srcs[j]:
                near_intra[i] = near_intra[j] = True

    rows = []
    for src, s in stats.items():
        idx = [k for k, x in enumerate(srcs) if x == src]
        n = len(idx) or 1
        rows.append({
            "source": src,
            "conservees": s["conservees"],
            "uniques": s["uniques"],
            "redondance_exacte": round(1 - s["uniques"] / s["conservees"], 4),
            "quasi90_intra": round(sum(near_intra[k] for k in idx) / n, 4),
            "quasi90_globale": round(sum(near_any[k] for k in idx) / n, 4),
        })
    out = pd.DataFrame(rows)
    out.to_csv("results/t25_near_duplicates.csv", index=False)
    print(out.to_string(index=False))
    print(f"\nPaires quasi-doublons (>= {THRESHOLD:.0%}) : {n_pairs}")


if __name__ == "__main__":
    main()