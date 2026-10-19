"""
Verification independante du taux de quasi-doublons a 90 % (T25).

Usage :
    python -m peplib.check_near_duplicates

Recalcule le taux SANS utiliser peplib.dedup, puis le compare a
data/processed/near_duplicates_t25.csv.

Methode : index de 3-mers (produit creux), filtre sur : deux sequences
d'identite >= 0.90 partagent au moins D - 3*d_max 3-mers distincts
(D = nb de 3-mers distincts de la plus longue, d_max = floor(0.1 * longueur max)),
puis distance de Levenshtein exacte (bandee) sur les paires retenues.
"""
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp

CORPUS = "data/processed/corpus_filtered_all.csv"
REFERENCE = "data/processed/near_duplicates_t25.csv"
THRESHOLD = 0.90
K = 3


def levenshtein_within(a, b, dmax):
    """True si Levenshtein(a, b) <= dmax (programmation dynamique bandee)."""
    if abs(len(a) - len(b)) > dmax:
        return False
    prev = list(range(len(b) + 1))
    for i, x in enumerate(a, 1):
        cur = [i] + [10**6] * len(b)
        lo, hi = max(1, i - dmax), min(len(b), i + dmax)
        for j in range(lo, hi + 1):
            cur[j] = min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (x != b[j - 1]))
        if min(cur[lo:hi + 1]) > dmax:
            return False
        prev = cur
    return prev[-1] <= dmax


def near_duplicate_stats(sequences):
    seqs = sorted(set(str(s) for s in sequences))
    n = len(seqs)
    lengths = np.array([len(s) for s in seqs])
    distinct = np.zeros(n, dtype=int)
    kmer_id, rows, cols = {}, [], []
    for i, s in enumerate(seqs):
        kmers = {s[p:p + K] for p in range(len(s) - K + 1)}
        distinct[i] = len(kmers)
        for km in kmers:
            rows.append(i)
            cols.append(kmer_id.setdefault(km, len(kmer_id)))
    x = sp.csr_matrix((np.ones(len(rows), dtype=np.int16), (rows, cols)),
                      shape=(n, len(kmer_id)))

    pairs, members = 0, set()
    for start in range(0, n, 500):
        shared = (x[start:start + 500] @ x.T).toarray()
        for r in range(shared.shape[0]):
            i = start + r
            js = np.nonzero(shared[r, i + 1:])[0] + i + 1
            if len(js) == 0:
                continue
            lmax = np.maximum(lengths[i], lengths[js])
            dmax = np.floor((1 - THRESHOLD) * lmax + 1e-9).astype(int)
            d_long = np.where(lengths[i] >= lengths[js], distinct[i], distinct[js])
            ok = (np.abs(lengths[i] - lengths[js]) <= dmax) & (shared[r, js] >= d_long - 3 * dmax)
            for j, dm in zip(js[ok], dmax[ok]):
                if levenshtein_within(seqs[i], seqs[j], int(dm)):
                    pairs += 1
                    members.update((i, j))
    return {"unique_sequences": n, "near_duplicate_pairs": pairs,
            "sequences_in_near_dups": len(members),
            "near_dup_rate": len(members) / n if n else 0.0}


def main():
    mine = near_duplicate_stats(pd.read_csv(CORPUS)["sequence"])
    ref = pd.read_csv(REFERENCE).iloc[0].to_dict()
    ok = True
    for key, value in mine.items():
        same = abs(float(ref[key]) - value) < 1e-9
        ok &= same
        print(f"{key:24s} verif={value:<22} equipe={ref[key]:<22} {'OK' if same else 'DIFFERENT'}")
    print("\nRESULTAT :", "identique a T25 (OK)" if ok else "DIFFERENCE avec T25")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())