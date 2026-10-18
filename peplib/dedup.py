import os
from typing import Dict, List, Set, Tuple
import pandas as pd


def exact_dedup_stats(
    sequences: List[str] | pd.Series,
) -> Dict[str, float | int]:
    """T22 — Calcule le nombre total, le nombre unique et le taux de redondance exacte."""
    seq_list = [str(s) for s in sequences]
    total = len(seq_list)
    unique = len(set(seq_list))
    duplicates = total - unique
    redundancy_rate = (duplicates / total) if total > 0 else 0.0
    return {
        "total": total,
        "unique": unique,
        "duplicates": duplicates,
        "redundancy_rate": redundancy_rate,
    }


def build_kmer_index(sequences: List[str], k: int = 3) -> Dict[str, Set[int]]:
    """T23 — Index inversé k-mer -> ensemble d'indices de séquences."""
    index: Dict[str, Set[int]] = {}
    for idx, seq in enumerate(sequences):
        for i in range(len(seq) - k + 1):
            kmer = seq[i : i + k]
            if kmer not in index:
                index[kmer] = set()
            index[kmer].add(idx)
    return index


def find_candidate_pairs(
    sequences: List[str], k: int = 3, min_shared: int = 1
) -> Set[Tuple[int, int]]:
    """T23 — Paires candidates partageant au moins min_shared k-mers."""
    index = build_kmer_index(sequences, k=k)
    pair_counts: Dict[Tuple[int, int], int] = {}
    for indices in index.values():
        sorted_idx = sorted(indices)
        n = len(sorted_idx)
        for a in range(n):
            i = sorted_idx[a]
            for b in range(a + 1, n):
                pair = (i, sorted_idx[b])
                pair_counts[pair] = pair_counts.get(pair, 0) + 1
    return {pair for pair, count in pair_counts.items() if count >= min_shared}


def levenshtein_distance(a: str, b: str) -> int:
    """T23 — Distance de Levenshtein exacte en Python pur."""
    if a == b:
        return 0
    if len(a) < len(b):
        a, b = b, a
    if len(b) == 0:
        return len(a)
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, start=1):
        curr = [i]
        for j, cb in enumerate(b, start=1):
            cost = 0 if ca == cb else 1
            curr.append(min(curr[-1] + 1, prev[j] + 1, prev[j - 1] + cost))
        prev = curr
    return prev[-1]


def sequence_identity(a: str, b: str) -> float:
    """T23 — Identité = 1 - d_Lev(a, b) / max(|a|, |b|)."""
    max_len = max(len(a), len(b))
    if max_len == 0:
        return 1.0
    return 1.0 - (levenshtein_distance(a, b) / max_len)


def compute_near_duplicates_90(
    sequences: List[str], threshold: float = 0.90
) -> Dict[str, float | int]:
    """T25 — Mesure le taux de quasi-doublons à 90% d'identité."""
    unique_seqs = sorted(set(str(s) for s in sequences))
    n_unique = len(unique_seqs)
    if n_unique == 0:
        return {
            "unique_sequences": 0,
            "near_duplicate_pairs": 0,
            "sequences_in_near_dups": 0,
            "near_dup_rate": 0.0,
        }

    candidates = find_candidate_pairs(unique_seqs, k=3, min_shared=3)
    near_dup_seqs: Set[int] = set()
    near_pairs = 0

    for i, j in candidates:
        s1, s2 = unique_seqs[i], unique_seqs[j]
        if min(len(s1), len(s2)) / max(len(s1), len(s2)) < threshold:
            continue
        if sequence_identity(s1, s2) >= threshold:
            near_pairs += 1
            near_dup_seqs.add(i)
            near_dup_seqs.add(j)

    return {
        "unique_sequences": n_unique,
        "near_duplicate_pairs": near_pairs,
        "sequences_in_near_dups": len(near_dup_seqs),
        "near_dup_rate": len(near_dup_seqs) / n_unique,
    }


CLUSTERING_MONTH2_JUSTIFICATION = (
    "Justification du report du clustering complet (union-find) au Mois 2 (T25) :\n"
    "1. En Semaine 1, l'objectif est d'auditer la redondance intrinsèque du corpus brut filtré sans altérer "
    "les distributions par source.\n"
    "2. Le regroupement par composantes connexes (union-find à 90% d'identité) et l'alignement global/local "
    "(Needleman-Wunsch / Smith-Waterman avec BLOSUM62) sont requis au Mois 2 pour construire les splits "
    "train/val/test (90/5/5) par clusters disjoints et garantir l'absence totale de fuite de données (data leakage) "
    "lors de l'entraînement des modèles génératifs."
)


def run_dedup_pipeline(
    processed_dir: str = "data/processed",
) -> Tuple[pd.DataFrame, Dict[str, float | int]]:
    sources = ["uniprot", "amp", "bactibase"]
    dfs = {
        s: pd.read_csv(os.path.join(processed_dir, f"{s}.csv")) for s in sources
    }
    rows = []
    all_seqs = []

    for s, df in dfs.items():
        seqs = df["sequence"].astype(str).tolist()
        all_seqs.extend(seqs)
        st = exact_dedup_stats(seqs)
        rows.append({"scope": s, **st})

    st_global = exact_dedup_stats(all_seqs)
    rows.append({"scope": "global_inter_sources", **st_global})

    df_exact = pd.DataFrame(rows)
    df_exact.to_csv(
        os.path.join(processed_dir, "exact_redundancy_t22.csv"), index=False
    )

    near_stats = compute_near_duplicates_90(all_seqs, threshold=0.90)
    pd.DataFrame([near_stats]).to_csv(
        os.path.join(processed_dir, "near_duplicates_t25.csv"), index=False
    )

    with open(
        os.path.join(processed_dir, "clustering_justification_t25.md"),
        "w",
        encoding="utf-8",
    ) as f:
        f.write(CLUSTERING_MONTH2_JUSTIFICATION)

    return df_exact, near_stats


if __name__ == "__main__":
    df_ex, near_st = run_dedup_pipeline()
    print("=== T22 : Doublons exacts ===")
    print(df_ex.to_string(index=False))
    print("\n=== T25 : Quasi-doublons à 90% ===")
    print(near_st)