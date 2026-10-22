from peplib.dedup import find_candidate_pairs, sequence_identity
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


class UnionFind:
    """Structure Union-Find avec compression de chemin et union par rang."""

    def __init__(self, n: int):
        self.parent = list(range(n))
        self.rank = [0] * n

    def find(self, x: int) -> int:
        """Retourne le représentant du groupe de x."""
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, a: int, b: int) -> bool:
        """Fusionne les groupes de a et b."""
        root_a = self.find(a)
        root_b = self.find(b)

        if root_a == root_b:
            return False

        if self.rank[root_a] < self.rank[root_b]:
            root_a, root_b = root_b, root_a

        self.parent[root_b] = root_a

        if self.rank[root_a] == self.rank[root_b]:
            self.rank[root_a] += 1

        return True


def cluster_sequences(
    sequences: list[str],
    threshold: float = 0.90,
) -> list[int]:
    """T35 : clustering avec filtrage rapide des longueurs."""

    if not 0 < threshold <= 1:
        raise ValueError("threshold doit être dans ]0, 1]")

    uf = UnionFind(len(sequences))

    candidates = find_candidate_pairs(
        sequences,
        k=3,
        min_shared=1,
    )

    total_candidates = len(candidates)
    filtered_pairs = 0
    compared_pairs = 0
    merged_pairs = 0

    print(f"Paires candidates : {total_candidates}", flush=True)

    for i, j in sorted(candidates):
        a = sequences[i]
        b = sequences[j]

        # Filtre de longueur : aucune paire admissible
        # ne peut être supprimée par cette condition.
        if min(len(a), len(b)) / max(len(a), len(b)) < threshold:
            filtered_pairs += 1
            continue

        compared_pairs += 1

        identity = sequence_identity(a, b)

        if identity >= threshold:
            if uf.union(i, j):
                merged_pairs += 1

        if compared_pairs % 10000 == 0:
            print(
                f"Comparaisons : {compared_pairs} | "
                f"Fusions : {merged_pairs}",
                flush=True,
            )

    roots = [uf.find(i) for i in range(len(sequences))]

    root_to_cluster = {
        root: idx
        for idx, root in enumerate(sorted(set(roots)))
    }

    print(f"Paires écartées par longueur : {filtered_pairs}")
    print(f"Paires comparées : {compared_pairs}")
    print(f"Fusions : {merged_pairs}")

    return [root_to_cluster[root] for root in roots]

def run_clustering(
    input_path="data/processed/corpus_filtered_all.csv",
    output_path="data/processed/clusters_90.csv",
    figure_path="results/cluster_sizes.png",
):
    """T35 : clustering à 90 % et export des résultats."""

    df = pd.read_csv(input_path)

    required = {"sequence", "source", "id"}
    if not required.issubset(df.columns):
        raise ValueError("Colonnes manquantes dans le corpus")

    if df[list(required)].isna().any().any():
        raise ValueError("Le corpus contient des valeurs manquantes")

    # Travailler sur les séquences uniques pour éviter
    # de comparer inutilement les doublons exacts.
    sequences = sorted(df["sequence"].unique().tolist())

    print(f"Séquences totales : {len(df)}")
    print(f"Séquences uniques : {len(sequences)}")

    cluster_ids = cluster_sequences(
        sequences,
        threshold=0.90,
    )

    # Une séquence identique reçoit toujours le même cluster.
    mapping = dict(zip(sequences, cluster_ids))

    result = pd.DataFrame({
        "sequence_id": df["id"],
        "cluster_id": df["sequence"].map(mapping),
        "source": df["source"],
    })

    Path(output_path).parent.mkdir(
        parents=True, exist_ok=True
    )

    result.to_csv(
        output_path,
        index=False,
        encoding="utf-8",
    )

    # Distribution des tailles en séquences uniques.
    unique_clusters = pd.DataFrame({
        "sequence": sequences,
        "cluster_id": cluster_ids,
    })

    sizes = unique_clusters.groupby("cluster_id").size()

    Path(figure_path).parent.mkdir(
        parents=True, exist_ok=True
    )

    plt.figure(figsize=(9, 5))
    plt.hist(sizes, bins=30, edgecolor="black")
    plt.xlabel("Nombre de séquences uniques par famille")
    plt.ylabel("Nombre de familles")
    plt.title("Distribution des tailles des clusters à 90 %")
    plt.tight_layout()
    plt.savefig(figure_path, dpi=150)
    plt.close()

    print(f"Nombre de clusters : {sizes.size}")
    print(f"Taille du plus grand cluster : {sizes.max()}")
    print(f"CSV généré : {output_path}")
    print(f"Figure générée : {figure_path}")

    return result