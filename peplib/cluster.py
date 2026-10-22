from peplib.dedup import find_candidate_pairs, sequence_identity
from pathlib import Path

import random

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
    stats: dict | None = None,
) -> list[int]:
    """T35 : clustering avec filtrage rapide des longueurs.

    Si `stats` est un dict fourni par l'appelant, il est rempli avec les
    compteurs finaux : `candidats`, `eliminees_longueur`, `comparees`,
    `fusions`. Le retour reste inchangé (liste d'identifiants de cluster).
    """

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

        # Filtre de longueur, sans perte : l'identité vérifie
        # id(a,b) = 1 - d/max <= min/max (car d >= | |a|-|b| |),
        # donc min/max < threshold implique identité < threshold :
        # toute paire écartée ici ne pouvait atteindre le seuil.
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

    if stats is not None:
        stats.update(
            candidats=total_candidates,
            eliminees_longueur=filtered_pairs,
            comparees=compared_pairs,
            fusions=merged_pairs,
        )

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


def run_sample_clustering(
    input_path: str = "data/processed/corpus_filtered_all.csv",
    n: int = 200,
    seed: int = 42,
    threshold: float = 0.90,
) -> tuple[list[int], dict]:
    """T35 début — run reproductible sur un échantillon semé du corpus.

    Échantillon : séquences uniques TRIÉES du corpus versionné, tirées
    par `random.Random(seed)`. c'est le point d'entrée à citer dans une
    description de PR : les chiffres ci-dessous sont les SEULS autorisés
    et sont verrouillés par un test (ancre v0.1, corpus 13 082 lignes) :

        candidats 3 326 | éliminées par longueur 2 584 | comparées 742
        fusions 3 | clusters 197

    tout autre chiffre dans une description = non reproductible (rejet).
    """
    df = pd.read_csv(input_path)

    sequences = sorted(df["sequence"].dropna().unique().tolist())

    # PRNG semé (S2245) — usage NON sécuritaire : aucune donnée
    # sensible ; la graine (42 par défaut) sert uniquement à la
    # reproductibilité de l'ancrage (règle projet : graine 42, AGENTS.md).
    # Un Mersenne-Twister semé est ici requis et inoffensif.
    sample = random.Random(seed).sample(sequences, min(n, len(sequences)))  # NOSONAR

    stats: dict = {}
    clusters = cluster_sequences(sample, threshold=threshold, stats=stats)

    return clusters, stats
