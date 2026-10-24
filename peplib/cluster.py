"""T35 — clustering par composantes connexes à 90 % d'identité.

Les clusters sont des composantes connexes : deux séquences peuvent
appartenir au même cluster sans que leur identité mutuelle atteigne 90 %
(chaîne de quasi-doublons). Mesure sur le corpus v0.1 : 229 des 822
clusters multi-membres (27,9 %) contiennent une paire sous 90 %, la plus
éloignée à 45,5 % (cluster 1892). T39 doit donc traiter le cluster comme
l'unité indivisible du split, jamais la paire.

L'artefact produit est `data/processed/clusters_v01_90.csv` : il porte la
version du corpus qui l'a produit et doit être recalculé sur tout autre
corpus (v0.2 après T37) avant d'alimenter un split.
"""

from fractions import Fraction
from pathlib import Path

import math
import random

import matplotlib.pyplot as plt
import pandas as pd

from peplib.dedup import find_candidate_pairs


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


def levenshtein_bounded(a: str, b: str, max_dist: int) -> int:
    """
    Distance de Levenshtein bornée.

    Retourne la distance exacte si elle est <= max_dist.
    Sinon retourne max_dist + 1.
    """

    if max_dist < 0:
        raise ValueError("max_dist doit être positif ou nul")

    n, m = len(a), len(b)
    inf = max_dist + 1

    # Cas des chaînes vides
    if n == 0:
        return min(m, inf)

    if m == 0:
        return min(n, inf)

    # Différence de longueur trop importante
    if abs(n - m) > max_dist:
        return inf

    if a == b:
        return 0

    prev = [
        j if j <= max_dist else inf
        for j in range(m + 1)
    ]

    for i in range(1, n + 1):
        curr = [inf] * (m + 1)

        if i <= max_dist:
            curr[0] = i

        start = max(1, i - max_dist)
        end = min(m, i + max_dist)

        for j in range(start, end + 1):
            cost = 0 if a[i - 1] == b[j - 1] else 1

            curr[j] = min(
                prev[j] + 1,          # Suppression
                curr[j - 1] + 1,      # Insertion
                prev[j - 1] + cost,   # Substitution
            )

        # Arrêt anticipé si aucune solution admissible
        if min(curr[start:end + 1], default=inf) > max_dist:
            return inf

        prev = curr

    return min(prev[m], inf)


def max_distance(threshold: float, length: int) -> int:
    """Distance maximale admissible pour une identité >= threshold.

    Calcul exact : au seuil 0,90, une longueur 10 autorise exactement une
    édition. Aucun epsilon flottant n'est utilisé.
    """
    if not 0 < threshold <= 1:
        raise ValueError("threshold doit être dans ]0, 1]")
    return math.floor((1 - Fraction(str(threshold))) * length)


def cluster_sequences(
    sequences: list[str],
    threshold: float = 0.90,
    stats: dict | None = None,
) -> list[int]:
    """T35 : clustering Union-Find à 90 % d'identité.

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

        # Éliminer les paires incapables d'atteindre le seuil.
        if min(len(a), len(b)) / max(len(a), len(b)) < threshold:
            filtered_pairs += 1
            continue

        compared_pairs += 1

        # Distance maximale autorisée, calculée exactement.
        max_dist = max_distance(threshold, max(len(a), len(b)))

        # Calcul optimisé de Levenshtein.
        distance = levenshtein_bounded(a, b, max_dist)

        # Fusionner les familles si identité >= 90 %.
        if distance <= max_dist:
            if uf.union(i, j):
                merged_pairs += 1

        # Afficher la progression.
        if compared_pairs % 10000 == 0:
            print(
                f"Comparaisons : {compared_pairs} | "
                f"Fusions : {merged_pairs}",
                flush=True,
            )

    # Représentant de chaque séquence.
    roots = [uf.find(i) for i in range(len(sequences))]

    # Attribuer des numéros de clusters.
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
    output_path="data/processed/clusters_v01_90.csv",
    figure_path="results/cluster_sizes_v01.png",
):
    """T35 : clustering à 90 % et export des résultats.

    Le CSV produit contient une ligne par séquence unique
    (`sequence, sequence_id, cluster_id, source`) : il est autosuffisant
    et ne dépend ni d'une jointure ni de l'ordre des lignes du corpus.
    """

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

    # Une ligne par séquence unique : l'artefact est autosuffisant.
    # `drop_duplicates` garde la première occurrence du corpus (identifiant
    # et source de provenance) ; la séquence porte le cluster.
    meta = df.drop_duplicates("sequence", keep="first")

    result = pd.DataFrame({
        "sequence": meta["sequence"].to_numpy(),
        "sequence_id": meta["id"].to_numpy(),
        "cluster_id": meta["sequence"].map(mapping).to_numpy(),
        "source": meta["source"].to_numpy(),
    })

    Path(output_path).parent.mkdir(
        parents=True, exist_ok=True
    )

    result.to_csv(
        output_path,
        index=False,
        encoding="utf-8",
        lineterminator="\n",
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


if __name__ == "__main__":
    run_clustering()

