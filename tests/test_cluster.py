import random
from pathlib import Path

import pandas as pd

from peplib.cluster import (
    UnionFind,
    cluster_sequences,
    levenshtein_bounded,
    max_distance,
    run_clustering,
    run_sample_clustering,
)
from peplib.dedup import (
    find_candidate_pairs,
    levenshtein_distance,
    sequence_identity,
)


def test_union_find_transitivity():
    uf = UnionFind(4)

    uf.union(0, 1)
    uf.union(1, 2)

    assert uf.find(0) == uf.find(2)
    assert uf.find(0) != uf.find(3)


def test_cluster_sequences():
    sequences = [
        "AAAAAAAAAA",
        "AAAAAAAACA",
        "AAAAAAACCA",
        "GGGGGGGGGG",
    ]

    clusters = cluster_sequences(sequences, threshold=0.90)

    assert clusters[0] == clusters[1]
    assert clusters[1] == clusters[2]
    assert clusters[3] != clusters[0]


def test_empty_corpus():
    assert cluster_sequences([]) == []
    
    
def test_clustering_export(tmp_path):
    input_file = tmp_path / "corpus.csv"
    output_file = tmp_path / "clusters.csv"
    figure_file = tmp_path / "clusters.png"

    df = pd.DataFrame({
        "sequence": [
            "AAAAAAAAAA",
            "AAAAAAAACA",
            "AAAAAAAAAA",
            "GGGGGGGGGG",
        ],
        "source": ["uniprot", "amp", "amp", "bactibase"],
        "id": ["p1", "p2", "p3", "p4"],
    })

    df.to_csv(input_file, index=False)

    result = run_clustering(
        input_path=input_file,
        output_path=output_file,
        figure_path=figure_file,
    )

    assert len(result) == 3
    assert output_file.exists()
    assert figure_file.exists()

    assert list(result.columns) == [
        "sequence",
        "sequence_id",
        "cluster_id",
        "source",
    ]
    assert result["sequence"].is_unique
    assert result.loc[0, "cluster_id"] == result.loc[1, "cluster_id"]
    assert result.loc[2, "cluster_id"] != result.loc[0, "cluster_id"]
    
def test_kmer_candidate_coverage():
    sequences = [
        "AAAAAAAAAA",
        "AAAAAAAACA",
        "AAAAAAACCA",
        "GGGGGGGGGG",
        "GGGGGGGGGA",
        "ACDEFGHIKL",
    ]

    threshold = 0.90

    # Méthode exhaustive : toutes les paires
    exact_pairs = set()

    for i in range(len(sequences)):
        for j in range(i + 1, len(sequences)):
            if sequence_identity(
                sequences[i], sequences[j]
            ) >= threshold:
                exact_pairs.add((i, j))

    # Méthode optimisée : index des 3-mers
    candidates = find_candidate_pairs(
        sequences, k=3, min_shared=1
    )

    # Aucun vrai quasi-doublon ne doit être oublié
    assert exact_pairs.issubset(candidates)    
    
def test_kmer_coverage_real_sample():
    df = pd.read_csv(
        "data/processed/corpus_filtered_all.csv"
    )

    sequences = sorted(df["sequence"].dropna().unique())

    rng = random.Random(42)
    sample = rng.sample(sequences, min(200, len(sequences)))

    exact_pairs = set()

    for i in range(len(sample)):
        for j in range(i + 1, len(sample)):
            if sequence_identity(sample[i], sample[j]) >= 0.90:
                exact_pairs.add((i, j))

    candidates = find_candidate_pairs(sample, k=3, min_shared=1)

    missed = exact_pairs - candidates

    print(f"Paires similaires : {len(exact_pairs)}")
    print(f"Paires manquées : {len(missed)}")

    assert not missed


def test_union_find_path_compression():
    # chaîne volontairement non compressée : 0 -> 1 -> 2 (racine)
    uf = UnionFind(5)
    uf.parent[0] = 1
    uf.parent[1] = 2

    assert uf.find(0) == 2

    # la compression de chemin doit raccourcir les deux maillons visités
    assert uf.parent[0] == 2
    assert uf.parent[1] == 2


def test_union_find_union_by_rank():
    uf = UnionFind(4)

    uf.union(0, 1)
    uf.union(2, 3)

    root_a = uf.find(0)
    root_b = uf.find(2)

    # deux arbres de même rang
    assert uf.rank[root_a] == uf.rank[root_b] == 1

    # union à rangs égaux : un arbre s'attache à l'autre, le rang monte à 2
    assert uf.union(1, 3) is True
    assert uf.find(0) == uf.find(2)
    assert uf.rank[uf.find(0)] == 2

    # déjà même groupe : aucune re-fusion, le rang ne bouge plus
    assert uf.union(1, 3) is False
    assert uf.rank[uf.find(0)] == 2


def test_sample_clustering_anchor_v01():
    # Ancre reproductible du PR #12 (T35 début) :
    # échantillon déterministe (graine 42) du corpus v0.1 versionné.
    # Résultats verrouillés (relecture indépendante du 22/10/2026) :
    # candidats 3 326 | éliminées (longueur) 2 584 | comparées 742
    # | fusions 3 | clusters 197.
    # Ce test échouera volontairement si le corpus change : s'applique
    # uniquement à data/processed/corpus_filtered_all.csv (v0.1).
    clusters, stats = run_sample_clustering(
        input_path="data/processed/corpus_filtered_all.csv",
        n=200,
        seed=42,
    )

    assert stats["candidats"] == 3326
    assert stats["eliminees_longueur"] == 2584
    assert stats["comparees"] == 742
    # vérifie en même temps l'arithmétique de la description du PR
    assert stats["comparees"] == stats["candidats"] - stats["eliminees_longueur"]
    assert stats["fusions"] == 3
    assert 200 - len(set(clusters)) == stats["fusions"]
    assert len(set(clusters)) == 197


def test_levenshtein_bounded():
    sequences = [
        "",
        "A",
        "AAAA",
        "AAAC",
        "ACDEFGHIKL",
        "GGGGGGGGGG",
    ]

    for a in sequences:
        for b in sequences:
            for limit in range(6):
                expected = min(
                    levenshtein_distance(a, b),
                    limit + 1,
                )

                assert levenshtein_bounded(a, b, limit) == expected


def test_max_distance_seuil_exact():
    for length in range(10, 51):
        assert max_distance(0.90, length) == length // 10


def test_seuil_identite_exactement_90():
    rng = random.Random(42)

    for length in range(10, 51):
        d = length // 10
        a = "".join(
            rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(length)
        )
        b = list(a)
        for pos in rng.sample(range(length), d):
            b[pos] = "A" if b[pos] != "A" else "C"
        b = "".join(b)

        clusters = cluster_sequences([a, b])
        assert clusters[0] == clusters[1], (length, a, b)


def test_artifact_clusters_v01():
    artifact = pd.read_csv(
        "data/processed/clusters_v01_90.csv",
        dtype=str,
    )
    corpus = pd.read_csv(
        "data/processed/corpus_filtered_all.csv",
        dtype=str,
    )

    assert list(artifact.columns) == [
        "sequence",
        "sequence_id",
        "cluster_id",
        "source",
    ]
    assert not artifact.isna().any().any()
    assert artifact["sequence"].is_unique
    assert len(artifact) == 9021
    assert artifact["cluster_id"].nunique() == 6855

    sizes = artifact.groupby("cluster_id").size()
    assert sizes.max() == 62
    assert sizes.median() == 1

    assert set(artifact["sequence"]) == set(corpus["sequence"].unique())

    raw = Path("data/processed/clusters_v01_90.csv").read_bytes()
    assert b"\r\n" not in raw
    assert raw.endswith(b"\n")
