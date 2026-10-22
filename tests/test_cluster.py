from peplib.cluster import UnionFind, cluster_sequences

import pandas as pd
from peplib.cluster import run_clustering
from peplib.dedup import (
    find_candidate_pairs,
    sequence_identity,
)
import random
import pandas as pd


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

    assert len(result) == 4
    assert output_file.exists()
    assert figure_file.exists()

    assert result.loc[0, "cluster_id"] == result.loc[1, "cluster_id"]
    assert result.loc[0, "cluster_id"] == result.loc[2, "cluster_id"]
    assert result.loc[3, "cluster_id"] != result.loc[0, "cluster_id"] 
    
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