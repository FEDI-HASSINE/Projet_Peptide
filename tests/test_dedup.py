import pytest

from peplib.dedup import (
    build_kmer_index,
    compute_near_duplicates_90,
    exact_dedup_stats,
    find_candidate_pairs,
    levenshtein_distance,
    sequence_identity,
)


def test_levenshtein_kitten_sitting():
    """Ancre obligatoire du spec : Levenshtein("kitten", "sitting") == 3."""
    assert levenshtein_distance("kitten", "sitting") == 3


def test_levenshtein_cas_limites():
    assert levenshtein_distance("ABC", "ABC") == 0
    assert levenshtein_distance("", "ABC") == 3
    assert levenshtein_distance("ABC", "") == 3
    assert levenshtein_distance("kitten", "sitting") == levenshtein_distance(
        "sitting", "kitten"
    )
    assert levenshtein_distance("ABC", "AXC") == levenshtein_distance(
        "AXC", "ABC"
    )


def test_exact_dedup_stats():
    st = exact_dedup_stats(["AAA", "aaa", "AAA"])
    assert st["total"] == 3
    assert st["unique"] == 2
    assert st["duplicates"] == 1
    assert st["redundancy_rate"] == pytest.approx(1 / 3)


def test_exact_dedup_stats_vides():
    st = exact_dedup_stats(["", "", ""])
    assert st["total"] == 3
    assert st["unique"] == 1
    assert st["duplicates"] == 2


def test_build_kmer_index():
    index = build_kmer_index(["ABCD", "BCDE"], k=3)
    assert index["BCD"] == {0, 1}


def test_find_candidate_pairs_regles():
    seqs = ["AAAAAAAAAA", "AAAAAAAABC", "WWWWWYYYYY"]

    # Sequences identiques (dupliquees apres dedup interne) : candidates.
    u = sorted(set(seqs))
    pairs = find_candidate_pairs(u, k=3, min_shared=1)
    # 'AAAAAAAAAA' et 'AAAAAAAABC' partagent des 3-mers.
    ia, ib = u.index("AAAAAAAAAA"), u.index("AAAAAAAABC")
    assert (ia, ib) in pairs or (ib, ia) in pairs

    # Sequences sans aucun 3-mer commun : NON candidates.
    assert "CCCCCCCCCC" not in u
    seqs2 = ["AAAAAAAAAA", "CCCCCCCCCC"]
    assert find_candidate_pairs(seqs2, k=3, min_shared=1) == set()

    # k parametrisable.
    for k in (2, 3, 4):
        assert isinstance(find_candidate_pairs(seqs2, k=k), set)


def test_sequence_identity():
    assert sequence_identity("ABC", "ABC") == 1.0
    assert sequence_identity("kitten", "sitting") == pytest.approx(1 - 3 / 7)
    assert sequence_identity("", "A") == 0.0


def test_compute_near_duplicates_90_mini_corpus():
    seqs = [
        "ACDEFGHIKLM",
        "ACDEFGHIKLM",
        "ACDEFGHIKLN",
        "WWWWWWWWWWW",
        "YYYYYYYYYYY",
        "ACDEFGHIKLMN",
    ]
    st = compute_near_duplicates_90(seqs, threshold=0.90)
    assert st["unique_sequences"] == 5
    assert st["near_duplicate_pairs"] == 3
    assert st["sequences_in_near_dups"] == 3
    assert st["near_dup_rate"] == pytest.approx(0.6)

def test_kmer_repetitions():
    # Deux peptides avec 90 % d'identite.
    a = "AAAAAAAAAA"
    b = "AAAAAAAACA"

    assert sequence_identity(a, b) == pytest.approx(0.90)

    # Le seul 3-mer distinct commun est AAA.
    pairs_m1 = find_candidate_pairs([a, b], k=3, min_shared=1)
    pairs_m3 = find_candidate_pairs([a, b], k=3, min_shared=3)

    assert (0, 1) in pairs_m1
    assert (0, 1) not in pairs_m3


def test_near_duplicates_kmers_repetes():
    # Verifier que la valeur par defaut m=1 conserve la paire.
    a = "AAAAAAAAAA"
    b = "AAAAAAAACA"

    stats = compute_near_duplicates_90([a, b])

    assert stats["near_duplicate_pairs"] == 1
    assert stats["sequences_in_near_dups"] == 2


def test_near_duplicates_min_shared_parametrable():
    # Verifier le fonctionnement du parametre min_shared.
    a = "AAAAAAAAAA"
    b = "AAAAAAAACA"

    stats_m1 = compute_near_duplicates_90(
        [a, b], min_shared=1
    )

    stats_m3 = compute_near_duplicates_90(
        [a, b], min_shared=3
    )

    assert stats_m1["near_duplicate_pairs"] == 1
    assert stats_m3["near_duplicate_pairs"] == 0



def test_kmer_rappel_echantillon():
    import random
    import pandas as pd

    # Lire le corpus de peptides v0.1
    df = pd.read_csv("data/processed/corpus_filtered_all.csv")

    # Tirage reproductible de 200 sequences uniques
    sequences = sorted(df["sequence"].dropna().unique().tolist())
    sample = random.Random(42).sample(sequences, 200)

    # Recherche rapide avec les 3-mers
    candidates = find_candidate_pairs(
        sample,
        k=3,
        min_shared=1,
    )

    # Recherche exhaustive sans aucun prefiltre
    vrais_similaires = set()

    for i in range(len(sample)):
        for j in range(i + 1, len(sample)):
            if sequence_identity(sample[i], sample[j]) >= 0.90:
                vrais_similaires.add((i, j))

    # Les vraies paires doivent toutes etre candidates
    manquantes = vrais_similaires - candidates

    print("\nPaires similaires :", len(vrais_similaires))
    print("Paires similaires manquees :", len(manquantes))

    assert len(vrais_similaires) > 0
    assert len(manquantes) == 0
