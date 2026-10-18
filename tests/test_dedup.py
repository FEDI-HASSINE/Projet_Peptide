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
