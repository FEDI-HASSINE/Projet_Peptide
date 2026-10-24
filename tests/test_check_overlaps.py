
"""Tests T38 : recouvrement exact entre sources."""

import pandas as pd

from peplib.check_overlaps import compute_overlaps
from peplib.build_corpus_v02 import SOURCES


def prepare_sources(tmp_path):
    sequences = {
        "uniprot": ["AAAAAAAAAA", "CCCCCCCCCC"],
        "amp": ["AAAAAAAAAA", "GGGGGGGGGG"],
        "bactibase": ["TTTTTTTTTT"],
        "apd3": ["AAAAAAAAAA", "CCCCCCCCCC"],
        "challenge": ["AAAAAAAAAA", "DDDDDDDDDD"],
    }

    for source, filename in SOURCES.items():
        pd.DataFrame({
            "sequence": sequences[source],
            "source": source,
            "id": [f"{source}_{i}" for i in range(len(sequences[source]))],
        }).to_csv(tmp_path / filename, index=False)


def test_ten_pairwise_comparisons(tmp_path):
    prepare_sources(tmp_path)
    result = compute_overlaps(tmp_path)

    assert len(result) == 10
    assert (tmp_path / "overlap_t38.csv").exists()


def test_overlap_and_jaccard(tmp_path):
    prepare_sources(tmp_path)
    result = compute_overlaps(tmp_path)

    row = result[
        (result["source_a"] == "uniprot")
        & (result["source_b"] == "apd3")
    ].iloc[0]

    assert row["overlap"] == 2
    assert row["jaccard"] == 1.0


def test_no_overlap(tmp_path):
    prepare_sources(tmp_path)
    result = compute_overlaps(tmp_path)

    row = result[
        (result["source_a"] == "uniprot")
        & (result["source_b"] == "bactibase")
    ].iloc[0]

    assert row["overlap"] == 0
    assert row["jaccard"] == 0.0


def test_unique_sequences_only(tmp_path):
    prepare_sources(tmp_path)

    path = tmp_path / "amp.csv"
    df = pd.read_csv(path)
    df = pd.concat([df, df.iloc[[0]]], ignore_index=True)
    df.to_csv(path, index=False)

    result = compute_overlaps(tmp_path)

    row = result[
        (result["source_a"] == "uniprot")
        & (result["source_b"] == "amp")
    ].iloc[0]

    assert row["unique_b"] == 2
    assert row["overlap"] == 1
