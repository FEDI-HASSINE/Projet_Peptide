
"""Tests des vues par source du corpus v0.2."""

import pandas as pd
import pytest

from peplib.source_views_v02 import (
    normalize_sources,
    build_source_views,
)


def test_normalize_single_source():
    assert normalize_sources("UniProt") == {"UniProt"}
    assert normalize_sources("APD3") == {"APD3"}
    assert normalize_sources("AMPChallenge") == {"AMPChallenge"}


def test_normalize_multiple_sources():
    assert normalize_sources("AMPChallenge+APD3") == {
        "AMPChallenge", "APD3"
    }

    assert normalize_sources("UniProt|APD3") == {
        "UniProt", "APD3"
    }


def test_unknown_source():
    with pytest.raises(ValueError, match="Source inconnue"):
        normalize_sources("UnknownDatabase")


def test_source_views(tmp_path):
    corpus = pd.DataFrame({
        "sequence": ["AAAA", "CCCC", "DDDD", "EEEE"],
        "sources": [
            "UniProt",
            "APD3",
            "AMPChallenge",
            "UniProt+AMPChallenge",
        ],
        "ids": ["U1", "A1", "C1", "U2|C2"],
    })

    corpus_path = tmp_path / "corpus.csv"
    corpus.to_csv(corpus_path, index=False)

    results = build_source_views(
        corpus_path=corpus_path,
        output_dir=tmp_path,
    )

    assert results == {
        "UniProt": 2,
        "APD3": 1,
        "AMPChallenge": 2,
    }

    uniprot = pd.read_csv(tmp_path / "view_uniprot_v02.csv")
    apd3 = pd.read_csv(tmp_path / "view_apd3_v02.csv")
    challenge = pd.read_csv(
        tmp_path / "view_ampchallenge_v02.csv"
    )

    assert len(uniprot) == 2
    assert len(apd3) == 1
    assert len(challenge) == 2

    # Une séquence multi-source apparaît dans les deux vues.
    assert "EEEE" in set(uniprot["sequence"])
    assert "EEEE" in set(challenge["sequence"])


def test_duplicate_sequences_rejected(tmp_path):
    corpus = pd.DataFrame({
        "sequence": ["AAAA", "AAAA"],
        "sources": ["UniProt", "APD3"],
    })

    corpus_path = tmp_path / "corpus.csv"
    corpus.to_csv(corpus_path, index=False)

    with pytest.raises(ValueError, match="dupliquées"):
        build_source_views(corpus_path, tmp_path)
