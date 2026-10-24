
"""Tests T37 : construction du corpus v0.2."""

import pandas as pd
import pytest

from peplib.build_corpus_v02 import build_corpus


def make_sources(tmp_path):
    """Crée cinq petites sources pour tester la fusion."""
    data = {
        "uniprot.csv": [
            ("AAAAAAAAAA", "uniprot", "U1"),
            ("CCCCCCCCCC", "uniprot", "U2"),
        ],
        "amp.csv": [
            ("AAAAAAAAAA", "AMP", "A1"),
        ],
        "bactibase.csv": [
            ("DDDDDDDDDD", "Bactibase", "B1"),
        ],
        "amp_apd3.csv": [
            ("AAAAAAAAAA", "APD3", "AP1"),
        ],
        "challenge.csv": [
            ("EEEEEEEEEE", "AMPChallenge", "M1"),
        ],
    }

    for filename, rows in data.items():
        pd.DataFrame(
            rows, columns=["sequence", "source", "id"]
        ).to_csv(tmp_path / filename, index=False)


def test_fusion_five_sources(tmp_path):
    make_sources(tmp_path)
    corpus, stats = build_corpus(tmp_path)

    assert stats["total"] == 6
    assert stats["unique"] == 4
    assert stats["duplicates"] == 2
    assert len(corpus) == 4


def test_no_duplicate_sequences(tmp_path):
    make_sources(tmp_path)
    corpus, _ = build_corpus(tmp_path)

    assert corpus["sequence"].is_unique


def test_provenance_preserved(tmp_path):
    make_sources(tmp_path)
    corpus, _ = build_corpus(tmp_path)

    row = corpus.set_index("sequence").loc["AAAAAAAAAA"]

    assert set(row["sources"].split("|")) == {
        "uniprot", "AMP", "APD3"
    }
    assert set(row["ids"].split("|")) == {
        "U1", "A1", "AP1"
    }
    assert row["n_occurrences"] == 3
    assert row["n_sources"] == 3


def test_output_created(tmp_path):
    make_sources(tmp_path)
    build_corpus(tmp_path)

    output = tmp_path / "corpus_v0.2.csv"
    assert output.exists()

    df = pd.read_csv(output)
    assert len(df) == 4


def test_missing_column_rejected(tmp_path):
    make_sources(tmp_path)

    pd.DataFrame({
        "sequence": ["AAAAAAAAAA"],
        "source": ["uniprot"],
    }).to_csv(tmp_path / "uniprot.csv", index=False)

    with pytest.raises(ValueError, match="Colonnes manquantes"):
        build_corpus(tmp_path)


def test_challenge_metadata_preserved(tmp_path):
    make_sources(tmp_path)

    pd.DataFrame([{
        "sequence": "EEEEEEEEEE",
        "source": "AMPChallenge",
        "id": "M1",
        "len": 10,
        "charge": 2.5,
        "disulfide": 0,
        "activity": "antibacterial",
        "dbs_amont": "APD",
    }]).to_csv(tmp_path / "challenge.csv", index=False)

    corpus, _ = build_corpus(tmp_path)
    row = corpus.set_index("sequence").loc["EEEEEEEEEE"]

    assert row["charge"] == "2.5"
    assert row["activity"] == "antibacterial"
    assert row["dbs_amont"] == "APD"
