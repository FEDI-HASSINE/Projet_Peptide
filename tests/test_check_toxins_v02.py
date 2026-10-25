
"""Tests du signalement de toxicité v0.2."""

import pandas as pd
import pytest

from peplib.check_toxins_v02 import check_toxins_v02


@pytest.fixture
def files(tmp_path):
    corpus = pd.DataFrame({
        "sequence": [
            "AAAAAAAAAA",
            "CCCCCCCCCC",
            "DDDDDDDDDD",
            "EEEEEEEEEE",
        ],
        "sources": [
            "UniProt", "UniProt", "APD3", "AMPChallenge"
        ],
        "ids": ["U1", "U2", "AP1", "M1"],
    })

    uniprot = pd.DataFrame({
        "Entry": ["U1", "U2", "U3"],
        "Sequence": [
            "AAAAAAAAAA",
            "CCCCCCCCCC",
            "AAAAAAAAAA",
        ],
        "Protein names": [
            "Conotoxin alpha",
            "Ordinary peptide",
            "Other peptide",
        ],
        "Keywords": ["Toxin", "", ""],
    })

    corpus_path = tmp_path / "corpus.csv"
    uniprot_path = tmp_path / "uniprot.tsv"
    output_path = tmp_path / "flags.csv"

    corpus.to_csv(corpus_path, index=False)
    uniprot.to_csv(uniprot_path, sep="\t", index=False)

    return corpus_path, uniprot_path, output_path


def test_total_and_output(files):
    corpus_path, uniprot_path, output_path = files

    result = check_toxins_v02(
        corpus_path, uniprot_path, output_path
    )

    assert len(result) == 4
    assert output_path.exists()
    assert result["sequence"].is_unique


def test_toxin_detected(files):
    result = check_toxins_v02(*files)
    row = result.set_index("sequence").loc["AAAAAAAAAA"]

    assert bool(row["is_toxin_annotation"])
    assert row["toxicity_status"] == "flagged"


def test_unflagged_is_not_safe(files):
    result = check_toxins_v02(*files)
    row = result.set_index("sequence").loc["CCCCCCCCCC"]

    assert row["toxicity_status"] == "not_flagged_by_uniprot"


def test_unknown_without_annotation(files):
    result = check_toxins_v02(*files)
    row = result.set_index("sequence").loc["DDDDDDDDDD"]

    assert row["toxicity_status"] == "unknown"


def test_duplicate_uniprot_annotation(files):
    result = check_toxins_v02(*files)
    row = result.set_index("sequence").loc["AAAAAAAAAA"]

    # Une annotation toxique reste positive même si une autre
    # annotation de la même séquence ne mentionne pas de toxine.
    assert bool(row["is_toxin_annotation"])


def test_missing_columns_rejected(files):
    corpus_path, uniprot_path, output_path = files

    pd.DataFrame({
        "sequence": ["AAAAAAAAAA"]
    }).to_csv(corpus_path, index=False)

    with pytest.raises(ValueError, match="Colonnes manquantes"):
        check_toxins_v02(corpus_path, uniprot_path, output_path)
