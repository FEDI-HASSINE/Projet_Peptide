"""Tests T38 : recouvrement exact entre les sources."""

import pandas as pd

from peplib.build_corpus_v02 import ATOMES, SOURCES
from peplib.check_overlaps import (
    compute_overlaps,
    compute_triple,
    parts_par_atome,
    plot_overlaps,
)


def prepare_sources(tmp_path):
    sequences = {
        "uniprot.csv": ["AAAAAAAAAA", "CCCCCCCCCC"],
        "amp.csv": ["AAAAAAAAAA", "GGGGGGGGGG"],
        "bactibase.csv": ["TTTTTTTTTT"],
        "amp_apd3.csv": ["AAAAAAAAAA", "CCCCCCCCCC"],
        "challenge.csv": ["AAAAAAAAAA", "DDDDDDDDDD"],
    }

    for filename in SOURCES:
        valeurs = sequences[filename]
        pd.DataFrame(
            {
                "sequence": valeurs,
                "id": [f"{filename}_{i}" for i in range(len(valeurs))],
            }
        ).to_csv(tmp_path / filename, index=False)


def test_trois_comparaisons(tmp_path):
    prepare_sources(tmp_path)
    resultat = compute_overlaps(tmp_path)

    assert len(resultat) == 3
    assert (tmp_path / "overlap_t38.csv").exists()
    assert set(resultat["source_a"]) | set(resultat["source_b"]) == set(ATOMES)


def test_overlap_et_jaccard(tmp_path):
    prepare_sources(tmp_path)
    resultat = compute_overlaps(tmp_path)

    ligne = resultat[
        (resultat["source_a"] == "UniProt") & (resultat["source_b"] == "APD3")
    ].iloc[0]

    assert ligne["unique_a"] == 4
    assert ligne["unique_b"] == 2
    assert ligne["overlap"] == 2
    assert ligne["jaccard"] == 0.5


def test_jaccard_intermediaire(tmp_path):
    prepare_sources(tmp_path)
    resultat = compute_overlaps(tmp_path)

    ligne = resultat[
        (resultat["source_a"] == "APD3")
        & (resultat["source_b"] == "AMPChallenge")
    ].iloc[0]

    assert ligne["overlap"] == 1
    assert ligne["jaccard"] == 0.333333


def test_triple_recouvrement(tmp_path):
    prepare_sources(tmp_path)
    triple = compute_triple(tmp_path)

    assert triple["sources"] == "UniProt+APD3+AMPChallenge"
    assert triple["overlap"] == 1
    assert triple["union"] == 5
    assert triple["jaccard"] == 0.2
    assert (tmp_path / "overlap_t38_triple.csv").exists()


def test_parts_par_atome(tmp_path):
    prepare_sources(tmp_path)
    parts = parts_par_atome(tmp_path)

    assert parts["UniProt"] == {"unique": 2, "deux": 1, "trois": 1}
    assert parts["APD3"] == {"unique": 0, "deux": 1, "trois": 1}
    assert parts["AMPChallenge"] == {"unique": 1, "deux": 0, "trois": 1}


def test_figure_generee(tmp_path):
    prepare_sources(tmp_path)
    sortie = tmp_path / "figures" / "overlap_t38.png"

    chemin = plot_overlaps(tmp_path, sortie)

    assert chemin.exists()
    assert chemin.stat().st_size > 0
