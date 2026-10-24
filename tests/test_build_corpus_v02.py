"""Tests T37 : fusion et deduplication exacte du corpus v0.2."""

import pandas as pd
import pytest

from peplib.build_corpus_v02 import SOURCES, build_corpus


def prepare_sources(tmp_path, sequences=None, colonnes=("sequence", "id")):
    sequences = sequences or {
        "uniprot.csv": ["AAAAAAAAAA", "CCCCCCCCCC", "AAAAAAAAAA"],
        "amp.csv": ["AAAAAAAAAA", "GGGGGGGGGG"],
        "bactibase.csv": ["TTTTTTTTTT"],
        "amp_apd3.csv": ["AAAAAAAAAA", "CCCCCCCCCC"],
        "challenge.csv": ["AAAAAAAAAA", "DDDDDDDDDD"],
    }

    for filename in SOURCES:
        valeurs = sequences[filename]
        donnees = {
            "sequence": valeurs,
            "id": [f"{filename}_{i}" for i in range(len(valeurs))],
            "source": [filename] * len(valeurs),
        }
        pd.DataFrame({c: donnees[c] for c in colonnes}).to_csv(
            tmp_path / filename, index=False
        )


def test_fusion_trois_atomes(tmp_path):
    prepare_sources(tmp_path)
    corpus, stats = build_corpus(tmp_path)

    assert stats["total"] == 10
    assert stats["unique"] == 5
    assert stats["duplicates"] == 5
    assert len(corpus) == 5


def test_sources_en_atomes(tmp_path):
    prepare_sources(tmp_path)
    corpus, _ = build_corpus(tmp_path)

    ligne = corpus[corpus["sequence"] == "AAAAAAAAAA"].iloc[0]
    assert ligne["sources"] == "AMPChallenge+APD3+UniProt"
    assert ligne["n_sources"] == 3


def test_amp_et_bactibase_ne_sont_pas_des_sources(tmp_path):
    prepare_sources(tmp_path)
    corpus, _ = build_corpus(tmp_path)

    ligne = corpus[corpus["sequence"] == "GGGGGGGGGG"].iloc[0]
    assert ligne["sources"] == "UniProt"
    assert ligne["n_sources"] == 1
    assert ligne["n_occurrences"] == 1


def test_doublon_interne_a_un_atome(tmp_path):
    prepare_sources(tmp_path)
    corpus, _ = build_corpus(tmp_path)

    ligne = corpus[corpus["sequence"] == "CCCCCCCCCC"].iloc[0]
    assert ligne["n_occurrences"] == 2
    assert ligne["n_sources"] == 2


def test_ids_regroupes(tmp_path):
    prepare_sources(tmp_path)
    corpus, _ = build_corpus(tmp_path)

    ligne = corpus[corpus["sequence"] == "AAAAAAAAAA"].iloc[0]
    assert len(ligne["ids"].split("|")) == 5


def test_longueur_calculee_pour_toutes_les_lignes(tmp_path):
    prepare_sources(tmp_path)
    corpus, _ = build_corpus(tmp_path)

    assert (corpus["longueur"] == corpus["sequence"].str.len()).all()


def test_fichier_ecrit_et_colonnes(tmp_path):
    prepare_sources(tmp_path)
    build_corpus(tmp_path)

    sortie = tmp_path / "corpus_v0.2.csv"
    assert sortie.exists()
    df = pd.read_csv(sortie)

    assert list(df.columns)[:6] == [
        "sequence",
        "sources",
        "ids",
        "n_occurrences",
        "n_sources",
        "longueur",
    ]
    assert len(df) == 5


def test_metadonnees_challenge_conservees(tmp_path):
    prepare_sources(tmp_path)
    challenge = pd.DataFrame(
        {
            "sequence": ["DDDDDDDDDD"],
            "id": ["ML1"],
            "source": ["AMPChallenge"],
            "len": ["10"],
            "charge": ["2.5"],
            "disulfide": ["1"],
            "activity": ["antibacterial"],
            "dbs_amont": ["APD"],
        }
    )
    challenge.to_csv(tmp_path / "challenge.csv", index=False)

    corpus, _ = build_corpus(tmp_path)
    ligne = corpus[corpus["sequence"] == "DDDDDDDDDD"].iloc[0]

    assert ligne["charge"] == "2.5"
    assert ligne["activity"] == "antibacterial"
    assert ligne["dbs_amont"] == "APD"


def test_sequence_vide_rejetee(tmp_path):
    prepare_sources(tmp_path, sequences={
        "uniprot.csv": ["   "],
        "amp.csv": ["AAAAAAAAAA"],
        "bactibase.csv": ["TTTTTTTTTT"],
        "amp_apd3.csv": ["CCCCCCCCCC"],
        "challenge.csv": ["DDDDDDDDDD"],
    })

    with pytest.raises(ValueError, match="vide"):
        build_corpus(tmp_path)


def test_alphabet_non_canonique_rejete(tmp_path):
    prepare_sources(tmp_path, sequences={
        "uniprot.csv": ["ACDXAACDXA"],
        "amp.csv": ["AAAAAAAAAA"],
        "bactibase.csv": ["TTTTTTTTTT"],
        "amp_apd3.csv": ["CCCCCCCCCC"],
        "challenge.csv": ["DDDDDDDDDD"],
    })

    with pytest.raises(ValueError, match="alphabet"):
        build_corpus(tmp_path)


def test_longueur_hors_fenetre_rejetee(tmp_path):
    prepare_sources(tmp_path, sequences={
        "uniprot.csv": ["AAAAA"],
        "amp.csv": ["AAAAAAAAAA"],
        "bactibase.csv": ["TTTTTTTTTT"],
        "amp_apd3.csv": ["CCCCCCCCCC"],
        "challenge.csv": ["DDDDDDDDDD"],
    })

    with pytest.raises(ValueError, match="fenetre"):
        build_corpus(tmp_path)


def test_colonne_manquante_rejetee(tmp_path):
    prepare_sources(tmp_path)
    partiel = pd.DataFrame({"sequence": ["AAAAAAAAAA"]})
    partiel.to_csv(tmp_path / "amp.csv", index=False)

    with pytest.raises(ValueError, match="Colonnes manquantes"):
        build_corpus(tmp_path)
