"""Tests des splits 90/5/5 par familles et du contrôle aléatoire (T39/A6)."""

import json
from pathlib import Path

import pandas as pd
import pytest

from peplib.splits import (
    assigner_familles,
    ecrire_manifeste,
    run,
    split_aleatoire_control,
)


def test_aucune_famille_a_cheval():
    """Étanchéité : chaque famille entière dans une seule partie."""
    counts = {i: 2 for i in range(100)}    # 100 familles de 2 -> 200
    affect = assigner_familles(counts, total=200)
    assert set(affect.values()) <= {"train", "val", "test"}
    assert all(cid in affect for cid in counts)


def test_proportions_dans_tolerance():
    counts = {i: 2 for i in range(100)}    # 200 séquences
    affect = assigner_familles(counts, total=200)
    par_part = {"train": 0, "val": 0, "test": 0}
    for cid, part in affect.items():
        par_part[part] += counts[cid]
    cibles = {"train": 180, "val": 10, "test": 10}
    for part in par_part:
        assert abs(par_part[part] - cibles[part]) <= 0.005 * 200, part


def test_grande_famille_toujours_en_train():
    """Famille de plus de 1 % du corpus : systématiquement train."""
    # 700 séquences : 1 grande famille de 500 (> 1 % = 7) + 100 familles de 2.
    counts = {0: 500, **{i: 2 for i in range(1, 101)}}
    affect = assigner_familles(counts, total=700)
    assert affect[0] == "train"


def test_determinisme_deux_appels():
    counts = {i: 7 + i for i in range(60)}
    assert assigner_familles(counts) == assigner_familles(counts)


def test_split_controle_memes_tailles():
    rows = pd.DataFrame(
        {
            "sequence": [f"S{i:04d}" for i in range(100)],
            "sources": ["UniProt"] * 100,
            "ids": [f"U{i}" for i in range(100)],
            "cluster_id": [f"C{i % 10}" for i in range(100)],
            "toxicity_status": ["unknown"] * 100,
        }
    )
    parts = split_aleatoire_control(rows, cibles={"train": 90, "val": 5, "test": 5})
    assert len(parts["train"]) == 90
    assert len(parts["val"]) == 5
    assert len(parts["test"]) == 5
    # seed 42 : deux appels identiques
    parts2 = split_aleatoire_control(rows, cibles={"train": 90, "val": 5, "test": 5})
    assert parts["train"]["sequence"].tolist() == parts2["train"]["sequence"].tolist()


def test_manifeste_anti_drift(tmp_path):
    """Le manifeste contient les empreintes et se relit à l'identique."""
    chemin = Path(tmp_path) / "splits_manifest.json"
    manifeste = ecrire_manifeste(
        {"train": 9, "val": 1, "test": 0},
        {"train": {"UniProt": 9}, "val": {"APD3": 1}, "test": {}},
        {"train": SHA_TRAIN, "val": "b" * 64, "test": "c" * 64},
        chemin=chemin,
    )
    relu = json.loads(chemin.read_text(encoding="utf-8"))
    assert relu == manifeste
    assert relu["comptes_par_partie"] == {"train": 9, "val": 1, "test": 0}


def test_run_sur_artefacts_synthetiques(tmp_path):
    """Tube complet sur petits artefacts : étanchéité + colonnes + manifeste.

    2000 « peptides » factices répartis en 1000 familles de 2 : le seuil
    grande-famille vaut 1 % × 2000 = 20 (aucune famille ne le dépasse) et
    les cibles val/test (100 ± 10) tombent juste avec des familles de 2.
    """
    sequences = [f"SEQ{i:05d}" for i in range(2000)]

    # 1000 familles de 2.
    cluster_rows = [
        {"sequence": sequences[i], "sequence_id": f"ID{i}", "cluster_id": f"C{i // 2}", "sources": "UniProt"}
        for i in range(2000)
    ]
    corpus_rows = [
        {"sequence": sequences[i], "sources": "UniProt", "ids": f"ID{i}"}
        for i in range(2000)
    ]
    toxin_rows = [{"sequence": s, "toxicity_status": "not_flagged_by_uniprot"} for s in sequences]

    clusters = pd.DataFrame(cluster_rows)
    corpus = pd.DataFrame(corpus_rows)
    toxin = pd.DataFrame(toxin_rows)
    cpath = Path(tmp_path) / "clusters_c2.csv"
    kpath = Path(tmp_path) / "corpus.csv"
    tpath = Path(tmp_path) / "toxins.csv"
    clusters.to_csv(cpath, index=False)
    corpus.to_csv(kpath, index=False)
    toxin.to_csv(tpath, index=False)

    manifeste, stats = run(
        clusters_path=str(cpath),
        corpus_path=str(kpath),
        toxin_path=str(tpath),
        out_dir=tmp_path / "splits",
        manifest=tmp_path / "splits_manifest.json",
        control_dir=tmp_path / "splits_control",
    )

    assert stats["total"] == 2000
    assert stats["comptes_parties"]["train"] == 1800
    assert stats["comptes_parties"]["val"] == 100
    assert stats["comptes_parties"]["test"] == 100
    # aucune famille à cheval : relecture des 3 CSV + vérification par famille
    familles_par_part = {"train": set(), "val": set(), "test": set()}
    for part in ("train", "val", "test"):
        df = pd.read_csv(tmp_path / "splits" / f"{part}_v02.csv", dtype=str)
        assert list(df.columns) == [
            "sequence", "sources", "ids", "cluster_id", "toxicity_status"
        ]
        familles_par_part[part] = set(df["cluster_id"])
    assert not (familles_par_part["train"] & familles_par_part["val"])
    assert not (familles_par_part["train"] & familles_par_part["test"])
    assert not (familles_par_part["val"] & familles_par_part["test"])
    # total de lignes conservé : 2000 après relecture des trois CSV
    total_rows = (
        len(pd.read_csv(tmp_path / "splits" / "train_v02.csv"))
        + len(pd.read_csv(tmp_path / "splits" / "val_v02.csv"))
        + len(pd.read_csv(tmp_path / "splits" / "test_v02.csv"))
    )
    assert total_rows == 2000
    # le manifeste existe et contient des empreintes SHA-256
    relu = json.loads((tmp_path / "splits_manifest.json").read_text(encoding="utf-8"))
    assert set(relu["sha256"]) == {"train", "val", "test"}
    assert "controle_aleatoire" in relu


SHA_TRAIN = "a" * 64
