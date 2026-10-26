"""Verrou d'artefact : clusters v0.2 au critère de couverture 0,90 (A5/T42).

Chiffres mesurés le 26/10/2026 par le pipeline déterministe de
``peplib/cluster_v02.py`` (gaps affines o = 11, e = 1, BLOSUM62, préfiltre
multiset borné — ``Docs/A6_preuve_borne_couverture.md``). Toute divergence
signifie que l'artefact n'est plus celui de la PR de mission.
"""

import pytest

import pandas as pd

CLUSTERS_PATH = "data/processed/clusters_v02_90.csv"
T13_PATH = "results/week2/clusters_v02_t13.csv"
N_CORPUS_V02 = 42_575


def _lit_clusters() -> pd.DataFrame:
    return pd.read_csv(CLUSTERS_PATH, dtype=str)


def test_clusters_v02_structure_et_taille():
    df = _lit_clusters()
    assert list(df.columns) == [
        "sequence", "sequence_id", "cluster_id", "sources"
    ]
    assert len(df) == N_CORPUS_V02          # une ligne par séquence unique
    assert df["sequence"].is_unique
    assert df["cluster_id"].notna().all()


def test_clusters_v02_numerotation_compacte():
    """cluster_id : entiers 0..K-1 sans trou (tri des racines, déterminisme)."""
    ids = sorted(set(_lit_clusters()["cluster_id"].astype(int)))
    assert ids == list(range(len(ids)))


def test_clusters_v02_tab13_verrouille():
    t13 = pd.read_csv(T13_PATH)
    assert len(t13) == 1
    assert t13["nb_familles"][0] == 20_673
    assert t13["taille_max"][0] == 4_541
    assert t13["singletons"][0] == 17_592
    assert t13["fraction_plus_grande"][0] == pytest.approx(0.106659)


def test_clusters_v02_grande_famille_depassela_regles():
    """La plus grande famille dépasse 1 % du corpus : la règle « familles
    > 1 % -> train » de ``peplib/splits.py`` sera déclenchée sur cet artefact."""
    df = _lit_clusters()
    plus_grande = df["cluster_id"].value_counts().max()
    assert plus_grande == 4_541
    assert plus_grande > 0.01 * N_CORPUS_V02
