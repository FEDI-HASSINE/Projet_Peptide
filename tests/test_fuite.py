"""Tests de la courbe de fuite (A6, livrable D3) : val contre train."""

import pandas as pd

from peplib.fuite import (
    _ecritures,
    compare_splits,
    fractions_F,
    meilleures_couvertures,
)

FRAGMENT = "KTIDWREEKGVVDAGYSQLC"
PEPTIDE_50 = "ADPQSIRNLDEVWAK" + FRAGMENT + "MHTFLKGQESRVDCP"


def test_meilleures_couvertures_fuite_detectee():
    """Val = fragment de son train : le couple (1,0 ; 1,0) est mesuré."""
    df = meilleures_couvertures([FRAGMENT], [PEPTIDE_50])

    assert list(df.columns) == [
        "sequence", "best_identite", "best_couverture", "best_min"
    ]
    assert df["best_identite"][0] == 1.0
    assert df["best_couverture"][0] == 1.0
    assert df["best_min"][0] == 1.0


def test_meilleures_couvertures_paire_intra_val_ignorée():
    """Les paires internes à val (ou à train) ne comptent pas dans la courbe."""
    val = [FRAGMENT, PEPTIDE_50]
    train = ["CGCGCGCGCGCG"]          # sans 3-mer commun avec les deux
    df = meilleures_couvertures(val, train)
    assert (df["best_min"] == 0.0).all()


def test_meilleures_couvertures_sans_candidat():
    val = ["CGCGCGCGCGCG"]            # sans aucun 3-mer commun avec train
    df = meilleures_couvertures(val, ["AAAAAAAAAAAA"])
    assert df["best_min"][0] == 0.0


def test_fractions_f():
    stats = pd.DataFrame(
        {
            "best_identite": [0.95, 0.8, 0.9],
            "best_couverture": [0.95, 0.9, 0.85],
            # best_min = min(identite, couverture) DE LA MEME paire
            "best_min": [0.95, 0.8, 0.85],
        }
    )
    f = fractions_F(stats)
    # (0,95 ; 0,95) passe les deux seuils ; (0,8 ; 0,9) et (0,9 ; 0,85)
    # ne passent que 0,70.
    assert f["F_70"] == 1.0
    assert f["F_90"] == 1 / 3


def test_best_min_invariant_vrai_jeu():
    """best_min est le min de la MEilleure paire, jamais un melange."""
    val = [PEPTIDE_50]
    train = [FRAGMENT, "KTIDWREEKGVVDAGYSQL", "CGCGCGCGCGCG"]
    df = meilleures_couvertures(val, train)
    # le fragment exact domine : meilleure paire = (1,0 ; 1,0)
    assert df["best_min"][0] == min(
        df["best_identite"][0], df["best_couverture"][0]
    ) == 1.0
    # invariant général sur un autre cas à plusieurs candidats
    val2 = ["ACDEFGHIKLNPQRSTVWYAACDEFGHIKLNPQRSTVWYA"]       # 40
    train2 = [
        "ACDEFGHIKLNPQRSTVWYA",                       # moitié exacte
        "ACDEFGHIKLNPQRSTVWYSSACDEFGHIKLNPQRSTVWY",   # quasi total
    ]
    df2 = meilleures_couvertures(val2, train2)
    assert (df2["best_min"] == df2["best_identite"].combine(
        df2["best_couverture"], min
    )).all()


def test_ecritures_fractions_supp_ordonnees(tmp_path):
    """Reprise de calcul : les fractions injectées fusionnent ET l'ordre
    du tableau final est stable (clusters d'abord, contrôle ensuite)."""
    detail = pd.DataFrame(
        {
            "split": ["controle"] * 3,
            "sequence": ["AECDE", "ADEEC", "AAEEE"],
            # best_min = min(identite, couverture) DE LA MÊME paire
            "best_identite": [1.0, 0.8, 0.95],
            "best_couverture": [0.95, 0.9, 0.95],
            "best_min": [0.95, 0.8, 0.95],
        }
    )
    sortie = tmp_path / "fuite.csv"
    sup = [
        {"split": "clusters", "fraction": "F_70", "valeur": 0.6071599045346062},
        {"split": "clusters", "fraction": "F_90", "valeur": 0.0},
    ]
    fractions = _ecritures(
        detail,
        courbes={"controle": detail["best_min"].tolist()},
        sortie_csv=str(sortie),
        figure=str(tmp_path / "courbe.png"),
        fractions_supp=sup,
    )
    # lignes ordonnées clusters d'abord, contrôle ensuite
    assert list(fractions["split"][:2].unique()) == ["clusters"]
    assert list(fractions["split"][2:].unique()) == ["controle"]
    assert fractions.equals(
        fractions.sort_values(["split", "fraction"]).reset_index(drop=True)
    )
    # controle : best_min = (0,95 ; 0,8 ; 0,95) -> F_70 = 3/3 ; F_90 = 2/3
    ligne70 = fractions[
        (fractions["split"] == "controle") & (fractions["fraction"] == "F_70")
    ]["valeur"].iloc[0]
    ligne90 = fractions[
        (fractions["split"] == "controle") & (fractions["fraction"] == "F_90")
    ]["valeur"].iloc[0]
    assert abs(ligne70 - 1.0) < 1e-12
    assert abs(ligne90 - 2 / 3) < 1e-12


def test_compare_splits_toy(tmp_path):
    """Le split par familles est étanche, le contrôle aléatoire fuit."""
    val_clusters = [PEPTIDE_50]
    train_clusters = ["CGCGCGCGCGCG"]              # sans lien avec val

    val_controle = [FRAGMENT]                      # fragment de train !
    train_controle = [PEPTIDE_50, "CGCGCGCGCGCG"]

    fractions, detail = compare_splits(
        val_clusters,
        train_clusters,
        val_controle,
        train_controle,
        sortie_csv=tmp_path / "fuite.csv",
        figure=tmp_path / "courbe_fuite.png",
    )

    assert (tmp_path / "fuite.csv").exists()
    assert (tmp_path / "courbe_fuite.png").exists()

    par_split = fractions.set_index(["split", "fraction"])["valeur"]
    assert par_split[("clusters", "F_90")] == 0.0
    assert par_split[("clusters", "F_70")] == 0.0
    assert par_split[("controle", "F_90")] == 1.0
    assert par_split[("controle", "F_70")] == 1.0
    # le CSV de détail porte une ligne par (séquence, split)
    assert len(detail) == 2
    assert set(detail["split"]) == {"clusters", "controle"}
