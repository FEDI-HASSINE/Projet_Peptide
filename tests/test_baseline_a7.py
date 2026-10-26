"""Tests du stub baseline_a7 (A7 — evaluation finale sur vrais splits).

Le tableau final est DIFFERE (les splits A6 ne sont pas encore sur main,
raison ecrite dans le corps de la PR). Ces tests verrouillent :
  - le refus sonore (FileNotFoundError) tant que les splits manquent ;
  - la production correcte du tableau sur des splits de repechage (toy),
    y compris l'ancre d'ordre 0 : entropie = ln 21 nats/token.
"""
import math
from pathlib import Path

import pytest

from peplib.baseline_a7 import run

SEQS_VALID = [
    "ACDEFGHIKLMNPQRSTVWYAAAA",   # 25 residus
    "GAVLISTWKRQENDCMLFHYPWWWW",  # 25 residus
    "AAAAAAAAAAAAAAAMILPWWWWWW",  # 25 residus
    "VVARWQPGHCSTNDEILKMFYYYYY",  # 25 residus
    "KQWNPVCFLMIYHTSGDERAWAAAA",  # 25 residus
    "FLMHKREQDVWSNCATPIGYWAAAAA", # 26 residus
]


@pytest.fixture
def splits_jouets(tmp_path):
    """Deux pseudo-splits (familles + controle) pour exercer run()."""
    import pandas as pd

    fam, ctrl = tmp_path / "familles", tmp_path / "aleatoire"
    fam.mkdir()
    ctrl.mkdir()
    for dossier in (fam, ctrl):
        pd.DataFrame({"sequence": SEQS_VALID[:4]}).to_csv(
            dossier / "train_v02.csv", index=False, lineterminator="\n"
        )
        pd.DataFrame({"sequence": SEQS_VALID[4:]}).to_csv(
            dossier / "val_v02.csv", index=False, lineterminator="\n"
        )
    return fam, ctrl


def test_refuse_s_sans_splits(tmp_path):
    with pytest.raises(FileNotFoundError):
        run(
            splits_dir=tmp_path / "vide1",
            control_dir=tmp_path / "vide2",
            out_csv=tmp_path / "inutile.csv",
        )


def test_run_sur_splits_jouets(splits_jouets, tmp_path):
    fam, ctrl = splits_jouets
    df = run(
        splits_dir=fam,
        control_dir=ctrl,
        out_csv=tmp_path / "baseline.csv",
    )
    # 4 ordres x 2 splits = 8 lignes, colonnes specifiees
    assert len(df) == 8
    assert list(df.columns) == [
        "split", "ordre", "alpha", "entropie_nats_token", "perplexite",
    ]
    assert set(df["split"]) == {"familles", "aleatoire"}
    assert sorted(df["ordre"].unique()) == [0, 1, 2, 3]
    # grille d'alpha imposee
    assert set(df["alpha"]) <= {0.01, 0.1, 1.0}
    # ancre : ordre 0 uniforme sur 21 cibles -> CE = ln 21 exact
    ordre0 = df[df["ordre"] == 0]
    for ce in ordre0["entropie_nats_token"]:
        assert ce == pytest.approx(math.log(21), rel=1e-12)
    # exp(CE) == perplexite, ligne par ligne
    for _, row in df.iterrows():
        assert row["perplexite"] == pytest.approx(
            math.exp(row["entropie_nats_token"]), rel=1e-12
        )
