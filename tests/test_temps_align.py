"""Test léger du module de mesure de temps (A9/T9) : la vraie mesure est
exécutée par le coordinateur en dehors de la suite (grande échelle)."""

from pathlib import Path

import pandas as pd

from peplib.temps_align import run


def test_run_tout_petit_fichier(tmp_path):
    sortie = Path(tmp_path) / "temps.csv"
    df = run(longueurs=(5,), n_pairs=3, out_path=sortie, seed=42)

    assert len(df) == 2          # 1 longueur × 2 méthodes
    assert list(df.columns) == [
        "longueur",
        "methode",
        "ms_par_alignement",
        "paires_tous_contre_tous",
        "extrapolation_heures",
    ]
    assert sortie.exists()
    # lecture de contrôle : le fichier est bien celui écrit
    relit = pd.read_csv(sortie)
    assert len(relit) == 2
