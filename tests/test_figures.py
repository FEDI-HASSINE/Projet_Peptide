import math
import os

import matplotlib
import matplotlib.figure
import numpy as np
import pandas as pd
import pytest

import peplib.figures as figmod
from peplib.figures import (
    compute_descriptors,
    load_corpus,
    plot_descriptor_distributions,
    plot_information_curve,
    plot_pwm_heatmap,
    run_figures_pipeline,
)

SEQS = [
    "ACDEFGHIKLM", "MNPQRSTVWYA", "ACDEFGHIKLMN", "LKGHRVDEWIC",
    "AAAAAAGGGGG", "KKKKKKDDDDD", "EEEEEEHHHHHH", "WWWWWWPPPPPP",
    "FILMVYCSTAG", "KRHDEQNCPST",
]


def _corpus(tmp_path):
    proc = tmp_path / "processed"
    proc.mkdir()
    for s in figmod.SOURCES:
        pd.DataFrame(
            {"sequence": SEQS, "source": [s] * len(SEQS),
             "id": [f"{s}_{i}" for i in range(len(SEQS))]}
        ).to_csv(os.path.join(proc, f"{s}.csv"), index=False)
    return proc


def test_load_corpus(tmp_path):
    proc = _corpus(tmp_path)
    corpus = load_corpus(str(proc))
    assert set(corpus.keys()) == {"uniprot", "amp", "bactibase"}
    for df in corpus.values():
        assert list(df.columns) == ["sequence", "source", "id"]
        assert len(df) == len(SEQS)


def test_compute_descriptors():
    df = pd.DataFrame({"sequence": SEQS})
    d = compute_descriptors(df)
    assert list(d.columns) == [
        "charge", "pi", "gravy", "mass", "hydrophobic_moment", "boman"
    ]
    assert not d.isna().to_numpy().any()
    assert np.isfinite(d.to_numpy()).all()
    assert d["pi"].between(0, 14).all()


def test_plot_descriptor_distributions(tmp_path):
    proc = _corpus(tmp_path)
    corpus = load_corpus(str(proc))
    d = {s: compute_descriptors(df) for s, df in corpus.items()}
    out = tmp_path / "desc.png"
    fig = plot_descriptor_distributions(d, str(out))
    assert isinstance(fig, matplotlib.figure.Figure)
    assert out.exists() and out.stat().st_size > 10 * 1024


def test_plot_pwm_heatmap(tmp_path):
    proc = _corpus(tmp_path)
    from peplib.pwm import pwm_matrix, position_information

    corpus = load_corpus(str(proc)).values()
    seqs = [s for df in corpus for s in df["sequence"]]
    pn = pwm_matrix(seqs, terminus="N")
    pc = pwm_matrix(seqs, terminus="C")
    out = tmp_path / "pwm.png"
    fig = plot_pwm_heatmap(pn, pc, str(out))
    assert isinstance(fig, matplotlib.figure.Figure)
    assert out.exists() and out.stat().st_size > 10 * 1024


def test_plot_information_curve(tmp_path):
    proc = _corpus(tmp_path)
    from peplib.pwm import position_information

    corpus = load_corpus(str(proc)).values()
    seqs = [s for df in corpus for s in df["sequence"]]
    info_n = position_information(seqs, terminus="N")
    info_c = position_information(seqs, terminus="C")
    out = tmp_path / "info.png"
    fig = plot_information_curve(info_n, info_c, str(out))
    assert isinstance(fig, matplotlib.figure.Figure)
    assert out.exists() and out.stat().st_size > 10 * 1024


def test_run_figures_pipeline(tmp_path):
    proc = _corpus(tmp_path)
    res = tmp_path / "results"
    paths = run_figures_pipeline(processed_dir=str(proc), results_dir=str(res))
    for name in (
        "descriptors_distribution", "pwm_heatmap", "information_curve"
    ):
        assert os.path.exists(paths[name])
    for csv in ("pwm_N.csv", "pwm_C.csv", "information_N.csv", "information_C.csv"):
        assert os.path.exists(res / csv)


def test_backend_agg():
    assert matplotlib.get_backend().lower() == "agg"
