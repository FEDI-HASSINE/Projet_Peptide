"""Tests de la regeneration des figures (P0-1 : pas de figure inventee)."""
import matplotlib

matplotlib.use("Agg")

import pandas as pd  # noqa: E402
import pytest  # noqa: E402

from scripts.fix_visuals import best_per_run, build_loss_figure, load_runs  # noqa: E402


def _runs_frame() -> pd.DataFrame:
    rows = []
    for model, seed, values in [
        ("lstm", 42, [3.2, 2.9, 2.72]),
        ("lstm", 43, [3.25, 2.95, 2.78]),
        ("transformer", 42, [3.4, 3.1, 2.95]),
    ]:
        for epoch, val_loss in enumerate(values, start=1):
            rows.append(
                {
                    "model": model,
                    "seed": seed,
                    "epoch": epoch,
                    "train_loss": val_loss - 0.15,
                    "val_loss": val_loss,
                    "val_pp": round(2.718281828 ** val_loss, 4),
                }
            )
    return pd.DataFrame(rows)


def test_best_per_run():
    best = best_per_run(_runs_frame())
    lstm = best[(best["model"] == "lstm") & (best["seed"] == 42)].iloc[0]
    assert lstm["epoch"] == 3
    assert lstm["val_loss"] == pytest.approx(2.72)


def test_build_loss_figure_creates_file(tmp_path):
    out_png = tmp_path / "courbes.png"
    best = build_loss_figure(_runs_frame(), out_png)

    assert out_png.exists()
    assert out_png.stat().st_size > 1000
    assert set(best["model"]) == {"lstm", "transformer"}


def test_load_runs_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_runs(tmp_path / "absent.csv")


def test_load_runs_missing_columns(tmp_path):
    path = tmp_path / "runs.csv"
    pd.DataFrame({"model": ["lstm"], "seed": [42]}).to_csv(path, index=False)
    with pytest.raises(ValueError):
        load_runs(path)


def test_load_runs_empty(tmp_path):
    path = tmp_path / "runs.csv"
    pd.DataFrame(
        columns=["model", "seed", "epoch", "train_loss", "val_loss", "val_pp"]
    ).to_csv(path, index=False)
    with pytest.raises(ValueError):
        load_runs(path)
