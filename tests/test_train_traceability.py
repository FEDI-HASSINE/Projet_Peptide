"""Tests du journal d'entrainement et de sa tracabilite (P0-1)."""
import csv

import pandas as pd

from peplib.train import _append_log, build_model_from_config
from peplib.traceability import config_sha256, run_metadata


def _row(**overrides) -> dict:
    row = {
        "timestamp_utc": "2026-11-23T10:00:00Z",
        "git_sha": "abcdef123456",
        "git_dirty": False,
        "config_sha256": "deadbeef",
        "seed": 42,
        "device": "cuda",
        "torch_version": "2.5.1+cu121",
        "python_version": "3.11.9",
        "gpu_name": "Tesla T4",
        "model": "lstm",
        "epoch": 1,
        "train_loss": 3.1,
        "val_loss": 2.9,
        "val_pp": 18.17,
        "elapsed_s": 41.2,
    }
    row.update(overrides)
    return row


def test_append_log_writes_header_then_rows(tmp_path):
    path = tmp_path / "runs.csv"
    _append_log(path, _row(epoch=1))
    _append_log(path, _row(epoch=2, val_loss=2.8))

    frame = pd.read_csv(path)
    assert list(frame["epoch"]) == [1, 2]
    assert list(frame["git_sha"].unique()) == ["abcdef123456"]


def test_append_log_rotates_on_schema_change(tmp_path):
    path = tmp_path / "runs.csv"
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["model", "seed", "epoch"])
        writer.writeheader()
        writer.writerow({"model": "lstm", "seed": 42, "epoch": 1})

    _append_log(path, _row())

    legacy = list(tmp_path.glob("runs_legacy_*.csv"))
    assert len(legacy) == 1
    assert list(pd.read_csv(legacy[0]).columns) == ["model", "seed", "epoch"]

    frame = pd.read_csv(path)
    assert len(frame) == 1
    assert "git_sha" in frame.columns


def test_run_metadata_matches_config_hash():
    config = {"model": {"type": "lstm"}, "train": {"lr": 0.001}}
    meta = run_metadata(config, seed=7, device="cpu")
    assert meta["config_sha256"] == config_sha256(config)
    assert meta["seed"] == 7


def test_build_model_from_config_lstm():
    config = {
        "model": {
            "type": "lstm",
            "vocab_size": 23,
            "d_embed": 16,
            "d_hidden": 32,
            "n_layers": 1,
            "dropout": 0.0,
        }
    }
    model = build_model_from_config(config)
    assert model.embed_dim == 16
