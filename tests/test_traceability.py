"""Tests de tracabilite des runs (P0-1)."""
import re
from pathlib import Path

import pandas as pd
import pytest

from peplib.traceability import (
    best_run_table,
    checkpoint_manifest,
    config_sha256,
    model_summary,
    render_week3_provenance,
    run_metadata,
    sha256_file,
    upsert_markdown_block,
    WEEK3_END,
    WEEK3_START,
)


def test_sha256_file_known_value(tmp_path):
    path = tmp_path / "abc.txt"
    path.write_bytes(b"abc")
    assert sha256_file(path) == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"


def test_config_hash_ignores_key_order():
    first = {"model": {"type": "lstm", "d_embed": 128}, "train": {"lr": 0.001}}
    second = {"train": {"lr": 0.001}, "model": {"d_embed": 128, "type": "lstm"}}
    assert config_sha256(first) == config_sha256(second)


def test_config_hash_changes_with_value():
    base = {"train": {"lr": 0.001}}
    other = {"train": {"lr": 0.002}}
    assert config_sha256(base) != config_sha256(other)


def test_checkpoint_manifest(tmp_path):
    (tmp_path / "lstm" / "seed_42").mkdir(parents=True)
    (tmp_path / "lstm" / "seed_42" / "best.pt").write_bytes(b"checkpoint")
    (tmp_path / "transformer" / "seed_43").mkdir(parents=True)
    (tmp_path / "transformer" / "seed_43" / "best.pt").write_bytes(b"autre")

    manifest = checkpoint_manifest(tmp_path)

    assert list(manifest.columns) == ["checkpoint", "sha256", "size_bytes"]
    assert manifest["checkpoint"].tolist() == [
        "lstm/seed_42/best.pt",
        "transformer/seed_43/best.pt",
    ]
    assert manifest.loc[0, "sha256"] == sha256_file(tmp_path / "lstm" / "seed_42" / "best.pt")
    assert manifest.loc[0, "size_bytes"] == len(b"checkpoint")


def test_checkpoint_manifest_empty(tmp_path):
    manifest = checkpoint_manifest(tmp_path)
    assert manifest.empty
    assert list(manifest.columns) == ["checkpoint", "sha256", "size_bytes"]


def test_run_metadata_keys():
    meta = run_metadata({"model": {"type": "lstm"}}, seed=42, device="cpu")
    expected = {
        "timestamp_utc",
        "git_sha",
        "git_dirty",
        "config_sha256",
        "seed",
        "device",
        "torch_version",
        "python_version",
        "gpu_name",
    }
    assert expected <= set(meta)
    assert re.fullmatch(r"[0-9a-f]{12}|unknown", meta["git_sha"])
    assert meta["seed"] == 42
    assert meta["device"] == "cpu"


def test_upsert_block_is_idempotent(tmp_path):
    path = tmp_path / "PROVENANCE.md"
    path.write_text("# Provenance\n\nTexte existant.\n", encoding="utf-8")

    upsert_markdown_block(path, "bloc A")
    first = path.read_text(encoding="utf-8")
    upsert_markdown_block(path, "bloc A")
    second = path.read_text(encoding="utf-8")

    assert first == second
    assert "Texte existant." in second
    assert second.count(WEEK3_START) == 1
    assert second.count(WEEK3_END) == 1


def test_upsert_block_replaces_previous_content(tmp_path):
    path = tmp_path / "PROVENANCE.md"
    path.write_text("# Provenance\n", encoding="utf-8")

    upsert_markdown_block(path, "version 1")
    upsert_markdown_block(path, "version 2")
    text = path.read_text(encoding="utf-8")

    assert "version 2" in text
    assert "version 1" not in text


def _fake_runs() -> pd.DataFrame:
    rows = []
    for model, seed, losses in [
        ("lstm", 42, [3.0, 2.8, 2.7]),
        ("lstm", 43, [3.1, 2.9, 2.75]),
        ("transformer", 42, [3.2, 2.95, 2.9]),
    ]:
        for epoch, loss in enumerate(losses, start=1):
            rows.append(
                {
                    "model": model,
                    "seed": seed,
                    "epoch": epoch,
                    "train_loss": loss - 0.1,
                    "val_loss": loss,
                    "val_pp": 15.0,
                    "git_sha": "abcdef123456",
                    "git_dirty": False,
                    "device": "cuda",
                    "gpu_name": "Tesla T4",
                    "torch_version": "2.5.1+cu121",
                    "python_version": "3.11.9",
                }
            )
    return pd.DataFrame(rows)


def test_best_run_table():
    best = best_run_table(_fake_runs())
    lstm42 = best[(best["model"] == "lstm") & (best["seed"] == 42)].iloc[0]
    assert lstm42["epoch"] == 3
    assert lstm42["val_loss"] == pytest.approx(2.7)


def test_model_summary():
    summary = model_summary(_fake_runs()).set_index("model")
    assert summary.loc["lstm", "n_seeds"] == 2
    assert summary.loc["lstm", "val_pp_mean"] == pytest.approx(15.0)
    assert summary.loc["transformer", "n_seeds"] == 1


def test_render_provenance_contains_gate(tmp_path):
    block = render_week3_provenance(
        runs=_fake_runs(),
        generation=pd.DataFrame(
            [{"model": "lstm", "seed": 42, "n_attempts": 1200, "n_valid": 1000}]
        ),
        baselines={"unigramme": 19.4874, "unigramme_positionnel": 18.6461},
        manifest_csv=tmp_path / "checkpoints_sha256.csv",
        generated_at="2026-11-23T12:00:00Z",
    )
    assert "18.6461" in block
    assert "Gate retenu" in block
    assert "lstm" in block
    assert "Aucun predicteur de structure" in block
