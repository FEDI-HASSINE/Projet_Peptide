"""Tests du script de generation/evaluation Semaine 3 (P0-2)."""
import json
import math

import pandas as pd
import pytest
import torch
import yaml

from scripts.generate_week3 import (
    copy_analysis,
    evaluate_checkpoint,
    load_train_sequences,
    property_summary,
    wilson_interval,
    write_summary_row,
)
from peplib.models.lstm import CausalLSTM
from peplib.seed import set_seed
from peplib.tokenizer import VOCAB


TRAIN = [
    "ACDEFGHIKLM",
    "GGGGGGGGGG",
    "KWKLFKKIEKVGQN",
    "MNPQRSTVWY",
]


def test_wilson_interval_bounds():
    low, high = wilson_interval(0, 10)
    assert low == 0.0
    assert 0.2 < high < 0.35

    low, high = wilson_interval(10, 10)
    assert high == pytest.approx(1.0)
    assert 0.65 < low < 0.8

    for successes, n in [(3, 20), (1, 3), (7, 7)]:
        low, high = wilson_interval(successes, n)
        assert 0.0 <= low <= high <= 1.0


def test_wilson_interval_zero_n():
    assert wilson_interval(0, 0) == (0.0, 0.0)


def test_copy_analysis_exact_and_prefix():
    generated = [
        "ACDEFGHIKLM",   # copie exacte
        "ACDEFGHIKLM",   # deuxieme copie exacte
        "ACDEFGHIKL",    # prefixe du train
        "WWWWWWWWWWW",   # nouvelle
    ]
    result = copy_analysis(generated, TRAIN)

    assert result["n_exact_train"] == 2
    assert result["n_prefix_train"] == 1
    assert result["taux_copie_exacte"] == pytest.approx(0.5)
    assert result["taux_copie_prefixe"] == pytest.approx(0.25)
    assert result["exemples_copie_prefixe"] == ["ACDEFGHIKL"]


def test_copy_analysis_empty():
    result = copy_analysis([], TRAIN)
    assert result["n_exact_train"] == 0
    assert result["taux_copie_exacte"] == 0.0


def test_load_train_sequences_drops_invalid(tmp_path):
    path = tmp_path / "train.csv"
    pd.DataFrame(
        {"sequence": ["ACDEFGHIKLM", "  acdefghiklm  ", "ACDE", "ACDXGHIKLM"]}
    ).to_csv(path, index=False)

    sequences, dropped = load_train_sequences(path)

    assert "ACDEFGHIKLM" in sequences
    assert all(s == s.strip().upper() for s in sequences)
    assert dropped == 2  # trop courte + lettre ambigue X


def test_property_summary_values():
    result = property_summary(["ACDEFGHIKLM", "GGGGGGGGGG"])
    assert set(result) == {
        "charge_moyenne",
        "gravy_moyen",
        "moment_hydrophobe_moyen",
        "pi_moyen",
        "boman_moyen",
    }
    for value in result.values():
        assert isinstance(value, float)
        assert math.isfinite(value)


def test_write_summary_row_deduplicates(tmp_path):
    path = tmp_path / "generation_summary.csv"
    row = {
        "model": "lstm",
        "seed": 42,
        "checkpoint_sha256": "abc",
        "config_sha256": "def",
        "n_requested": 100,
        "temperature": 1.0,
        "top_k": 0,
        "top_p": 1.0,
        "timestamp_utc": "2026-11-23T10:00:00Z",
        "n_valid": 100,
    }
    write_summary_row(path, row)
    write_summary_row(path, {**row, "timestamp_utc": "2026-11-23T11:00:00Z", "n_valid": 99})

    frame = pd.read_csv(path)
    assert len(frame) == 1
    assert frame.loc[0, "n_valid"] == 99
    assert frame.loc[0, "timestamp_utc"] == "2026-11-23T11:00:00Z"


def test_write_summary_row_keeps_distinct_keys(tmp_path):
    path = tmp_path / "generation_summary.csv"
    base = {
        "model": "lstm",
        "seed": 42,
        "checkpoint_sha256": "abc",
        "config_sha256": "def",
        "n_requested": 100,
        "temperature": 1.0,
        "top_k": 0,
        "top_p": 1.0,
    }
    write_summary_row(path, {**base, "n_valid": 100})
    write_summary_row(path, {**base, "seed": 43, "n_valid": 98})

    assert len(pd.read_csv(path)) == 2


def _write_tiny_checkpoint(tmp_path):
    config = {
        "model": {
            "type": "lstm",
            "vocab_size": len(VOCAB),
            "d_embed": 16,
            "d_hidden": 32,
            "n_layers": 1,
            "dropout": 0.0,
        }
    }
    config_path = tmp_path / "lstm_tiny.yaml"
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")

    set_seed(0)
    model = CausalLSTM(
        vocab_size=len(VOCAB), embed_dim=16, hidden_dim=32, n_layers=1, dropout=0.0
    )
    checkpoint_path = tmp_path / "best.pt"
    torch.save(
        {
            "epoch": 4,
            "model_state": model.state_dict(),
            "optimizer_state": {},
            "val_loss": 2.7,
            "val_pp": math.exp(2.7),
            "config": config,
            "seed": 42,
        },
        checkpoint_path,
    )

    train_path = tmp_path / "train.csv"
    pd.DataFrame({"sequence": TRAIN}).to_csv(train_path, index=False)
    return config_path, checkpoint_path, train_path


def test_evaluate_checkpoint_end_to_end(tmp_path):
    config_path, checkpoint_path, train_path = _write_tiny_checkpoint(tmp_path)
    out_dir = tmp_path / "results"

    report = evaluate_checkpoint(
        config_path=config_path,
        checkpoint_path=checkpoint_path,
        seed=42,
        n=3,
        temperature=1.0,
        max_attempts=400,
        n_prox=0,
        train_csv=train_path,
        out_dir=out_dir,
        device="cpu",
    )

    # --- JSON complet et relisible ---
    report_files = list((out_dir / "generation").glob("*.json"))
    assert len(report_files) == 1
    payload = json.loads(report_files[0].read_text(encoding="utf-8"))

    assert payload["metadata"]["checkpoint_sha256"]
    assert payload["metadata"]["config_sha256"]
    assert payload["raw"]["n_raw"] == payload["raw"]["n_attempts"]
    assert payload["raw"]["n_valid"] + payload["raw"]["n_invalid"] == payload["raw"]["n_raw"]
    assert sum(payload["raw"]["rejections"].values()) == payload["raw"]["n_invalid"]
    assert isinstance(payload["sequences"], list)
    assert len(payload["sequences"]) == payload["raw"]["n_valid"]
    assert payload["metadata"]["device"] == "cpu"

    # --- CSV de synthese ---
    summary = pd.read_csv(out_dir / "generation_summary.csv")
    assert len(summary) == 1
    assert summary.loc[0, "n_attempts"] == payload["raw"]["n_attempts"]
    assert summary.loc[0, "taux_copie_exacte"] == pytest.approx(
        payload["copies"]["taux_copie_exacte"]
    )


def test_evaluate_checkpoint_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        evaluate_checkpoint(
            config_path=tmp_path / "absent.yaml",
            checkpoint_path=tmp_path / "absent.pt",
            out_dir=tmp_path,
        )
