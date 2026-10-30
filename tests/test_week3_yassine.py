import pytest
import torch
import yaml

from peplib.config import load_config
from peplib.seed import set_seed
from peplib.datasets import PeptideDataset, make_loader_s3


def test_configs():
    for name in ("lstm", "transformer"):
        cfg = load_config(f"configs/{name}_base.yaml")
        assert cfg.model["type"] == name
        assert cfg.model["vocab_size"] == 23


def test_unknown_config_key(tmp_path):
    with open("configs/lstm_base.yaml", encoding="utf-8-sig") as f:
        config = yaml.safe_load(f)

    config["model"]["unknown_key"] = 123
    path = tmp_path / "invalid.yaml"
    path.write_text(yaml.safe_dump(config), encoding="utf-8")

    with pytest.raises(ValueError, match="inconnues"):
        load_config(path)


def test_missing_config_key(tmp_path):
    with open("configs/lstm_base.yaml", encoding="utf-8-sig") as f:
        config = yaml.safe_load(f)

    del config["train"]["batch_size"]
    path = tmp_path / "missing.yaml"
    path.write_text(yaml.safe_dump(config), encoding="utf-8")

    with pytest.raises(ValueError, match="manquantes"):
        load_config(path)


def test_seed_reproducibility():
    set_seed(42)
    first = torch.rand(10)
    set_seed(42)
    second = torch.rand(10)
    assert torch.equal(first, second)


def test_loader_padding():
    dataset = PeptideDataset("data/splits/train_v02.csv")
    batch = next(iter(make_loader_s3(dataset, batch_size=4, seed=42)))

    assert batch["input_ids"].shape == batch["labels"].shape
    assert batch["pad_mask"].dtype == torch.bool
    assert torch.equal(
        batch["labels"] == -100,
        ~batch["pad_mask"],
    )


def test_loader_reproducibility():
    dataset = PeptideDataset("data/splits/train_v02.csv")

    first = next(iter(make_loader_s3(
        dataset, batch_size=8, shuffle=True, seed=42
    )))
    second = next(iter(make_loader_s3(
        dataset, batch_size=8, shuffle=True, seed=42
    )))

    for key in first:
        assert torch.equal(first[key], second[key])
