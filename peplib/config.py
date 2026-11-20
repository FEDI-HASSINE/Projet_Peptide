"""Chargement et validation stricte des configurations S3."""
from dataclasses import dataclass
from pathlib import Path
import yaml

MODEL_KEYS = {
    "lstm": {"type", "vocab_size", "d_embed", "d_hidden", "n_layers", "dropout", "tie_weights"},
    "transformer": {"type", "vocab_size", "d_model", "n_heads", "n_layers", "d_ff", "dropout"},
}
TRAIN_KEYS = {
    "seed", "batch_size", "max_epochs", "optimizer", "lr",
    "weight_decay", "grad_clip", "early_stopping",
    "checkpoints_dir", "results_dir",
}
DATA_KEYS = {"train", "val", "max_len"}
SAMPLING_KEYS = {"temperature", "top_k", "top_p"}


@dataclass(frozen=True)
class ExperimentConfig:
    model: dict
    train: dict
    data: dict
    sampling_grid: dict


def _check_keys(section, value, expected):
    if not isinstance(value, dict):
        raise ValueError(f"{section} doit etre un dictionnaire")
    missing = expected - value.keys()
    unknown = value.keys() - expected
    if missing or unknown:
        raise ValueError(
            f"{section}: cles manquantes={sorted(missing)}, "
            f"cles inconnues={sorted(unknown)}"
        )


def load_config(path: str | Path) -> ExperimentConfig:
    path = Path(path)
    with path.open("r", encoding="utf-8-sig") as f:
        raw = yaml.safe_load(f)

    _check_keys(
        "racine", raw, {"model", "train", "data", "sampling_grid"}
    )

    model = raw["model"]
    if not isinstance(model, dict):
        raise ValueError("model doit etre un dictionnaire")

    model_type = model.get("type")
    if model_type not in MODEL_KEYS:
        raise ValueError(f"Type de modele inconnu: {model_type}")

    _check_keys("model", model, MODEL_KEYS[model_type])
    # Validation train : cles requises + optionnelles (warmup_steps)
    train_keys = raw["train"]
    missing = TRAIN_KEYS - train_keys.keys()
    if missing:
        raise ValueError(f"train: cles manquantes={sorted(missing)}")
    unknown = train_keys.keys() - (TRAIN_KEYS | {"warmup_steps"})
    if unknown:
        raise ValueError(f"train: cles inconnues={sorted(unknown)}")
    _check_keys("data", raw["data"], DATA_KEYS)
    _check_keys("sampling_grid", raw["sampling_grid"], SAMPLING_KEYS)
    _check_keys(
        "train.early_stopping",
        raw["train"]["early_stopping"],
        {"metric", "patience"},
    )

    if model["vocab_size"] != 23:
        raise ValueError("vocab_size doit etre 23")
    if raw["data"]["max_len"] != 50:
        raise ValueError("max_len doit etre 50")
    if model["n_layers"] < 1:
        raise ValueError("n_layers doit etre >= 1")
    if not 0 <= model["dropout"] < 1:
        raise ValueError("dropout doit etre entre 0 et 1")
    if raw["train"]["batch_size"] < 1:
        raise ValueError("batch_size doit etre positif")
    if raw["train"]["max_epochs"] < 1:
        raise ValueError("max_epochs doit etre positif")
    if raw["train"]["lr"] <= 0:
        raise ValueError("lr doit etre positif")
    if raw["train"]["early_stopping"]["patience"] < 1:
        raise ValueError("patience doit etre positif")
    if model_type == "transformer":
        if model["d_model"] % model["n_heads"] != 0:
            raise ValueError("d_model doit etre divisible par n_heads")
        # warmup_steps est optionnel mais si present doit etre >= 0
        ws = raw["train"].get("warmup_steps", 0)
        if not isinstance(ws, int) or ws < 0:
            raise ValueError("warmup_steps doit etre un entier >= 0")

    return ExperimentConfig(**raw)
