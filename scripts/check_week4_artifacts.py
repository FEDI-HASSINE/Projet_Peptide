"""Contrôles locaux de provenance. Aucune lecture du split test."""
import csv
import json
from pathlib import Path

import pandas as pd
import yaml

from peplib.traceability import checkpoint_manifest, config_sha256, sha256_file


def verify_splits(root):
    manifest = json.loads((root / "splits_manifest.json").read_text(encoding="utf-8"))
    for split in ("train", "val"):
        path = root / "data/splits" / f"{split}_v02.csv"
        with path.open(encoding="utf-8", newline="") as stream:
            size = sum(1 for _ in csv.DictReader(stream))
        if size != manifest["comptes_par_partie"][split] or sha256_file(path) != manifest["sha256"][split]:
            raise ValueError(f"Split incohérent ({split})")


def verify_runs(runs, config_path):
    config = yaml.safe_load(Path(config_path).read_text(encoding="utf-8"))
    required = {"seed", "epoch", "T", "schedule", "nelbo_val", "n_params", "elapsed_s",
                "git_commit", "config_sha256"}
    if required - set(runs.columns) or runs[list(required)].isna().any().any():
        raise ValueError("Métadonnées de run incomplètes")
    if set(runs.seed) != {42, 43, 44} or runs.duplicated(["seed", "epoch"]).any():
        raise ValueError("Graines ou epochs incohérentes")
    if not runs.config_sha256.eq(config_sha256(config)).all():
        raise ValueError("Empreinte de configuration différente")


def verify_checkpoints(manifest_path, checkpoint_dir):
    expected = pd.read_csv(manifest_path)
    actual = checkpoint_manifest(checkpoint_dir)
    if actual.empty or expected.empty:
        raise ValueError("Manifeste ou checkpoints absents")
    columns = list(actual.columns)
    pd.testing.assert_frame_equal(actual, expected[columns].sort_values("checkpoint").reset_index(drop=True))


if __name__ == "__main__":
    verify_splits(Path(__file__).resolve().parents[1])
    print("train et val conformes au manifeste v0.2")
