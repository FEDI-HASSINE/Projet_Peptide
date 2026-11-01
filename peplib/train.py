"""
peplib/train.py — Boucle d'entraînement générique LSTM / Transformer (T53).

Usage
-----
    # Entraîner le LSTM avec la graine 42 :
    python -m peplib.train --config configs/lstm_base.yaml --seed 42

    # Reprendre depuis le dernier checkpoint :
    python -m peplib.train --config configs/lstm_base.yaml --seed 42 --resume

Fonctionnalités
---------------
- Charge la config YAML (T49).
- Crée les splits si absents (appelle peplib.make_splits).
- Instancie le bon modèle (CausalLSTM ou CausalTransformer).
- Boucle epochs avec :
    * CrossEntropyLoss(ignore_index=-100)
    * Adam + gradient clipping
    * Warmup linéaire optionnel (transformer)
    * Early stopping (patience sur val_loss)
    * Checkpoint best.pt + epoch_N.pt
- Log CSV : results/week3/runs.csv (append)
- Affiche la durée du 1er epoch (décision CPU consignée).
"""
from __future__ import annotations

import argparse
import csv
import math
import time
from pathlib import Path
from typing import Any

import torch
import torch.nn as nn
from torch.utils.data import DataLoader

# Tente d'importer yaml ; repli sur json si absent
try:
    import yaml  # type: ignore
    def _load_config(path: Path) -> dict:
        with open(path) as f:
            return yaml.safe_load(f)
except ImportError:
    import json
    def _load_config(path: Path) -> dict:  # type: ignore[misc]
        # Conversion YAML minimal → dict (supporte uniquement clé: valeur simple)
        cfg: dict = {}
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                if ":" in line:
                    k, _, v = line.partition(":")
                    k = k.strip()
                    v = v.strip()
                    # Conversion de type basique
                    if v.lower() in ("true", "false"):
                        cfg[k] = v.lower() == "true"
                    else:
                        try:
                            cfg[k] = int(v)
                        except ValueError:
                            try:
                                cfg[k] = float(v)
                            except ValueError:
                                cfg[k] = v
        return cfg


from peplib.seed import set_seed
from peplib.make_splits import make_splits
from peplib.datasets import PeptideDataset
from peplib.tokenizer import MAX_LEN


# -----------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------

def _build_model(cfg: dict) -> nn.Module:
    """Instancie le bon modèle selon cfg['model']."""
    model_name = cfg.get("model", "lstm").lower()
    if model_name == "lstm":
        from peplib.models.lstm import CausalLSTM
        return CausalLSTM(
            vocab_size=cfg.get("vocab_size", 23),
            embed_dim=cfg.get("embed_dim", 128),
            hidden_dim=cfg.get("hidden_dim", 256),
            n_layers=cfg.get("n_layers", 2),
            dropout=cfg.get("dropout", 0.2),
        )
    elif model_name == "transformer":
        from peplib.models.transformer import CausalTransformer
        return CausalTransformer(
            vocab_size=cfg.get("vocab_size", 23),
            embed_dim=cfg.get("embed_dim", 128),
            n_heads=cfg.get("n_heads", 4),
            n_layers=cfg.get("n_layers", 4),
            ffn_dim=cfg.get("ffn_dim", 256),
            dropout=cfg.get("dropout", 0.1),
            max_len=cfg.get("max_len", MAX_LEN),
        )
    else:
        raise ValueError(f"Modèle inconnu dans la config : {model_name!r}")


def _make_loaders(cfg: dict) -> tuple[DataLoader, DataLoader]:
    """Crée les DataLoaders train et val."""
    splits_dir = Path(cfg.get("splits_dir", "data/processed/splits"))
    pad_to     = cfg.get("max_len", MAX_LEN)
    batch_size = cfg.get("batch_size", 64)

    train_ds = PeptideDataset(splits_dir / "train.csv", pad_to=pad_to)
    val_ds   = PeptideDataset(splits_dir / "val.csv",   pad_to=pad_to)

    train_loader = DataLoader(
        train_ds, batch_size=batch_size, shuffle=True,
        num_workers=0, drop_last=False,
    )
    val_loader = DataLoader(
        val_ds, batch_size=batch_size, shuffle=False,
        num_workers=0, drop_last=False,
    )
    return train_loader, val_loader


def _warmup_lambda(step: int, warmup_steps: int) -> float:
    """Facteur de multiplication du LR avec warmup linéaire."""
    if warmup_steps <= 0:
        return 1.0
    return min(1.0, step / warmup_steps)


# -----------------------------------------------------------------------
# Epoch train / eval
# -----------------------------------------------------------------------

def _run_epoch_train(
    model: nn.Module,
    loader: DataLoader,
    optimizer: torch.optim.Optimizer,
    criterion: nn.CrossEntropyLoss,
    grad_clip: float,
    scheduler: Any | None,
    device: torch.device,
    model_name: str,
) -> float:
    model.train()
    total_loss = 0.0
    total_tokens = 0

    for batch in loader:
        input_ids = batch["input_ids"].to(device)  # (B, L)
        labels    = batch["labels"].to(device)     # (B, L)

        optimizer.zero_grad()

        if model_name == "lstm":
            logits, _ = model(input_ids)
        else:
            logits = model(input_ids)
        # logits: (B, L, V)

        # Reshape pour CrossEntropyLoss : (B*L, V) vs (B*L,)
        B, L, V = logits.shape
        loss = criterion(logits.view(B * L, V), labels.view(B * L))
        loss.backward()

        if grad_clip > 0:
            nn.utils.clip_grad_norm_(model.parameters(), grad_clip)

        optimizer.step()
        if scheduler is not None:
            scheduler.step()

        # Compter uniquement les tokens non-ignorés pour la loss moyenne
        n_tokens = (labels != -100).sum().item()
        total_loss   += loss.item() * n_tokens
        total_tokens += n_tokens

    return total_loss / max(total_tokens, 1)


@torch.no_grad()
def _run_epoch_val(
    model: nn.Module,
    loader: DataLoader,
    criterion: nn.CrossEntropyLoss,
    device: torch.device,
    model_name: str,
) -> float:
    model.eval()
    total_loss = 0.0
    total_tokens = 0

    for batch in loader:
        input_ids = batch["input_ids"].to(device)
        labels    = batch["labels"].to(device)

        if model_name == "lstm":
            logits, _ = model(input_ids)
        else:
            logits = model(input_ids)

        B, L, V = logits.shape
        loss = criterion(logits.view(B * L, V), labels.view(B * L))

        n_tokens = (labels != -100).sum().item()
        total_loss   += loss.item() * n_tokens
        total_tokens += n_tokens

    return total_loss / max(total_tokens, 1)


# -----------------------------------------------------------------------
# Fonctions de checkpoint
# -----------------------------------------------------------------------

def _save_checkpoint(
    path: Path,
    model: nn.Module,
    optimizer: torch.optim.Optimizer,
    epoch: int,
    val_loss: float,
    val_pp: float,
    cfg: dict,
    seed: int,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    torch.save(
        {
            "epoch": epoch,
            "model_state": model.state_dict(),
            "optimizer_state": optimizer.state_dict(),
            "val_loss": val_loss,
            "val_pp": val_pp,
            "config": cfg,
            "seed": seed,
        },
        path,
    )


def _load_checkpoint(
    path: Path,
    model: nn.Module,
    optimizer: torch.optim.Optimizer,
    device: torch.device,
) -> tuple[int, float]:
    """Charge un checkpoint. Retourne (start_epoch, best_val_loss)."""
    ckpt = torch.load(path, map_location=device)
    model.load_state_dict(ckpt["model_state"])
    optimizer.load_state_dict(ckpt["optimizer_state"])
    return ckpt["epoch"] + 1, ckpt["val_loss"]


# -----------------------------------------------------------------------
# Log CSV
# -----------------------------------------------------------------------

def _append_log(
    csv_path: Path,
    row: dict,
) -> None:
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    write_header = not csv_path.exists()
    with open(csv_path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(row.keys()))
        if write_header:
            writer.writeheader()
        writer.writerow(row)


# -----------------------------------------------------------------------
# Point d'entrée principal
# -----------------------------------------------------------------------

def train_model(
    config_path: str | Path,
    seed: int = 42,
    out_dir: str | Path | None = None,
    resume: bool = False,
) -> dict:
    """Entraîne un modèle selon la config YAML et la graine donnée.

    Parameters
    ----------
    config_path : chemin du fichier YAML de configuration.
    seed        : graine aléatoire (reproductibilité bit-à-bit sur CPU).
    out_dir     : répertoire de checkpoints (écrase cfg si fourni).
    resume      : si True et best.pt existe, reprend l'entraînement.

    Returns
    -------
    dict avec les clés : model, seed, best_val_loss, best_val_pp, best_epoch.
    """
    cfg = _load_config(Path(config_path))
    set_seed(seed)

    model_name = cfg.get("model", "lstm").lower()
    device     = torch.device("cpu")  # CPU uniquement (pas de CUDA détecté)

    # --- Splits ---
    splits_dir = Path(cfg.get("splits_dir", "data/processed/splits"))
    make_splits(splits_dir=splits_dir)

    # --- DataLoaders ---
    train_loader, val_loader = _make_loaders(cfg)

    # --- Modèle ---
    model = _build_model(cfg).to(device)
    n_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print(f"\n[train] Modèle : {model_name.upper()} | seed={seed} | params={n_params:,}")
    print(f"[train] Train : {len(train_loader.dataset):,} séq. | Val : {len(val_loader.dataset):,} séq.")

    # --- Critère + Optimiseur ---
    criterion = nn.CrossEntropyLoss(ignore_index=cfg.get("label_ignore_index", -100))
    optimizer = torch.optim.Adam(model.parameters(), lr=cfg.get("lr", 1e-3))

    # Warmup optionnel (transformer)
    warmup_steps = cfg.get("warmup_steps", 0)
    scheduler = None
    if warmup_steps > 0:
        scheduler = torch.optim.lr_scheduler.LambdaLR(
            optimizer,
            lr_lambda=lambda step: _warmup_lambda(step, warmup_steps),
        )

    # --- Checkpoints ---
    ckpt_base = Path(out_dir) if out_dir else Path(cfg.get("checkpoints_dir", f"checkpoints/{model_name}"))
    ckpt_dir  = ckpt_base / f"seed_{seed}"
    ckpt_dir.mkdir(parents=True, exist_ok=True)
    best_ckpt = ckpt_dir / "best.pt"

    # --- Log CSV ---
    results_dir = Path(cfg.get("results_dir", "results/week3"))
    runs_csv    = results_dir / "runs.csv"

    # --- Reprise ---
    start_epoch  = 0
    best_val_loss = math.inf
    best_epoch    = 0
    no_improve    = 0
    patience      = cfg.get("patience", 5)

    if resume and best_ckpt.exists():
        start_epoch, best_val_loss = _load_checkpoint(best_ckpt, model, optimizer, device)
        print(f"[train] Reprise depuis epoch {start_epoch} (best_val_loss={best_val_loss:.4f})")

    max_epochs = cfg.get("max_epochs", 30)
    grad_clip  = cfg.get("grad_clip", 1.0)

    # --- Boucle epochs ---
    first_epoch_time: float | None = None
    for epoch in range(start_epoch, max_epochs):
        t0 = time.perf_counter()

        train_loss = _run_epoch_train(
            model, train_loader, optimizer, criterion, grad_clip,
            scheduler, device, model_name,
        )
        val_loss = _run_epoch_val(model, val_loader, criterion, device, model_name)
        val_pp   = math.exp(min(val_loss, 20.0))  # cap pour éviter overflow

        elapsed = time.perf_counter() - t0

        if epoch == start_epoch:
            first_epoch_time = elapsed
            device_note = "CPU (pas de CUDA détecté)"
            print(
                f"\n[train] Epoch 0 : {elapsed:.1f} s -- {device_note}\n"
                f"        Decision : entrainement sur {device_note}"
            )

        print(
            f"[train] Epoch {epoch+1:>3}/{max_epochs} | "
            f"train_loss={train_loss:.4f} | "
            f"val_loss={val_loss:.4f} | "
            f"val_pp={val_pp:.2f} | "
            f"{elapsed:.1f}s"
        )

        # --- Log ---
        _append_log(
            runs_csv,
            {
                "model": model_name,
                "seed": seed,
                "epoch": epoch + 1,
                "train_loss": round(train_loss, 6),
                "val_loss": round(val_loss, 6),
                "val_pp": round(val_pp, 4),
                "elapsed_s": round(elapsed, 2),
            },
        )

        # --- Checkpoint epoch ---
        _save_checkpoint(
            ckpt_dir / f"epoch_{epoch+1:03d}.pt",
            model, optimizer, epoch, val_loss, val_pp, cfg, seed,
        )

        # --- Best model ---
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_epoch    = epoch + 1
            no_improve    = 0
            _save_checkpoint(best_ckpt, model, optimizer, epoch, val_loss, val_pp, cfg, seed)
            print(f"[train]   >> Nouveau best model sauve (val_pp={val_pp:.2f})")
        else:
            no_improve += 1
            if no_improve >= patience:
                print(f"[train] Early stopping - {patience} epochs sans amelioration.")
                break

    best_val_pp = math.exp(min(best_val_loss, 20.0))
    print(
        f"\n[train] DONE : {model_name.upper()} seed={seed} | "
        f"best epoch={best_epoch} | val_pp={best_val_pp:.2f}"
    )

    return {
        "model": model_name,
        "seed": seed,
        "best_val_loss": best_val_loss,
        "best_val_pp": best_val_pp,
        "best_epoch": best_epoch,
        "first_epoch_s": first_epoch_time,
        "device": str(device),
    }


# -----------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Entraîne LSTM ou Transformer sur les peptides.")
    parser.add_argument("--config", required=True, help="Chemin du YAML de config")
    parser.add_argument("--seed",   type=int, default=42)
    parser.add_argument("--out-dir", default=None, help="Répertoire de checkpoints (override config)")
    parser.add_argument("--resume", action="store_true", help="Reprendre depuis best.pt si existant")
    args = parser.parse_args()

    result = train_model(args.config, args.seed, args.out_dir, args.resume)
    print("\nRésumé :", result)
