"""
peplib/train.py -- Boucle d'entrainement generique LSTM / Transformer (T53).

Corrections v2 (suite a revision PR #25) :
  - Schema de config imbrique (harmonise avec PR #24, Yassine).
  - Splits verrouilles : data/splits/{train,val}_v02.csv (jamais make_splits).
  - Optimiseur : AdamW (weight_decay=0.01) conforme a la PR #24.
  - DataLoader avec torch.Generator seede (reproductibilite bit-a-bit).
  - Aucun caractere Unicode dans les print() (compatibilite cp1252 Windows).

Usage
-----
    python -m peplib.train --config configs/lstm_base.yaml --seed 42
    python -m peplib.train --config configs/transformer_base.yaml --seed 42 --resume
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
import yaml

from peplib.seed import set_seed
from peplib.datasets import PeptideDataset
from peplib.tokenizer import MAX_LEN


# -----------------------------------------------------------------------
# Chargement de la config (schema imbrique PR #24)
# -----------------------------------------------------------------------

def _load_config(path: Path) -> dict:
    with open(path, encoding="utf-8") as f:
        return yaml.safe_load(f)


def _cfg(cfg: dict, *keys, default=None):
    """Acces securise a une cle imbriquee : _cfg(cfg, 'train', 'lr')."""
    obj = cfg
    for k in keys:
        if not isinstance(obj, dict) or k not in obj:
            return default
        obj = obj[k]
    return obj


# -----------------------------------------------------------------------
# Construction du modele
# -----------------------------------------------------------------------

def _build_model(cfg: dict) -> nn.Module:
    mtype = _cfg(cfg, "model", "type", default="lstm").lower()
    if mtype == "lstm":
        from peplib.models.lstm import CausalLSTM
        return CausalLSTM(
            vocab_size=_cfg(cfg, "model", "vocab_size", default=23),
            embed_dim=_cfg(cfg, "model", "d_embed", default=128),
            hidden_dim=_cfg(cfg, "model", "d_hidden", default=256),
            n_layers=_cfg(cfg, "model", "n_layers", default=2),
            dropout=_cfg(cfg, "model", "dropout", default=0.2),
        )
    elif mtype == "transformer":
        from peplib.models.transformer import CausalTransformer
        return CausalTransformer(
            vocab_size=_cfg(cfg, "model", "vocab_size", default=23),
            embed_dim=_cfg(cfg, "model", "d_model", default=128),
            n_heads=_cfg(cfg, "model", "n_heads", default=4),
            n_layers=_cfg(cfg, "model", "n_layers", default=4),
            ffn_dim=_cfg(cfg, "model", "d_ff", default=512),
            dropout=_cfg(cfg, "model", "dropout", default=0.1),
            max_len=MAX_LEN,
        )
    raise ValueError(f"model.type inconnu : {mtype!r}")


# -----------------------------------------------------------------------
# DataLoaders avec generateur seede (reproductibilite bit-a-bit)
# -----------------------------------------------------------------------

def _make_loaders(cfg: dict, seed: int) -> tuple[DataLoader, DataLoader]:
    """Cree train/val DataLoaders a partir des splits verrouilles v02."""
    train_path = Path(_cfg(cfg, "data", "train", default="data/splits/train_v02.csv"))
    val_path   = Path(_cfg(cfg, "data", "val",   default="data/splits/val_v02.csv"))
    max_len    = _cfg(cfg, "data", "max_len", default=50) + 2  # +BOS+EOS
    batch_size = _cfg(cfg, "train", "batch_size", default=128)

    if not train_path.exists():
        raise FileNotFoundError(
            f"Split train introuvable : {train_path}\n"
            "Verifier que les artefacts data/splits/*_v02.csv sont bien presents (PR #20)."
        )
    if not val_path.exists():
        raise FileNotFoundError(f"Split val introuvable : {val_path}")

    train_ds = PeptideDataset(train_path, pad_to=max_len)
    val_ds   = PeptideDataset(val_path,   pad_to=max_len)

    # Generateur seede pour le DataLoader train (shuffle reproductible)
    g = torch.Generator()
    g.manual_seed(seed)

    train_loader = DataLoader(
        train_ds, batch_size=batch_size, shuffle=True,
        num_workers=0, drop_last=False, generator=g,
    )
    val_loader = DataLoader(
        val_ds, batch_size=batch_size, shuffle=False,
        num_workers=0, drop_last=False,
    )
    return train_loader, val_loader


# -----------------------------------------------------------------------
# Warmup lineaire (optionnel, transformer)
# -----------------------------------------------------------------------

def _warmup_lambda(step: int, warmup_steps: int) -> float:
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
    mtype: str,
) -> float:
    model.train()
    total_loss, total_tokens = 0.0, 0

    for batch in loader:
        input_ids = batch["input_ids"].to(device)
        labels    = batch["labels"].to(device)

        optimizer.zero_grad()
        if mtype == "lstm":
            logits, _ = model(input_ids)
        else:
            logits = model(input_ids)

        B, L, V = logits.shape
        loss = criterion(logits.view(B * L, V), labels.view(B * L))
        loss.backward()

        if grad_clip > 0:
            nn.utils.clip_grad_norm_(model.parameters(), grad_clip)
        optimizer.step()
        if scheduler is not None:
            scheduler.step()

        n_tok = (labels != -100).sum().item()
        total_loss   += loss.item() * n_tok
        total_tokens += n_tok

    return total_loss / max(total_tokens, 1)


@torch.no_grad()
def _run_epoch_val(
    model: nn.Module,
    loader: DataLoader,
    criterion: nn.CrossEntropyLoss,
    device: torch.device,
    mtype: str,
) -> float:
    model.eval()
    total_loss, total_tokens = 0.0, 0

    for batch in loader:
        input_ids = batch["input_ids"].to(device)
        labels    = batch["labels"].to(device)

        if mtype == "lstm":
            logits, _ = model(input_ids)
        else:
            logits = model(input_ids)

        B, L, V = logits.shape
        loss = criterion(logits.view(B * L, V), labels.view(B * L))

        n_tok = (labels != -100).sum().item()
        total_loss   += loss.item() * n_tok
        total_tokens += n_tok

    return total_loss / max(total_tokens, 1)


# -----------------------------------------------------------------------
# Checkpoint
# -----------------------------------------------------------------------

def _save_ckpt(path: Path, model, optimizer, epoch, val_loss, val_pp, cfg, seed):
    path.parent.mkdir(parents=True, exist_ok=True)
    torch.save({
        "epoch": epoch,
        "model_state": model.state_dict(),
        "optimizer_state": optimizer.state_dict(),
        "val_loss": val_loss,
        "val_pp": val_pp,
        "config": cfg,
        "seed": seed,
    }, path)


def _load_ckpt(path: Path, model, optimizer, device) -> tuple[int, float]:
    ckpt = torch.load(path, map_location=device)
    model.load_state_dict(ckpt["model_state"])
    optimizer.load_state_dict(ckpt["optimizer_state"])
    return ckpt["epoch"] + 1, ckpt["val_loss"]


# -----------------------------------------------------------------------
# Log CSV
# -----------------------------------------------------------------------

def _append_log(csv_path: Path, row: dict) -> None:
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    write_header = not csv_path.exists()
    with open(csv_path, "a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(row.keys()))
        if write_header:
            writer.writeheader()
        writer.writerow(row)


# -----------------------------------------------------------------------
# Point d'entree
# -----------------------------------------------------------------------

def train_model(
    config_path: str | Path,
    seed: int = 42,
    resume: bool = False,
) -> dict:
    """Entraine un modele selon la config YAML et la graine donnee.

    Parameters
    ----------
    config_path : chemin du YAML de configuration (schema imbrique PR #24).
    seed        : graine pour la reproductibilite bit-a-bit sur CPU.
    resume      : si True et best.pt existe, reprend depuis ce checkpoint.

    Returns
    -------
    dict : model, seed, best_val_loss, best_val_pp, best_epoch.
    """
    cfg = _load_config(Path(config_path))
    set_seed(seed)

    mtype  = _cfg(cfg, "model", "type", default="lstm").lower()
    device = torch.device("cpu")

    # --- Splits verrouilles (aucune creation automatique) ---
    train_loader, val_loader = _make_loaders(cfg, seed)

    # --- Modele ---
    model  = _build_model(cfg).to(device)
    n_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print(f"\n[train] {mtype.upper()} | seed={seed} | params={n_params:,} | device=CPU")
    print(f"[train] train={len(train_loader.dataset):,} | val={len(val_loader.dataset):,}")

    # --- Critere ---
    criterion = nn.CrossEntropyLoss(ignore_index=-100)

    # --- Optimiseur : AdamW (aligne PR #24) ---
    lr           = _cfg(cfg, "train", "lr", default=1e-3)
    weight_decay = _cfg(cfg, "train", "weight_decay", default=0.01)
    optimizer = torch.optim.AdamW(
        model.parameters(), lr=lr, weight_decay=weight_decay
    )

    # Warmup lineaire optionnel
    warmup_steps = _cfg(cfg, "train", "warmup_steps", default=0)
    scheduler    = None
    if warmup_steps > 0:
        scheduler = torch.optim.lr_scheduler.LambdaLR(
            optimizer,
            lr_lambda=lambda s: _warmup_lambda(s, warmup_steps),
        )

    # --- Chemins (dans train.* pour compatibilite schema PR #24) ---
    ckpt_base = Path(_cfg(cfg, "train", "checkpoints_dir", default=f"checkpoints/{mtype}"))
    ckpt_dir  = ckpt_base / f"seed_{seed}"
    ckpt_dir.mkdir(parents=True, exist_ok=True)
    best_ckpt = ckpt_dir / "best.pt"

    results_dir = Path(_cfg(cfg, "train", "results_dir", default="results/week3"))
    runs_csv    = results_dir / "runs.csv"

    # --- Reprise ---
    start_epoch   = 0
    best_val_loss = math.inf
    best_epoch    = 0
    no_improve    = 0
    patience      = _cfg(cfg, "train", "early_stopping", "patience", default=5)
    max_epochs    = _cfg(cfg, "train", "max_epochs", default=30)
    grad_clip     = _cfg(cfg, "train", "grad_clip", default=1.0)

    if resume and best_ckpt.exists():
        start_epoch, best_val_loss = _load_ckpt(best_ckpt, model, optimizer, device)
        print(f"[train] Reprise depuis epoch {start_epoch} | best_val_loss={best_val_loss:.4f}")

    # --- Boucle d'entrainement ---
    first_epoch_time: float | None = None

    for epoch in range(start_epoch, max_epochs):
        t0 = time.perf_counter()

        train_loss = _run_epoch_train(
            model, train_loader, optimizer, criterion, grad_clip,
            scheduler, device, mtype,
        )
        val_loss = _run_epoch_val(model, val_loader, criterion, device, mtype)
        val_pp   = math.exp(min(val_loss, 20.0))
        elapsed  = time.perf_counter() - t0

        if epoch == start_epoch:
            first_epoch_time = elapsed
            print(f"[train] Epoch 0 : {elapsed:.1f}s -- CPU (pas de CUDA)")

        print(
            f"[train] Epoch {epoch+1:>3}/{max_epochs} | "
            f"train_loss={train_loss:.4f} | val_loss={val_loss:.4f} | "
            f"val_pp={val_pp:.2f} | {elapsed:.1f}s"
        )

        # Log CSV
        _append_log(runs_csv, {
            "model":      mtype,
            "seed":       seed,
            "epoch":      epoch + 1,
            "train_loss": round(train_loss, 6),
            "val_loss":   round(val_loss, 6),
            "val_pp":     round(val_pp, 4),
            "elapsed_s":  round(elapsed, 2),
        })

        # Checkpoint epoch
        _save_ckpt(
            ckpt_dir / f"epoch_{epoch+1:03d}.pt",
            model, optimizer, epoch, val_loss, val_pp, cfg, seed,
        )

        # Best model
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_epoch    = epoch + 1
            no_improve    = 0
            _save_ckpt(best_ckpt, model, optimizer, epoch, val_loss, val_pp, cfg, seed)
            print(f"[train]   >> best model sauve (val_pp={val_pp:.2f})")
        else:
            no_improve += 1
            if no_improve >= patience:
                print(f"[train] Early stopping -- {patience} epochs sans amelioration.")
                break

    best_val_pp = math.exp(min(best_val_loss, 20.0))
    print(
        f"\n[train] DONE : {mtype.upper()} seed={seed} | "
        f"best epoch={best_epoch} | val_pp={best_val_pp:.2f}"
    )
    return {
        "model":          mtype,
        "seed":           seed,
        "best_val_loss":  best_val_loss,
        "best_val_pp":    best_val_pp,
        "best_epoch":     best_epoch,
        "first_epoch_s":  first_epoch_time,
        "device":         str(device),
    }


# -----------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Entraine LSTM ou Transformer.")
    parser.add_argument("--config", required=True)
    parser.add_argument("--seed",   type=int, default=42)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    result = train_model(args.config, args.seed, args.resume)
    print("Resultat :", result)
