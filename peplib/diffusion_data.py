"""Lots de données S4. Un temps entier est partagé par toutes les lignes du lot."""
import csv
from pathlib import Path

import torch
from torch.utils.data import Dataset

from peplib.tokenizer_diffusion import BOS_ID, MASK_ID, MAX_RESIDUES, PAD_ID, encode


class DiffusionDataset(Dataset):
    """Choix explicite train/val. Chaque ligne contient 50 emplacements + EOS."""

    def __init__(self, split, directory="data/splits"):
        if split not in ("train", "val"):
            raise ValueError("Split autorisé attendu (train ou val)")
        path = Path(directory) / f"{split}_v02.csv"
        rows = []
        with path.open(encoding="utf-8", newline="") as stream:
            reader = csv.DictReader(stream)
            if "sequence" not in (reader.fieldnames or []):
                raise ValueError("Colonne sequence absente")
            for row in reader:
                ids = encode(row["sequence"])[1:]
                rows.append(ids + [PAD_ID] * (MAX_RESIDUES + 1 - len(ids)))
        if not rows:
            raise ValueError("Split vide")
        self.ids = torch.tensor(rows, dtype=torch.long)

    def __len__(self):
        return len(self.ids)

    def __getitem__(self, index):
        return self.ids[index]


def mask_batch(x0, t=None, generator=None, *, T=50):
    """Retourne (x_t, masque booléen, temps entier). EOS est masquable."""
    if x0.dtype != torch.long or x0.ndim != 2 or not x0.numel():
        raise ValueError("x0 doit être un tenseur int64 non vide de forme (B, L)")
    if ((x0 < 0) | (x0 >= MASK_ID)).any():
        raise ValueError("x0 contient des identifiants invalides ou déjà masqués")
    if type(T) is not int or T < 1:
        raise ValueError("T doit être un entier positif")
    if t is None:
        t = torch.randint(1, T + 1, (), device=x0.device, generator=generator).item()
    if type(t) is not int or not 0 <= t <= T:
        raise ValueError("t doit être un entier entre 0 et T")
    eligible = (x0 != PAD_ID) & (x0 != BOS_ID)
    mask = torch.zeros_like(eligible)
    draws = torch.rand(int(eligible.sum()), generator=generator, device=x0.device)
    mask[eligible] = draws < t / T
    return x0.masked_fill(mask, MASK_ID), mask, t
