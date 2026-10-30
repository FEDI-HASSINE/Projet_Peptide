"""
Wrapper PyTorch pour les splits de peptides (A8).

Charge un fichier CSV (colonne 'sequence'), encode chaque sequence via
peplib.tokenizer (BOS + residus + EOS + PAD), et retourne des tenseurs
adaptes a l'entrainement causal (next-token prediction).

Structure d'un batch :
    input_ids : LongTensor (L,)   BOS aa_1 ... aa_n EOS [PAD ...]
    labels    : LongTensor (L,)   aa_1 ... aa_n EOS [PAD ...] (-100 sur PAD)

La valeur -100 sur les positions PAD est la convention PyTorch
(torch.nn.CrossEntropyLoss ignore_index=-100 par defaut).

Compatibilite : torch >= 2.0, pandas >= 1.5.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import pandas as pd
import torch
from torch.utils.data import Dataset

from peplib.tokenizer import MAX_LEN, PAD_ID, encode, pad

LABEL_IGNORE_INDEX: int = -100  # convention CrossEntropyLoss


class PeptideDataset(Dataset):
    """Dataset PyTorch pour les sequences de peptides.

    Parameters
    ----------
    csv_path : str | Path
        Chemin vers un CSV contenant au moins une colonne 'sequence'.
    pad_to : int
        Longueur de remplissage (tokens). Doit etre >= MAX_LEN (52).
        Par defaut MAX_LEN.
    seq_col : str
        Nom de la colonne de sequences dans le CSV.

    Raises
    ------
    FileNotFoundError
        Si le fichier CSV est introuvable.
    ValueError
        Si la colonne de sequences est absente, ou si pad_to < MAX_LEN.
    """

    def __init__(
        self,
        csv_path: str | Path,
        pad_to: int = MAX_LEN,
        seq_col: str = "sequence",
    ) -> None:
        csv_path = Path(csv_path)
        if not csv_path.exists():
            raise FileNotFoundError(f"Fichier introuvable : {csv_path}")
        if pad_to < MAX_LEN:
            raise ValueError(
                f"pad_to ({pad_to}) < MAX_LEN ({MAX_LEN}) : "
                "le remplissage ne peut pas tronquer les sequences encodees."
            )

        df = pd.read_csv(csv_path, dtype=str)
        if seq_col not in df.columns:
            raise ValueError(
                f"Colonne '{seq_col}' absente du CSV. "
                f"Colonnes disponibles : {list(df.columns)}"
            )

        sequences = df[seq_col].dropna().astype(str).tolist()

        # Encode et remplit toutes les sequences a la meme longueur.
        self._pad_to = pad_to
        self._ids: list[list[int]] = []

        skipped = 0
        for seq in sequences:
            try:
                ids = encode(seq)
                ids = pad(ids, pad_to)
                self._ids.append(ids)
            except ValueError:
                # Sequence invalide (trop longue, alphabet) -> ignoree silencieusement.
                skipped += 1

        if skipped:
            import warnings
            warnings.warn(
                f"{skipped} sequence(s) invalide(s) ignoree(s) dans {csv_path.name}.",
                stacklevel=2,
            )

        if not self._ids:
            raise ValueError(f"Aucune sequence valide dans {csv_path}")

    # ------------------------------------------------------------------
    # Interface Dataset
    # ------------------------------------------------------------------

    def __len__(self) -> int:
        return len(self._ids)

    def __getitem__(self, idx: int) -> dict[str, torch.Tensor]:
        """Retourne un dictionnaire avec input_ids et labels.

        input_ids : tokens[:-1]  (de BOS au dernier AA inclus)
        labels    : tokens[1:]   avec -100 sur les positions PAD

        Pour un modele causal, le token i predit le token i+1.
        Les positions PAD en labels sont masquees (-100) : la loss
        n'y est pas calculee.
        """
        ids = self._ids[idx]
        input_ids = torch.tensor(ids[:-1], dtype=torch.long)
        labels_raw = torch.tensor(ids[1:], dtype=torch.long)
        # Masquer les PAD dans les labels
        labels = torch.where(labels_raw == PAD_ID, LABEL_IGNORE_INDEX, labels_raw)
        return {"input_ids": input_ids, "labels": labels}

    # ------------------------------------------------------------------
    # Utilitaires
    # ------------------------------------------------------------------

    @property
    def pad_to(self) -> int:
        """Longueur de remplissage utilisee."""
        return self._pad_to

    def __repr__(self) -> str:  # pragma: no cover
        return (
            f"PeptideDataset("
            f"n={len(self)}, pad_to={self._pad_to})"
        )


# ----------------------------------------------------------------------
# Fonction utilitaire : creer les trois DataLoaders train/val/test
# ----------------------------------------------------------------------

def make_loader(
    dataset: "PeptideDataset",
    batch_size: int = 64,
    shuffle: bool = False,
    seed: int = 42,
    num_workers: int = 0,
    drop_last: bool = False,
) -> "torch.utils.data.DataLoader":
    """Construit un DataLoader avec un generateur seede (reproductible).

    Le generateur (torch.Generator().manual_seed(seed)) gouverne le
    melange des batchs : deux loaders construits avec le meme seed
    produisent exactement les memes batchs dans le meme ordre (A8 :
    « meme seed, memes batchs »).
    """
    from torch.utils.data import DataLoader

    generator = torch.Generator()
    generator.manual_seed(seed)
    return DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=shuffle,
        num_workers=num_workers,
        drop_last=drop_last,
        generator=generator,
    )


def load_splits(
    splits_dir: str | Path = "data/processed/splits",
    pad_to: int = MAX_LEN,
    batch_size: int = 64,
    num_workers: int = 0,
    seed: int = 42,
) -> dict[str, "torch.utils.data.DataLoader"]:
    """Cree un DataLoader pour chacun des fichiers train/val/test.

    Attend les fichiers ``{train,val,test}.csv`` dans ``splits_dir``.
    Retourne un dict ``{"train": ..., "val": ..., "test": ...}``.
    Leve FileNotFoundError si un split est manquant.

    Parameters
    ----------
    splits_dir : str | Path
        Repertoire contenant train.csv, val.csv, test.csv.
    pad_to : int
        Longueur de remplissage passee a PeptideDataset.
    batch_size : int
        Taille du batch.
    num_workers : int
        Nombre de workers DataLoader (0 = monothread, recommande sur Windows).
    seed : int
        Graine du generateur de melange, passee a make_loader (A8 :
        « meme seed, memes batchs »).
    """
    from torch.utils.data import DataLoader

    splits_dir = Path(splits_dir)
    loaders: dict[str, DataLoader] = {}

    for split in ("train", "val", "test"):
        path = splits_dir / f"{split}.csv"
        if not path.exists():
            raise FileNotFoundError(
                f"Split '{split}' introuvable : {path}. "
                "Verifier que A6 (splits.py de Sofiene) a ete execute."
            )
        dataset = PeptideDataset(path, pad_to=pad_to)
        shuffle = split == "train"
        loaders[split] = make_loader(
            dataset,
            batch_size=batch_size,
            shuffle=shuffle,
            seed=seed,
            num_workers=num_workers,
        )

    return loaders

def collate_peptides(batch):
    """Padding dynamique, masque et longueurs pour un batch."""
    from torch.nn.utils.rnn import pad_sequence

    inputs = []
    targets = []
    lengths = []

    for item in batch:
        x = item["input_ids"]
        y = item["labels"]

        # Les labels PAD sont deja remplaces par -100.
        length = int((y != LABEL_IGNORE_INDEX).sum().item())
        inputs.append(x[:length])
        targets.append(y[:length])
        lengths.append(length)

    input_ids = pad_sequence(
        inputs, batch_first=True, padding_value=PAD_ID
    )
    labels = pad_sequence(
        targets, batch_first=True, padding_value=LABEL_IGNORE_INDEX
    )
    lengths = torch.tensor(lengths, dtype=torch.long)
    pad_mask = (
        torch.arange(input_ids.shape[1]).unsqueeze(0)
        < lengths.unsqueeze(1)
    )

    return {
        "input_ids": input_ids,
        "labels": labels,
        "lengths": lengths,
        "pad_mask": pad_mask,
    }


def make_loader_s3(
    dataset,
    batch_size=128,
    shuffle=False,
    seed=42,
    num_workers=0,
):
    """Chargeur S3 reproductible avec padding dynamique."""
    from torch.utils.data import DataLoader

    generator = torch.Generator().manual_seed(seed)

    return DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=shuffle,
        num_workers=num_workers,
        generator=generator,
        collate_fn=collate_peptides,
    )
