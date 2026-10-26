"""A8 — Tests de peplib.datasets (PeptideDataset et load_splits).

Couverture :
  - Chargement d'un CSV minimal, longueur, forme des tenseurs.
  - input_ids commence par BOS, labels termine par EOS.
  - Positions PAD masquees a -100 dans les labels.
  - Refus de pad_to < MAX_LEN.
  - Refus de colonne absente, fichier absent.
  - Aller-retour tokenizer : decode(input_ids[1:]) redonne la sequence.
  - Test sur corpus reel (skip si absent).
"""
import tempfile
from pathlib import Path

import pandas as pd
import pytest
import torch

from peplib.datasets import LABEL_IGNORE_INDEX, PeptideDataset
from peplib.tokenizer import BOS_ID, EOS_ID, MAX_LEN, PAD_ID, decode, encode

ROOT = Path(__file__).resolve().parents[1]

# Sequences de test (longueurs variees, toutes valides)
_TEST_SEQS = [
    "ACDEFGHIKL",        # longueur 10 (minimum)
    "ACDEFGHIKLMNPQRST", # longueur 18
    "LLLLKKKRRR",        # longueur 10
    "A" * 50,            # longueur 50 (maximum)
]


# ---------- fixture CSV temporaire ----------

@pytest.fixture
def tmp_csv(tmp_path):
    """Cree un CSV temporaire avec les sequences de test."""
    path = tmp_path / "test_seqs.csv"
    pd.DataFrame({"sequence": _TEST_SEQS}).to_csv(path, index=False)
    return path


# ---------- construction ----------

def test_dataset_len(tmp_csv):
    ds = PeptideDataset(tmp_csv)
    assert len(ds) == len(_TEST_SEQS)


def test_dataset_returns_dict(tmp_csv):
    ds = PeptideDataset(tmp_csv)
    item = ds[0]
    assert isinstance(item, dict)
    assert "input_ids" in item and "labels" in item


def test_tensors_are_longtensors(tmp_csv):
    ds = PeptideDataset(tmp_csv)
    item = ds[0]
    assert item["input_ids"].dtype == torch.long
    assert item["labels"].dtype == torch.long


def test_input_ids_and_labels_same_length(tmp_csv):
    ds = PeptideDataset(tmp_csv)
    for i in range(len(ds)):
        item = ds[i]
        assert item["input_ids"].shape == item["labels"].shape


def test_input_ids_length_equals_pad_to_minus_one(tmp_csv):
    # pad_to tokens au total ; input = tokens[:-1], labels = tokens[1:]
    ds = PeptideDataset(tmp_csv, pad_to=MAX_LEN)
    item = ds[0]
    assert item["input_ids"].shape[0] == MAX_LEN - 1
    assert item["labels"].shape[0] == MAX_LEN - 1


# ---------- contenu des tokens ----------

def test_input_ids_starts_with_bos(tmp_csv):
    ds = PeptideDataset(tmp_csv)
    for i in range(len(ds)):
        assert ds[i]["input_ids"][0].item() == BOS_ID


def test_labels_contain_eos(tmp_csv):
    """L'EOS doit apparaitre dans les labels (apres les residus)."""
    ds = PeptideDataset(tmp_csv)
    for i in range(len(ds)):
        labels = ds[i]["labels"].tolist()
        non_pad_labels = [l for l in labels if l != LABEL_IGNORE_INDEX]
        assert non_pad_labels[-1] == EOS_ID


def test_pad_positions_are_ignored_in_labels(tmp_csv):
    """Toutes les positions PAD dans les labels doivent valoir LABEL_IGNORE_INDEX."""
    ds = PeptideDataset(tmp_csv)
    for i in range(len(ds)):
        item = ds[i]
        # Les positions PAD dans input_ids[1:] correspondent aux labels masques
        labels = item["labels"].tolist()
        input_ids_shifted = item["input_ids"][1:].tolist()  # labels = tokens[1:]
        for lbl, inp in zip(labels[1:], input_ids_shifted):
            if inp == PAD_ID:
                assert lbl == LABEL_IGNORE_INDEX


# ---------- aller-retour tokenizer ----------

def test_round_trip_via_input_ids(tmp_csv):
    """decode reconstruit la sequence originale depuis les tenseurs.

    Pour les sequences de longueur < 50, EOS est present dans input_ids.
    Pour la sequence de 50 AA, pad_to=MAX_LEN=52 => tokens a 52 positions ;
    input_ids = tokens[:-1] = 51 positions = [BOS, 50 AA, EOS] -> EOS y est.
    labels = tokens[1:] = 51 positions = [50 AA, EOS] -> EOS en derniere position.
    On reconstruit [BOS ... EOS] a partir des labels si EOS absent d'input_ids.
    """
    ds = PeptideDataset(tmp_csv, pad_to=MAX_LEN)
    for i, seq in enumerate(_TEST_SEQS):
        item = ds[i]
        input_ids = item["input_ids"].tolist()
        labels = item["labels"].tolist()

        if EOS_ID in input_ids:
            eos_pos = input_ids.index(EOS_ID)
            full = input_ids[:eos_pos + 1]
        else:
            # Sequence exactement de longueur MAX_LEN-2 (50 AA) :
            # input_ids = [BOS, aa_1..aa_50, EOS] sur 52 positions tronquees a 51
            # => EOS est le dernier element de labels (non masque)
            non_ignore = [l for l in labels if l != LABEL_IGNORE_INDEX]
            assert non_ignore and non_ignore[-1] == EOS_ID, (
                f"EOS attendu en fin de labels pour la sequence {i}"
            )
            full = [BOS_ID] + non_ignore[:-1] + [EOS_ID]

        assert decode(full) == seq, f"Aller-retour echoue pour la sequence {i}"


# ---------- gestion des erreurs ----------

def test_rejects_missing_file():
    with pytest.raises(FileNotFoundError):
        PeptideDataset("data/processed/nonexistent_12345.csv")


def test_rejects_missing_column(tmp_path):
    bad_csv = tmp_path / "bad.csv"
    pd.DataFrame({"seq": ["ACDE"]}).to_csv(bad_csv, index=False)
    with pytest.raises(ValueError, match="sequence"):
        PeptideDataset(bad_csv)


def test_rejects_pad_to_too_small(tmp_csv):
    with pytest.raises(ValueError):
        PeptideDataset(tmp_csv, pad_to=MAX_LEN - 1)


def test_custom_seq_col(tmp_path):
    csv_path = tmp_path / "custom.csv"
    pd.DataFrame({"seq": _TEST_SEQS}).to_csv(csv_path, index=False)
    ds = PeptideDataset(csv_path, seq_col="seq")
    assert len(ds) == len(_TEST_SEQS)


# ---------- pad_to personnalise ----------

def test_larger_pad_to(tmp_csv):
    ds = PeptideDataset(tmp_csv, pad_to=60)
    item = ds[0]
    assert item["input_ids"].shape[0] == 59    # 60 - 1
    assert item["labels"].shape[0] == 59


# ---------- make_loader : generateur seede (A8) ----------

def test_make_loader_same_seed_same_batches_bitwise(tmp_csv):
    """Deux parcours d'un loader de meme seed : batchs identiques bit a bit.

    A8 : « meme seed, memes batchs ». On compare chaque batch (input_ids
    et labels) entre deux loaders construits avec le meme seed et de
    plusieurs seeds differents pour verifier que le melange change.
    """
    from peplib.datasets import make_loader

    dataset = PeptideDataset(tmp_csv, pad_to=MAX_LEN)

    # _TEST_SEQS doit contenir assez d'elements pour > 1 batch.
    batch_size = max(2, len(dataset) // 2)

    batches_a = list(make_loader(dataset, batch_size=batch_size, shuffle=True, seed=42))
    batches_b = list(make_loader(dataset, batch_size=batch_size, shuffle=True, seed=42))

    # Deux loaders de meme seed : ordre et contenu strictement identiques.
    assert len(batches_a) == len(batches_b)
    for ba, bb in zip(batches_a, batches_b):
        assert torch.equal(ba["input_ids"], bb["input_ids"])
        assert torch.equal(ba["labels"], bb["labels"])

    # Un seed different doit melanger autrement (continuité de la convention).
    batches_c = list(make_loader(dataset, batch_size=batch_size, shuffle=True, seed=43))
    at_least_one_differs = any(
        not torch.equal(ba["input_ids"], bc["input_ids"])
        for ba, bc in zip(batches_a, batches_c)
    )
    assert at_least_one_differs


# ---------- corpus reel ----------

def test_dataset_on_real_corpus():
    for name in ("corpus_v0.2.csv", "corpus_filtered_all.csv"):
        path = ROOT / "data" / "processed" / name
        if path.exists():
            ds = PeptideDataset(path)
            assert len(ds) >= 1000
            item = ds[0]
            assert item["input_ids"][0].item() == BOS_ID
            return
    pytest.skip("Aucun corpus disponible pour le test reel")
