"""T40 - Tokenisation V = 23 : ids stables, aller-retour 100 %, refus des cas invalides."""
import random
from pathlib import Path

import pandas as pd
import pytest

from peplib.constants import AA_LIST
from peplib.tokenizer import (
    AA_TO_ID, BOS_ID, EOS_ID, FORBIDDEN_OUTPUT_IDS, ID_TO_TOKEN, PAD_ID, VOCAB,
    decode, encode, pad,
)

ROOT = Path(__file__).resolve().parents[1]


# ---------- vocabulaire ----------
def test_vocab_size_is_23():
    assert len(VOCAB) == 23 and len(ID_TO_TOKEN) == 23


def test_special_ids():
    assert (PAD_ID, BOS_ID, EOS_ID) == (0, 1, 2)
    assert VOCAB["<PAD>"] == 0 and VOCAB["<BOS>"] == 1 and VOCAB["<EOS>"] == 2


def test_amino_acid_ids_are_fixed_alphabetical():
    # ids ecrits en dur : si l'ordre derive un jour, ce test devient rouge
    expected = dict(zip("ACDEFGHIKLMNPQRSTVWY", range(3, 23)))
    assert AA_TO_ID == expected
    assert AA_TO_ID["A"] == 3 and AA_TO_ID["G"] == 8 and AA_TO_ID["Y"] == 22
    assert "".join(AA_LIST) == "ACDEFGHIKLMNPQRSTVWY"


def test_ids_are_contiguous_and_unique():
    assert sorted(VOCAB.values()) == list(range(23))


def test_forbidden_output_ids():
    assert set(FORBIDDEN_OUTPUT_IDS) == {PAD_ID, BOS_ID}


# ---------- encode ----------
def test_encode_by_hand():
    assert encode("AG") == [1, 3, 8, 2]
    assert encode("ACDEFGHIKLMNPQRSTVWY") == [1] + list(range(3, 23)) + [2]


def test_encode_adds_bos_and_eos():
    ids = encode("KKLLK")
    assert ids[0] == BOS_ID and ids[-1] == EOS_ID and len(ids) == 7
    assert PAD_ID not in ids and ids[1:-1].count(BOS_ID) == 0


def test_encode_empty_sequence():
    assert encode("") == [1, 2]
    assert decode([1, 2]) == ""


@pytest.mark.parametrize("seq", ["AAXA", "AAUA", "AABA", "AAZA", "aaaa", "AA-A", "AA A", "AA*A"])
def test_encode_rejects_non_canonical(seq):
    with pytest.raises(ValueError):
        encode(seq)


def test_encode_accepts_exactly_50_residues():
    ids = encode("A" * 50)
    assert len(ids) == 52  # BOS + 50 residus + EOS (max_len)
    assert decode(ids) == "A" * 50
    assert pad(ids, 52) == ids


def test_encode_rejects_51_residues():
    # convention max_len = 52 (50 residus) : aucune troncature silencieuse
    with pytest.raises(ValueError):
        encode("A" * 51)
    with pytest.raises(ValueError):
        encode("A" * 60)


# ---------- decode ----------
def test_decode_by_hand():
    assert decode([1, 3, 8, 2]) == "AG"


def test_decode_accepts_trailing_pad():
    assert decode(pad(encode("AG"), 10)) == "AG"


@pytest.mark.parametrize("ids", [
    [],                    # vide
    [3, 8, 2],             # pas de BOS
    [1, 3, 8],             # pas de EOS
    [1, 3, 23, 2],         # id inconnu (hors vocabulaire)
    [1, 3, -1, 2],         # id negatif
    [1, 3, 0, 8, 2],       # PAD au milieu
    [1, 3, 1, 8, 2],       # BOS au milieu
    [1, 3, 2, 8],          # contenu apres EOS
    [1, 3, 2, 2],          # deuxieme EOS
])
def test_decode_rejects_invalid(ids):
    with pytest.raises(ValueError):
        decode(ids)


# ---------- pad ----------
def test_pad_right_and_no_truncation():
    assert pad([1, 3, 2], 5) == [1, 3, 2, 0, 0]
    assert pad([1, 3, 2], 3) == [1, 3, 2]
    with pytest.raises(ValueError):
        pad([1, 3, 2], 2)


# ---------- aller-retour 100 % ----------
def _random_sequences(n=1000):
    rnd = random.Random(42)
    return ["".join(rnd.choice(AA_LIST) for _ in range(rnd.randint(10, 50))) for _ in range(n)]


def test_round_trip_1000_random_sequences():
    seqs = _random_sequences()
    assert len(seqs) == 1000
    assert all(decode(encode(s)) == s for s in seqs)


def test_round_trip_with_padding_1000_sequences():
    seqs = _random_sequences()
    assert all(decode(pad(encode(s), 52)) == s for s in seqs)


def test_round_trip_on_real_corpus():
    for name in ("corpus_v0.2.csv", "corpus_filtered_all.csv"):
        path = ROOT / "data" / "processed" / name
        if path.exists():
            seqs = pd.read_csv(path)["sequence"].dropna().astype(str).tolist()
            break
    else:
        pytest.skip("aucun corpus dans data/processed/")
    assert len(seqs) >= 1000
    assert all(decode(encode(s)) == s for s in seqs)
