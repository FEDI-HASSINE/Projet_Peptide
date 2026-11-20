"""
tests/test_generate.py -- Tests unitaires de peplib.generate (correction PR #25).

Verifie :
  1. generate_sequences retourne des sequences valides (alphabet, 10-50 aa)
  2. PAD et BOS ne sont jamais dans les sequences generees
  3. Reproductibilite avec seed fixe
  4. Fonctionne pour LSTM et Transformer
  5. Respect de max_len
"""
import pytest
import torch

from peplib.generate import generate_sequences
from peplib.models.lstm import CausalLSTM
from peplib.models.transformer import CausalTransformer
from peplib.seed import set_seed
from peplib.tokenizer import VOCAB, PAD_ID, BOS_ID, EOS_ID
from peplib.filters import classify_sequence

VOCAB_SIZE = len(VOCAB)  # 23


@pytest.fixture
def small_lstm():
    set_seed(0)
    return CausalLSTM(
        vocab_size=VOCAB_SIZE, embed_dim=16, hidden_dim=32,
        n_layers=1, dropout=0.0,
    ).eval()


@pytest.fixture
def small_transformer():
    set_seed(0)
    return CausalTransformer(
        vocab_size=VOCAB_SIZE, embed_dim=16, n_heads=4,
        n_layers=2, ffn_dim=32, dropout=0.0, max_len=52,
    ).eval()


# -----------------------------------------------------------------------
# Test 1 : sequences valides (alphabet + longueur)
# -----------------------------------------------------------------------

def test_lstm_generates_valid_sequences(small_lstm):
    """Les sequences generees par le LSTM doivent passer classify_sequence."""
    seqs = generate_sequences(small_lstm, n=10, temperature=1.0, seed=42)
    assert len(seqs) > 0, "Aucune sequence valide generee par le LSTM"
    for s in seqs:
        reason = classify_sequence(s)
        assert reason is None, f"Sequence invalide generee : {s!r} -> {reason}"


def test_transformer_generates_valid_sequences(small_transformer):
    """Les sequences generees par le Transformer doivent passer classify_sequence."""
    seqs = generate_sequences(small_transformer, n=10, temperature=1.0, seed=42)
    assert len(seqs) > 0, "Aucune sequence valide generee par le Transformer"
    for s in seqs:
        reason = classify_sequence(s)
        assert reason is None, f"Sequence invalide generee : {s!r} -> {reason}"


# -----------------------------------------------------------------------
# Test 2 : PAD et BOS absents
# -----------------------------------------------------------------------

def test_no_pad_or_bos_in_output(small_lstm):
    """Les sequences retournees ne doivent contenir ni PAD ni BOS."""
    seqs = generate_sequences(small_lstm, n=20, temperature=1.0, seed=1)
    for s in seqs:
        assert "<PAD>" not in s, f"PAD dans la sequence : {s!r}"
        assert "<BOS>" not in s, f"BOS dans la sequence : {s!r}"
        assert "<EOS>" not in s, f"EOS dans la sequence : {s!r}"
        # Verification que les caracteres sont bien des AA canoniques
        assert all(c in "ACDEFGHIKLMNPQRSTVWY" for c in s), (
            f"Caractere non canonique dans : {s!r}"
        )


# -----------------------------------------------------------------------
# Test 3 : Reproductibilite
# -----------------------------------------------------------------------

def test_reproducibility_lstm(small_lstm):
    """Meme seed -> memes sequences."""
    set_seed(42)
    seqs1 = generate_sequences(small_lstm, n=5, temperature=1.0, seed=42)
    set_seed(42)
    seqs2 = generate_sequences(small_lstm, n=5, temperature=1.0, seed=42)
    assert seqs1 == seqs2, "Les sequences ne sont pas reproductibles avec la meme seed."


# -----------------------------------------------------------------------
# Test 4 : max_len respectee
# -----------------------------------------------------------------------

def test_max_len_respected(small_lstm):
    """Les sequences generees ne depassent pas max_len - 2 residus (sans BOS/EOS)."""
    max_len = 20  # 20 tokens au total (BOS + max 18 AA + EOS)
    seqs = generate_sequences(small_lstm, n=20, temperature=1.0, seed=5, max_len=max_len)
    for s in seqs:
        # max_len - 2 car BOS et EOS comptent chacun pour 1
        assert len(s) <= max_len - 2, (
            f"Sequence trop longue ({len(s)} > {max_len-2}) : {s!r}"
        )


# -----------------------------------------------------------------------
# Test 5 : n=0 retourne liste vide
# -----------------------------------------------------------------------

def test_n_zero_returns_empty(small_lstm):
    seqs = generate_sequences(small_lstm, n=0, seed=0)
    assert seqs == []
