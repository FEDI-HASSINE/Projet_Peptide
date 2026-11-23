"""Tests de l'initialisation du Transformer (P0-3a : weight tying).

Le bug corrige : `_init_weights` reinitialisait `head.weight`, qui partage
son tenseur avec `token_emb.weight`. L'embedding etait donc ecrase par une
init Xavier (std ~0.115 au lieu de 0.02) et la ligne PAD n'etait plus nulle.
"""
import pytest
import torch

from peplib.models.transformer import CausalTransformer
from peplib.seed import set_seed
from peplib.tokenizer import PAD_ID, VOCAB


@pytest.fixture
def model():
    set_seed(0)
    return CausalTransformer(
        vocab_size=len(VOCAB),
        embed_dim=32,
        n_heads=4,
        n_layers=2,
        ffn_dim=64,
        dropout=0.0,
    )


def test_weight_tying_shares_storage(model):
    """head.weight et token_emb.weight sont le meme tenseur."""
    assert model.head.weight.data_ptr() == model.token_emb.weight.data_ptr()


def test_pad_embedding_is_zero(model):
    """La ligne PAD reste nulle apres initialisation."""
    assert torch.count_nonzero(model.token_emb.weight[PAD_ID]).item() == 0


def test_embedding_std_is_0_02(model):
    """Les lignes non-PAD gardent l'init N(0, 0.02), pas Xavier."""
    rows = torch.cat([model.token_emb.weight[:PAD_ID], model.token_emb.weight[PAD_ID + 1:]])
    std = rows.std().item()
    assert abs(std - 0.02) < 0.005, f"std={std:.4f} : embedding ecrase par une init Xavier"


def test_head_bias_is_zero(model):
    assert torch.count_nonzero(model.head.bias).item() == 0


def test_reinitialisation_keeps_pad_zero(model):
    """Le bug : un second _init_weights ecrasait l'embedding via head.weight."""
    model._init_weights()
    assert torch.count_nonzero(model.token_emb.weight[PAD_ID]).item() == 0
    rows = torch.cat([model.token_emb.weight[:PAD_ID], model.token_emb.weight[PAD_ID + 1:]])
    assert abs(rows.std().item() - 0.02) < 0.005


def test_initialisation_is_reproducible():
    set_seed(123)
    first = CausalTransformer(vocab_size=len(VOCAB), embed_dim=16, n_heads=4, n_layers=1, ffn_dim=32)
    set_seed(123)
    second = CausalTransformer(vocab_size=len(VOCAB), embed_dim=16, n_heads=4, n_layers=1, ffn_dim=32)
    assert torch.equal(first.token_emb.weight, second.token_emb.weight)
