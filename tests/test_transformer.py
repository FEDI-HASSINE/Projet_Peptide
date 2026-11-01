"""
tests/test_transformer.py — Tests unitaires du Transformer causal (T55).

Vérifie :
  1. Forme de la sortie (B, L, vocab_size)
  2. Masque causal : logit[i] ne dépend pas des tokens > i
  3. Gradient passe sans erreur (backward)
  4. Reproductibilité : même seed → même sortie
  5. Comportement avec séquences paddées (mask de padding actif)
"""
import pytest
import torch

from peplib.models.transformer import CausalTransformer
from peplib.seed import set_seed
from peplib.tokenizer import VOCAB, PAD_ID, BOS_ID, EOS_ID


VOCAB_SIZE = len(VOCAB)  # 23
B, L = 4, 10


@pytest.fixture
def model() -> CausalTransformer:
    set_seed(0)
    return CausalTransformer(
        vocab_size=VOCAB_SIZE,
        embed_dim=32,
        n_heads=4,
        n_layers=2,
        ffn_dim=64,
        dropout=0.0,   # dropout=0 pour tests déterministes
        max_len=52,
    )


@pytest.fixture
def synthetic_input() -> torch.Tensor:
    set_seed(1)
    ids = torch.randint(3, VOCAB_SIZE, (B, L))
    ids[:, 0] = BOS_ID
    ids[:, -1] = EOS_ID
    return ids


# -----------------------------------------------------------------------
# Test 1 : Forme de la sortie
# -----------------------------------------------------------------------

def test_output_shape(model, synthetic_input):
    """La sortie doit être (B, L, vocab_size)."""
    model.eval()
    with torch.no_grad():
        logits = model(synthetic_input)
    assert logits.shape == (B, L, VOCAB_SIZE), (
        f"Forme attendue ({B}, {L}, {VOCAB_SIZE}), obtenu {logits.shape}"
    )


def test_output_dtype(model, synthetic_input):
    """Les logits doivent être des float32."""
    model.eval()
    with torch.no_grad():
        logits = model(synthetic_input)
    assert logits.dtype == torch.float32


# -----------------------------------------------------------------------
# Test 2 : Masque causal
# -----------------------------------------------------------------------

def test_causal_mask(model):
    """Le logit à la position i ne doit pas dépendre des tokens aux positions > i.

    Stratégie : on compare les logits d'une séquence originale avec une copie
    où l'on a modifié le token à la position T+1. Les logits aux positions 0..T
    doivent être identiques.
    """
    model.eval()
    T = 4  # position de référence (les positions <= T ne doivent pas changer)

    set_seed(3)
    ids_orig = torch.randint(3, VOCAB_SIZE, (1, L))
    ids_orig[0, 0] = BOS_ID

    ids_mod = ids_orig.clone()
    # Modifier le token à la position T+1 (future)
    future_pos = T + 1
    new_tok = (ids_orig[0, future_pos].item() % (VOCAB_SIZE - 3)) + 3
    ids_mod[0, future_pos] = new_tok

    with torch.no_grad():
        logits_orig = model(ids_orig)  # (1, L, V)
        logits_mod  = model(ids_mod)   # (1, L, V)

    # Les logits aux positions 0..T doivent être IDENTIQUES
    assert torch.allclose(logits_orig[0, :T+1], logits_mod[0, :T+1], atol=1e-4), (
        f"Violation du masque causal : modifier la position {future_pos} "
        f"a changé les logits aux positions <= {T}."
    )


def test_causal_mask_generated_correctly():
    """Le masque causal doit avoir -inf au-dessus de la diagonale et 0 ailleurs."""
    L_test = 5
    mask = CausalTransformer._causal_mask(L_test, torch.device("cpu"))
    assert mask.shape == (L_test, L_test)
    # Diagonale et au-dessous → 0
    for i in range(L_test):
        for j in range(i + 1):
            assert mask[i, j] == 0.0 or not mask[i, j].isinf(), (
                f"Position ({i},{j}) devrait être 0, obtenu {mask[i,j]}"
            )
    # Au-dessus de la diagonale → -inf
    for i in range(L_test):
        for j in range(i + 1, L_test):
            assert mask[i, j] == float("-inf"), (
                f"Position ({i},{j}) devrait être -inf, obtenu {mask[i,j]}"
            )


# -----------------------------------------------------------------------
# Test 3 : Gradient
# -----------------------------------------------------------------------

def test_gradient_flows(model, synthetic_input):
    """Le backward doit passer sans erreur et produire des gradients non nuls."""
    model.train()
    logits = model(synthetic_input)  # (B, L, V)

    target = synthetic_input[:, 1:]  # (B, L-1)
    loss = torch.nn.functional.cross_entropy(
        logits[:, :-1].reshape(-1, VOCAB_SIZE),
        target.reshape(-1),
        ignore_index=PAD_ID,
    )
    loss.backward()

    # Gradient sur l'embedding token (exclut les lignes de poids partagés avec head)
    grad = model.token_emb.weight.grad
    assert grad is not None, "Pas de gradient sur token_emb."
    assert grad.abs().sum().item() > 0, "Gradient nul sur token_emb."


# -----------------------------------------------------------------------
# Test 4 : Reproductibilité
# -----------------------------------------------------------------------

def test_reproducibility():
    """Même seed → même forward sur CPU."""
    ids = torch.randint(3, VOCAB_SIZE, (2, 8))
    ids[:, 0] = BOS_ID

    set_seed(42)
    m1 = CausalTransformer(
        vocab_size=VOCAB_SIZE, embed_dim=16, n_heads=4,
        n_layers=2, ffn_dim=32, dropout=0.0, max_len=52,
    )
    m1.eval()
    with torch.no_grad():
        out1 = m1(ids)

    set_seed(42)
    m2 = CausalTransformer(
        vocab_size=VOCAB_SIZE, embed_dim=16, n_heads=4,
        n_layers=2, ffn_dim=32, dropout=0.0, max_len=52,
    )
    m2.eval()
    with torch.no_grad():
        out2 = m2(ids)

    assert torch.allclose(out1, out2, atol=1e-7), (
        "Deux modèles avec la même seed donnent des sorties différentes."
    )


# -----------------------------------------------------------------------
# Test 5 : Masque de padding
# -----------------------------------------------------------------------

def test_padding_mask(model):
    """Les positions PAD ne doivent pas influencer les logits des positions valides
    d'une même séquence (côté clé — key_padding_mask actif)."""
    model.eval()

    # Séquence A : 5 tokens valides + 5 PAD
    ids_a = torch.full((1, L), PAD_ID, dtype=torch.long)
    ids_a[0, :5] = torch.tensor([BOS_ID, 5, 6, 7, EOS_ID])

    # Séquence B : identique mais avec des tokens différents dans les positions PAD
    ids_b = ids_a.clone()
    ids_b[0, 5:] = 8  # valeurs différentes là où c'est PAD dans A

    with torch.no_grad():
        logits_a = model(ids_a)
        logits_b = model(ids_b)

    # Les logits aux positions 0..4 (valides) doivent être identiques
    assert torch.allclose(logits_a[0, :5], logits_b[0, :5], atol=1e-4), (
        "Les positions PAD influencent les logits des positions valides."
    )


# -----------------------------------------------------------------------
# Test 6 : count_parameters + embed_dim divisible par n_heads
# -----------------------------------------------------------------------

def test_count_parameters(model):
    n = model.count_parameters()
    assert n > 0
    assert 1_000 < n < 10_000_000, f"Paramètres suspects : {n}"


def test_embed_dim_not_divisible_raises():
    """embed_dim non divisible par n_heads doit lever ValueError."""
    with pytest.raises(ValueError, match="divisible"):
        CausalTransformer(vocab_size=23, embed_dim=10, n_heads=3)
