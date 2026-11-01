"""
tests/test_lstm.py — Tests unitaires du LSTM causal (T52).

Vérifie :
  1. Forme de la sortie (B, L, vocab_size)
  2. Causalité : modifier le token t ne change pas les logits aux positions < t
  3. Gradient passe sans erreur (backward)
  4. Reproductibilité : même seed → même sortie
  5. generate_step : forme et cohérence avec forward
"""
import pytest
import torch

from peplib.models.lstm import CausalLSTM
from peplib.seed import set_seed
from peplib.tokenizer import VOCAB, PAD_ID, BOS_ID, EOS_ID


VOCAB_SIZE = len(VOCAB)  # 23
B, L = 4, 10             # batch_size, seq_len pour les tests


@pytest.fixture
def model() -> CausalLSTM:
    set_seed(0)
    return CausalLSTM(
        vocab_size=VOCAB_SIZE,
        embed_dim=32,
        hidden_dim=64,
        n_layers=2,
        dropout=0.0,  # dropout=0 pour tests déterministes
    )


@pytest.fixture
def synthetic_input() -> torch.Tensor:
    """Batch synthétique : BOS + résidus AA valides + EOS + PAD."""
    set_seed(1)
    ids = torch.randint(3, VOCAB_SIZE, (B, L))  # AA ids (3..22)
    ids[:, 0] = BOS_ID    # BOS en position 0
    ids[:, -1] = EOS_ID   # EOS en dernière position
    return ids


# -----------------------------------------------------------------------
# Test 1 : Forme de la sortie
# -----------------------------------------------------------------------

def test_output_shape(model, synthetic_input):
    """La sortie doit être (B, L, vocab_size)."""
    model.eval()
    with torch.no_grad():
        logits, hidden = model(synthetic_input)
    assert logits.shape == (B, L, VOCAB_SIZE), (
        f"Forme attendue ({B}, {L}, {VOCAB_SIZE}), obtenu {logits.shape}"
    )


def test_hidden_state_shape(model, synthetic_input):
    """L'état caché retourné doit avoir la bonne forme."""
    model.eval()
    with torch.no_grad():
        _, (h, c) = model(synthetic_input)
    assert h.shape == (model.n_layers, B, model.hidden_dim)
    assert c.shape == (model.n_layers, B, model.hidden_dim)


# -----------------------------------------------------------------------
# Test 2 : Causalité
# -----------------------------------------------------------------------

def test_causality(model):
    """Modifier le token t ne doit pas changer les logits aux positions 0..t-1.

    Stratégie : on compare les logits d'une séquence originale avec les logits
    d'une copie où l'on a modifié le token à la position T. Les positions < T
    doivent être identiques ; les positions >= T peuvent changer.
    """
    model.eval()
    T = 5  # position modifiée

    set_seed(2)
    ids_orig = torch.randint(3, VOCAB_SIZE, (1, L))
    ids_orig[0, 0] = BOS_ID

    ids_mod = ids_orig.clone()
    # Modifier le token à la position T (choisir un id différent)
    new_tok = (ids_orig[0, T].item() % (VOCAB_SIZE - 3)) + 3
    ids_mod[0, T] = new_tok

    with torch.no_grad():
        logits_orig, _ = model(ids_orig)
        logits_mod,  _ = model(ids_mod)

    # Les positions 0..T-1 doivent être IDENTIQUES
    assert torch.allclose(logits_orig[0, :T], logits_mod[0, :T], atol=1e-5), (
        "Violation de causalité : les logits avant la position modifiée ont changé."
    )


# -----------------------------------------------------------------------
# Test 3 : Gradient
# -----------------------------------------------------------------------

def test_gradient_flows(model, synthetic_input):
    """Le backward doit passer sans erreur et produire des gradients non nuls."""
    model.train()
    logits, _ = model(synthetic_input)  # (B, L, V)

    # Cible : décalage d'un token vers la droite (next-token prediction)
    target = synthetic_input[:, 1:]  # (B, L-1)
    loss = torch.nn.functional.cross_entropy(
        logits[:, :-1].reshape(-1, VOCAB_SIZE),
        target.reshape(-1),
        ignore_index=PAD_ID,
    )
    loss.backward()

    # Vérifier que les gradients de l'embedding ne sont pas tous nuls
    grad = model.embedding.weight.grad
    assert grad is not None, "Pas de gradient sur l'embedding."
    assert grad.abs().sum().item() > 0, "Gradient nul sur l'embedding."


# -----------------------------------------------------------------------
# Test 4 : Reproductibilité
# -----------------------------------------------------------------------

def test_reproducibility():
    """Même seed → même sortie (forward deterministe sur CPU)."""
    ids = torch.randint(3, VOCAB_SIZE, (2, 8))
    ids[:, 0] = BOS_ID

    set_seed(42)
    m1 = CausalLSTM(vocab_size=VOCAB_SIZE, embed_dim=16, hidden_dim=32, n_layers=1, dropout=0.0)
    m1.eval()
    with torch.no_grad():
        out1, _ = m1(ids)

    set_seed(42)
    m2 = CausalLSTM(vocab_size=VOCAB_SIZE, embed_dim=16, hidden_dim=32, n_layers=1, dropout=0.0)
    m2.eval()
    with torch.no_grad():
        out2, _ = m2(ids)

    assert torch.allclose(out1, out2, atol=1e-7), (
        "Les deux modèles avec la même seed donnent des sorties différentes."
    )


# -----------------------------------------------------------------------
# Test 5 : generate_step
# -----------------------------------------------------------------------

def test_generate_step_shape(model):
    """generate_step doit retourner (B, 1, vocab_size) et un état caché valide."""
    model.eval()
    token = torch.tensor([[BOS_ID]] * B)  # (B, 1)
    with torch.no_grad():
        logits, hidden = model.generate_step(token, hidden=None)
    assert logits.shape == (B, 1, VOCAB_SIZE)
    assert hidden is not None


def test_generate_step_consistent_with_forward(model):
    """generate_step sur chaque token doit donner les mêmes logits que forward complet."""
    model.eval()
    ids = torch.randint(3, VOCAB_SIZE, (1, L))
    ids[0, 0] = BOS_ID

    with torch.no_grad():
        logits_full, _ = model(ids)  # (1, L, V)

        # Parcours token par token
        hidden = None
        for t in range(L):
            token = ids[:, t:t+1]  # (1, 1)
            logits_step, hidden = model.generate_step(token, hidden)
            assert torch.allclose(logits_full[:, t:t+1], logits_step, atol=1e-5), (
                f"Incohérence entre forward et generate_step à la position {t}."
            )


# -----------------------------------------------------------------------
# Test 6 : count_parameters
# -----------------------------------------------------------------------

def test_count_parameters(model):
    n = model.count_parameters()
    assert n > 0, "Le modèle doit avoir des paramètres entraînables."
    # Vérification grossière de l'ordre de grandeur
    # embed: 23*32 + LSTM + head
    assert 1_000 < n < 10_000_000, f"Nombre de paramètres suspect : {n}"
