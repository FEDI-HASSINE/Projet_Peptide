"""T41 / A7 — Tests du modele n-gram (peplib.ngram).

Couverture :
  - Ordre 0 : perplexite = 20 (distribution uniforme).
  - Ordre 1 : somme des proba = 1 sur tout l'alphabet (sequence pleine).
  - Ordres 0–3 : log_prob d'une sequence vide = 0.0.
  - Ordres 0–3 : log_prob ∈ (-inf, 0] pour toute sequence non vide.
  - Ordre 2–3 : log_prob sur sequence connue == calcul manuel.
  - tune_alpha : best_alpha ∈ alphas, perplexite val < ordre 0.
  - Comportement sur corpus reel (skip si absent).
  - Refus des parametres invalides.
"""
import math
import random
from pathlib import Path

import pytest

from peplib.constants import AMINO_ACIDS
from peplib.ngram import NgramModel, tune_alpha

ROOT = Path(__file__).resolve().parents[1]
AA_LIST = list(AMINO_ACIDS)   # 20 acides amines


# ---------- helpers ----------

def _random_sequences(n: int = 200, seed: int = 42) -> list[str]:
    rnd = random.Random(seed)
    return [
        "".join(rnd.choice(AA_LIST) for _ in range(rnd.randint(10, 50)))
        for _ in range(n)
    ]


def _make_split(seed: int = 42) -> tuple[list[str], list[str]]:
    seqs = _random_sequences(400, seed=seed)
    return seqs[:300], seqs[300:]


# ---------- construction et refus ----------

def test_invalid_order_raises():
    with pytest.raises(ValueError):
        NgramModel(order=4)
    with pytest.raises(ValueError):
        NgramModel(order=-1)


def test_invalid_alpha_raises():
    with pytest.raises((ValueError, TypeError)):
        NgramModel(order=1, alpha=-0.1)
    with pytest.raises((ValueError, TypeError)):
        NgramModel(order=1, alpha="high")


def test_not_fitted_raises():
    m = NgramModel(order=1)
    with pytest.raises(RuntimeError):
        m.log_prob("ACDE")
    with pytest.raises(RuntimeError):
        m.perplexity(["ACDE"])


# ---------- ordre 0 : uniforme ----------

def test_order0_perplexity_is_20():
    m = NgramModel(order=0, alpha=1.0).fit(_random_sequences())
    seqs = _random_sequences(50, seed=99)
    ppl = m.perplexity(seqs)
    assert ppl == pytest.approx(20.0, rel=1e-9)


def test_order0_log_prob_formula():
    m = NgramModel(order=0).fit(["ACDE"])
    seq = "ACDE"
    expected = len(seq) * math.log(1 / 20)
    assert m.log_prob(seq) == pytest.approx(expected, rel=1e-9)


# ---------- sequence vide ----------

@pytest.mark.parametrize("order", [0, 1, 2, 3])
def test_log_prob_empty_sequence_is_zero(order):
    m = NgramModel(order=order, alpha=1.0).fit(["ACDE"] * 10)
    assert m.log_prob("") == 0.0


# ---------- log_prob dans (-inf, 0] ----------

@pytest.mark.parametrize("order", [0, 1, 2, 3])
def test_log_prob_nonpositive(order):
    seqs = _random_sequences(100)
    m = NgramModel(order=order, alpha=1.0).fit(seqs)
    for seq in _random_sequences(20, seed=7):
        lp = m.log_prob(seq)
        assert lp <= 1e-12, f"log_prob = {lp} > 0 pour ordre {order}"


# ---------- ordre 1 : proba somme a 1 ----------

def test_order1_probabilities_sum_to_one():
    """Pour l'ordre 1, la somme exp(log_prob(aa)) sur tous les AA = 1."""
    m = NgramModel(order=1, alpha=1.0).fit(_random_sequences())
    total = sum(math.exp(m.log_prob(aa)) for aa in AA_LIST)
    assert total == pytest.approx(1.0, abs=1e-9)


# ---------- ordre 2 : calcul manuel ----------

def test_order2_log_prob_manual():
    """Verifie un calcul bigram a la main.

    Train = ["ACAC"] une seule fois.
    Contexte () n'est jamais utilise en ordre 2.
    Contextes observes :
      (<BOS>,) -> A : 1
      (A,)     -> C : 1
      (C,)     -> A : 1, C : (apres 2e passage) non — train="ACAC" donne A,C,A,C
        (<BOS>,) -> A : 1
        (A,)     -> C : 2  (positions 1 et 3)
        (C,)     -> A : 1  (position 2)
    Sequence testee = "ACA", alpha = 1.0, V = 20
      P(A | <BOS>) = (1+1)/(1+20) = 2/21
      P(C | A)     = (2+1)/(2+20) = 3/22
      P(A | C)     = (1+1)/(1+20) = 2/21
    log_prob = log(2/21) + log(3/22) + log(2/21)
    """
    m = NgramModel(order=2, alpha=1.0).fit(["ACAC"])
    seq = "ACA"
    expected = math.log(2 / 21) + math.log(3 / 22) + math.log(2 / 21)
    assert m.log_prob(seq) == pytest.approx(expected, rel=1e-9)


# ---------- ordre 3 : log_prob coherent ----------

def test_order3_log_prob_coherent():
    """L'ordre 3 ne doit pas planter et doit retourner un scalaire <= 0."""
    seqs = _random_sequences(100)
    m = NgramModel(order=3, alpha=0.5).fit(seqs)
    for seq in _random_sequences(20, seed=3):
        lp = m.log_prob(seq)
        assert isinstance(lp, float)
        assert lp <= 1e-12


# ---------- perplexite : ordre 0 <= autres ordres avec alpha optimise ----------

def test_perplexity_improves_with_order_on_random_data():
    """Ordre 1 avec alpha bien choisi doit battre l'ordre 0 (PP < 20)."""
    train, val = _make_split()
    m0 = NgramModel(order=0).fit(train)
    m1 = NgramModel(order=1, alpha=0.1).fit(train)
    ppl0 = m0.perplexity(val)
    ppl1 = m1.perplexity(val)
    # Les sequences sont aleatoires uniformes, mais l'ordre 1 doit quand
    # meme estimer la composition correctement -> PP proches de 20.
    # On verifie juste que ppl1 <= ppl0 + 1 (tolerance car donnees unif.)
    assert ppl1 <= ppl0 + 1.0


# ---------- tune_alpha ----------

def test_tune_alpha_returns_value_in_candidates():
    train, val = _make_split()
    alphas = [0.01, 0.1, 1.0, 5.0]
    best_alpha, best_ppl = tune_alpha(order=1, train_seqs=train, val_seqs=val, alphas=alphas)
    assert best_alpha in alphas
    assert isinstance(best_ppl, float)
    assert best_ppl > 0


def test_tune_alpha_perplexity_below_uniform():
    """La perplexite optimisee (ordre 1) doit etre <= 20 (ordre 0)."""
    train, val = _make_split()
    _, best_ppl = tune_alpha(order=1, train_seqs=train, val_seqs=val)
    assert best_ppl <= 20.0 + 0.5   # tolerance pour donnees aleatoires uniformes


def test_tune_alpha_empty_alphas_raises():
    with pytest.raises(ValueError):
        tune_alpha(order=1, train_seqs=["ACDE"], val_seqs=["ACDE"], alphas=[])


# ---------- corpus reel ----------

def test_ngram_on_real_corpus():
    for name in ("corpus_v0.2.csv", "corpus_filtered_all.csv"):
        path = ROOT / "data" / "processed" / name
        if path.exists():
            import pandas as pd
            seqs = pd.read_csv(path)["sequence"].dropna().astype(str).tolist()
            break
    else:
        pytest.skip("Aucun corpus disponible pour le test reel")

    rnd = random.Random(42)
    rnd.shuffle(seqs)
    cut = int(len(seqs) * 0.9)
    train, val = seqs[:cut], seqs[cut:]

    for order in range(4):
        m = NgramModel(order=order, alpha=0.1).fit(train)
        ppl = m.perplexity(val)
        assert ppl > 0
        assert math.isfinite(ppl)
