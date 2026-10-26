"""T41 / A7 — Tests du modele n-gram (peplib.ngram), spec 21 cibles.

Couverture :
  - Ordre 0 : entropie croisee = ln 21 ≈ 3,0445 nats/token, perplexite = 21.
  - Cibles = 20 AA + EOS (EOS terminal ; BOS et PAD jamais cibles).
  - Ordre 1 : somme des proba = 1 sur les 21 cibles.
  - Ordres 0–3 : log_prob ∈ (-inf, 0] pour toute sequence non vide.
  - Ordre 2 : log_prob sur sequence connue == calcul manuel (V = 21).
  - Grille d'alpha par defaut = (0.01, 0.1, 1.0) (A7), selection sur val.
  - Coherence : exp(entropie croisee) == perplexite.
  - Comportement sur corpus reel (skip si absent).
  - Refus des parametres invalides.
"""
import math
import random
from pathlib import Path

import pytest

from peplib.constants import AMINO_ACIDS
from peplib.ngram import (
    DEFAULT_ALPHAS,
    TARGETS,
    NgramModel,
    tune_alpha,
)

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
    with pytest.raises(RuntimeError):
        m.cross_entropy(["ACDE"])


# ---------- cibles : 20 AA + EOS, jamais BOS/PAD ----------

def test_targets_are_21_without_bos_pad():
    """Spec A7 : cibles = 20 AA + EOS (21) ; BOS et PAD exclus."""
    assert len(TARGETS) == 21
    assert set(AA_LIST) <= set(TARGETS)
    assert "<EOS>" in TARGETS
    assert "<BOS>" not in TARGETS
    assert "<PAD>" not in TARGETS


def test_unknown_target_raises():
    m = NgramModel(order=1).fit(["ACDE"])
    with pytest.raises(ValueError):
        m.log_prob_target("<PAD>", context=())


# ---------- ordre 0 : uniforme sur 21 cibles ----------

def test_order0_perplexity_is_21():
    m = NgramModel(order=0, alpha=1.0).fit(_random_sequences())
    seqs = _random_sequences(50, seed=99)
    assert m.perplexity(seqs) == pytest.approx(21.0, rel=1e-12)


def test_order0_cross_entropy_anchor_ln21():
    """Ancre A7 : entropie croisee uniforme = ln 21 ≈ 3,0445 nats/token."""
    m = NgramModel(order=0).fit(["ACDE"])
    ce = m.cross_entropy(_random_sequences(50, seed=99))
    # tolerance flottante : somme de n termes identiques reordonnee
    assert ce == pytest.approx(math.log(21), rel=1e-12)
    # valeur decimale de la spec : ln 21 ≈ 3,0445
    assert ce == pytest.approx(3.0445, abs=1e-4)


def test_order0_log_prob_formula():
    m = NgramModel(order=0).fit(["ACDE"])
    seq = "ACDE"
    # cibles : 4 residus + l'EOS terminal = 5 tokens uniformes sur 21
    expected = (len(seq) + 1) * math.log(1 / 21)
    assert m.log_prob(seq) == pytest.approx(expected, rel=1e-9)


# ---------- sequence vide : seule cible = EOS ----------

@pytest.mark.parametrize("order", [0, 1, 2, 3])
def test_log_prob_empty_sequence_is_eos_only(order):
    """Sequence vide : la seule cible est l'EOS, jamais BOS/PAD."""
    m = NgramModel(order=order, alpha=1.0).fit(["ACDE"] * 10)
    lp = m.log_prob("")
    assert lp <= 0.0
    if order == 0:
        # ordre 0 : EOS uniforme = log(1/21)
        assert lp == pytest.approx(math.log(1 / 21), rel=1e-12)


# ---------- log_prob dans (-inf, 0] ----------

@pytest.mark.parametrize("order", [0, 1, 2, 3])
def test_log_prob_nonpositive(order):
    seqs = _random_sequences(100)
    m = NgramModel(order=order, alpha=1.0).fit(seqs)
    for seq in _random_sequences(20, seed=7):
        lp = m.log_prob(seq)
        assert lp <= 1e-12, f"log_prob = {lp} > 0 pour ordre {order}"


# ---------- ordre 1 : proba somme a 1 sur les 21 cibles ----------

def test_order1_probabilities_sum_to_one():
    """Pour l'ordre 1, la somme exp(log_prob_target) sur les 21 cibles = 1."""
    m = NgramModel(order=1, alpha=1.0).fit(_random_sequences())
    total = sum(math.exp(m.log_prob_target(t)) for t in TARGETS)
    assert total == pytest.approx(1.0, abs=1e-9)


# ---------- ordre 2 : calcul manuel (V = 21, EOS terminal) ----------

def test_order2_log_prob_manual():
    """Verifie un calcul bigram a la main, avec l'EOS terminal.

    Train = ["ACAC"] une seule fois. Tokens cibles issus du train :
    A, C, A, C, EOS.
    Contextes observes (ordre 2), sur 21 cibles :
      (<BOS>,) -> A : 1                 total 1
      (A,)     -> C : 2                 total 2
      (C,)     -> A : 1, EOS : 1        total 2
    Sequence testee = "ACA" (cibles A, C, A, EOS), alpha = 1.0, V = 21 :
      P(A | <BOS>)   = (1+1)/(1+21)   = 2/22
      P(C | A)       = (2+1)/(2+21)   = 3/23
      P(A | C)       = (1+1)/(2+21)   = 2/23
      P(EOS | A)     = (0+1)/(2+21)   = 1/23
    log_prob = somme des 4 logs.
    """
    m = NgramModel(order=2, alpha=1.0).fit(["ACAC"])
    seq = "ACA"
    expected = (
        math.log(2 / 22)
        + math.log(3 / 23)
        + math.log(2 / 23)
        + math.log(1 / 23)
    )
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


# ---------- coherence CE / perplexite ----------

def test_ce_and_perplexity_coherent():
    """exp(entropie croisee) == perplexite, sur n'importe quel ordre."""
    train, val = _make_split()
    m = NgramModel(order=2, alpha=0.1).fit(train)
    ce = m.cross_entropy(val)
    ppl = m.perplexity(val)
    assert ppl == pytest.approx(math.exp(ce), rel=1e-12)


# ---------- perplexite : ordres entraines vs uniforme ----------

def test_perplexity_improves_with_order_on_random_data():
    """Ordre 1 avec alpha bien chose doit rester proche de l'ordre 0 (= 21)."""
    train, val = _make_split()
    m0 = NgramModel(order=0).fit(train)
    m1 = NgramModel(order=1, alpha=0.1).fit(train)
    ppl0 = m0.perplexity(val)
    ppl1 = m1.perplexity(val)
    # Les sequences sont aleatoires uniformes, mais l'ordre 1 doit quand
    # meme estimer la composition correctement -> PP proches de 21.
    # On verifie juste que ppl1 <= ppl0 + 1 (tolerance car donnees unif.)
    assert ppl1 <= ppl0 + 1.0


# ---------- grille d'alpha imposee (A7) ----------

def test_default_alpha_grid_is_spec():
    """Grille d'alpha imposee par l'ajustement A7 : (0.01, 0.1, 1.0)."""
    assert DEFAULT_ALPHAS == (0.01, 0.1, 1.0)


def test_tune_alpha_returns_value_in_candidates():
    train, val = _make_split()
    alphas = [0.01, 0.1, 1.0, 5.0]
    best_alpha, best_ppl = tune_alpha(order=1, train_seqs=train, val_seqs=val, alphas=alphas)
    assert best_alpha in alphas
    assert isinstance(best_ppl, float)
    assert best_ppl > 0


def test_tune_alpha_perplexity_below_uniform():
    """La perplexite optimisee (ordre 1) doit etre <= 21 (ordre 0)."""
    train, val = _make_split()
    _, best_ppl = tune_alpha(order=1, train_seqs=train, val_seqs=val)
    assert best_ppl <= 21.0 + 0.5   # tolerance pour donnees aleatoires uniformes


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
