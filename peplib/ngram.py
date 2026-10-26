"""
Modele n-gram pour les peptides (A7 / T41).

Ordres supportes : 0, 1, 2, 3.
  - Ordre 0 : distribution uniforme sur les 20 AA -> perplexite = 20.
  - Ordre 1 : unigram P(aa), lissage add-alpha.
  - Ordre 2 : bigram  P(aa | aa_{-1}), lissage add-alpha.
  - Ordre 3 : trigram P(aa | aa_{-2}, aa_{-1}), lissage add-alpha.

Lissage add-alpha (Laplace generalise) :
    P(aa | contexte) = (N(contexte, aa) + alpha) / (N(contexte) + alpha * V)
    V = 20 (taille de l'alphabet canonique).

Le contexte manquant en debut de sequence est rempli par un sentinel '<BOS>'
(ordre 2 : 1 sentinel ; ordre 3 : 2 sentinels). Les sentinels ne sont jamais
produits ni comptes comme des residus valides : ils servent uniquement a
initialiser le contexte.

Utilisation type :
    model = NgramModel(order=2, alpha=0.1).fit(train_seqs)
    ppl = model.perplexity(val_seqs)

    best_alpha, best_ppl = tune_alpha(order=2, train_seqs=train, val_seqs=val)
"""

import math
from collections import defaultdict
from typing import Sequence

from peplib.constants import AMINO_ACIDS

VOCAB_SIZE: int = len(AMINO_ACIDS)   # 20
_BOS = "<BOS>"                        # sentinel de debut de sequence


class NgramModel:
    """Modele n-gram sur les sequences de peptides.

    Parameters
    ----------
    order : int
        0 = uniforme, 1 = unigram, 2 = bigram, 3 = trigram.
    alpha : float
        Parametre de lissage add-alpha (>= 0). alpha=0 = maximum de
        vraisemblance pur (attention : log 0 possible sur sequences inedites).
    """

    def __init__(self, order: int, alpha: float = 1.0) -> None:
        if order not in (0, 1, 2, 3):
            raise ValueError(f"order doit etre 0, 1, 2 ou 3 (recu : {order})")
        if not isinstance(alpha, (int, float)) or isinstance(alpha, bool):
            raise TypeError("alpha doit etre un nombre reel")
        if alpha < 0:
            raise ValueError("alpha doit etre >= 0")
        self.order = order
        self.alpha = float(alpha)
        # _counts[context_tuple][aa] = int
        self._counts: dict[tuple, dict[str, int]] = defaultdict(
            lambda: defaultdict(int)
        )
        self._fitted: bool = False

    # ------------------------------------------------------------------
    # Ajustement
    # ------------------------------------------------------------------

    def fit(self, sequences: Sequence[str]) -> "NgramModel":
        """Compte les n-grams sur les sequences d'entrainement.

        Chaque sequence est traitee independamment ; les sentinels BOS
        sont preposes pour initialiser le contexte (ils ne contribuent
        jamais au numerateur, seulement au denominateur des contextes
        vus en debut de sequence).
        """
        self._counts.clear()

        if self.order == 0:
            # Pas de comptage necessaire : distribution uniforme pure.
            self._fitted = True
            return self

        for seq in sequences:
            if not seq:
                continue
            # Prepend (order - 1) sentinels pour les contextes de debut.
            context_window = [_BOS] * (self.order - 1)
            for aa in seq:
                context = tuple(context_window[-(self.order - 1):]) if self.order > 1 else ()
                self._counts[context][aa] += 1
                # Faire glisser la fenetre.
                context_window.append(aa)

        self._fitted = True
        return self

    # ------------------------------------------------------------------
    # Log-vraisemblance
    # ------------------------------------------------------------------

    def _log_prob_aa(self, aa: str, context: tuple) -> float:
        """Log-probabilite (base e) de l'acide amine aa etant donne le contexte."""
        if self.order == 0:
            return math.log(1.0 / VOCAB_SIZE)

        ctx_counts = self._counts.get(context, {})
        aa_count = ctx_counts.get(aa, 0)
        ctx_total = sum(ctx_counts.values())

        # add-alpha smoothing
        prob = (aa_count + self.alpha) / (ctx_total + self.alpha * VOCAB_SIZE)
        return math.log(prob)

    def log_prob(self, sequence: str) -> float:
        """Log-vraisemblance (base e) de toute la sequence.

        Retourne 0.0 pour une sequence vide.
        Leve RuntimeError si le modele n'a pas ete ajuste.
        """
        if not self._fitted:
            raise RuntimeError("NgramModel.fit() doit etre appele avant log_prob()")
        if not sequence:
            return 0.0

        total = 0.0
        context_window = [_BOS] * (self.order - 1)

        for aa in sequence:
            context = tuple(context_window[-(self.order - 1):]) if self.order > 1 else ()
            total += self._log_prob_aa(aa, context)
            context_window.append(aa)

        return total

    # ------------------------------------------------------------------
    # Perplexite
    # ------------------------------------------------------------------

    def perplexity(self, sequences: Sequence[str]) -> float:
        """Perplexite du modele sur un ensemble de sequences.

        PP = exp( -sum_log_prob / N_total_residus )

        Les sequences vides sont ignorees dans le decompte des residus
        (elles contribuent 0 au log-prob et 0 au total de tokens).
        Leve ValueError si aucun residu n'est present dans sequences.
        """
        if not self._fitted:
            raise RuntimeError("NgramModel.fit() doit etre appele avant perplexity()")

        total_log_prob = 0.0
        total_tokens = 0

        for seq in sequences:
            if not seq:
                continue
            total_log_prob += self.log_prob(seq)
            total_tokens += len(seq)

        if total_tokens == 0:
            raise ValueError("Aucun residu dans les sequences fournies.")

        return math.exp(-total_log_prob / total_tokens)

    # ------------------------------------------------------------------
    # Representation
    # ------------------------------------------------------------------

    def __repr__(self) -> str:  # pragma: no cover
        status = "fitted" if self._fitted else "not fitted"
        return f"NgramModel(order={self.order}, alpha={self.alpha}, {status})"


# ----------------------------------------------------------------------
# Optimisation de alpha sur l'ensemble de validation
# ----------------------------------------------------------------------

_DEFAULT_ALPHAS = [1e-4, 1e-3, 0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]


def tune_alpha(
    order: int,
    train_seqs: Sequence[str],
    val_seqs: Sequence[str],
    alphas: list[float] | None = None,
) -> tuple[float, float]:
    """Selectionne le meilleur alpha (perplexite minimale sur val).

    Parameters
    ----------
    order : int
        Ordre du modele n-gram (0–3).
    train_seqs : list[str]
        Sequences d'entrainement.
    val_seqs : list[str]
        Sequences de validation.
    alphas : list[float] | None
        Valeurs candidates ; _DEFAULT_ALPHAS par defaut.

    Returns
    -------
    (best_alpha, best_perplexity) : tuple[float, float]
        La valeur d'alpha qui minimise la perplexite sur val et la
        perplexite correspondante.
    """
    if alphas is None:
        alphas = _DEFAULT_ALPHAS

    if not alphas:
        raise ValueError("La liste alphas ne peut pas etre vide.")

    best_alpha: float = alphas[0]
    best_ppl: float = float("inf")

    for alpha in alphas:
        model = NgramModel(order=order, alpha=alpha)
        model.fit(train_seqs)
        ppl = model.perplexity(val_seqs)
        if ppl < best_ppl:
            best_ppl = ppl
            best_alpha = alpha

    return best_alpha, best_ppl


# ----------------------------------------------------------------------
# Point d'entree : evaluation sur corpus v0.2
# ----------------------------------------------------------------------

def run_ngram_eval(
    corpus_path: str = "data/processed/corpus_v0.2.csv",
    val_fraction: float = 0.1,
    seed: int = 42,
) -> None:
    """Evalue les 4 ordres (alpha optimise) sur le corpus.

    Le corpus est partitionne en train/val (1 - val_fraction / val_fraction)
    de facon reproductible. Utilise les splits de Sofiene (A6) si disponibles.
    """
    import random
    from pathlib import Path

    import pandas as pd

    # Charger les sequences
    splits_path = Path("data/processed/splits")
    if splits_path.exists():
        train_path = splits_path / "train.csv"
        val_path = splits_path / "val.csv"
        if train_path.exists() and val_path.exists():
            train_seqs = pd.read_csv(train_path)["sequence"].dropna().astype(str).tolist()
            val_seqs = pd.read_csv(val_path)["sequence"].dropna().astype(str).tolist()
        else:
            raise FileNotFoundError("Splits A6 incomplets dans data/processed/splits/")
    else:
        # Repli : split interne reproductible
        all_seqs = pd.read_csv(corpus_path)["sequence"].dropna().astype(str).tolist()
        rnd = random.Random(seed)
        rnd.shuffle(all_seqs)
        cut = int(len(all_seqs) * (1 - val_fraction))
        train_seqs, val_seqs = all_seqs[:cut], all_seqs[cut:]

    print(f"Train : {len(train_seqs)} sequences | Val : {len(val_seqs)} sequences\n")

    for order in range(4):
        best_alpha, best_ppl = tune_alpha(order, train_seqs, val_seqs)
        print(
            f"Ordre {order} | best_alpha={best_alpha:.4g} | "
            f"perplexite_val={best_ppl:.4f}"
        )


if __name__ == "__main__":  # pragma: no cover
    run_ngram_eval()
