"""
PWM alignée depuis les extrémités N et C, sans troncature (T20).

Formules du spec (Semaine 1, §T20) :

- PWM(i, a) = log2[ p̂ᵢ(a) / p_bg(a) ], p_bg = composition globale.
- H(i) = - Σ_a p̂ᵢ(a) · log2 p̂ᵢ(a), avec 0 · log 0 = 0 (masquage exact).
- I(i) = log2(20) - [ H(i) + e(nᵢ) ],
  e(n) = (s - 1) / (2 · ln 2 · n), s = 20, n = nᵢ
  (correction de petit échantillon, Schneider & Stephens 1990).

Conventions :
- Base log2 pour la PWM, l'entropie et l'information ; la formule e(n) est
  explicitement en nats (2 · ln 2 au dénominateur) — ne pas unifier.
- Aucune troncature : nᵢ = nombre de séquences couvrant la position i et
  varie avec i.
- Entrées nulles de la PWM masquées (NaN), jamais -inf ; pas de
  pseudo-count.
- Les séquences doivent être nettoyées en amont (alphabet canonique) ;
  tout caractère hors alphabet lève ValueError.
"""

import math
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd

from peplib.stats import AA_ORDER, compute_composition

S_ALPHABET = 20  # taille de l'alphabet (s du spec)
LOG2_S = math.log2(S_ALPHABET)  # log2 20 = 4.3219...
LN_2 = math.log(2.0)


def small_sample_correction(n) -> np.ndarray:
    """Correction de petit échantillon e(n) = (s-1)/(2·ln2·n), s = 20.

    Schneider & Stephens (1990). Voir docstring du module.
    """
    n = np.asarray(n, dtype=np.float64)
    if np.any(n <= 0):
        raise ValueError("n doit être strictement positif.")
    return (S_ALPHABET - 1) / (2.0 * LN_2 * n)


def build_counts(sequences, terminus: str = "N") -> Tuple[np.ndarray, np.ndarray]:
    """Compte les résidus par position depuis N ou depuis C.

    Position 0 = premier résidu (N) ou dernier résidu (C). Aucune
    troncature : le tableau s'étend jusqu'à la plus longue séquence et
    nᵢ = nombre de séquences couvrant la position i.

    Retourne (counts, n) : counts de forme (n_positions, 20) dans
    l'ordre AA_ORDER, n de forme (n_positions,).
    """
    if terminus not in ("N", "C"):
        raise ValueError("terminus doit valoir 'N' ou 'C'.")

    seqs: List[str] = []
    for seq in sequences:
        if seq is None or (isinstance(seq, float) and math.isnan(seq)):
            continue
        s = str(seq).upper()
        for ch in s:
            if ch not in AA_ORDER:
                raise ValueError(f"Caractère hors alphabet : {ch!r}")
        seqs.append(s)

    if not seqs:
        raise ValueError("Aucune séquence valide pour la PWM.")

    max_len = max(len(s) for s in seqs)
    counts = np.zeros((max_len, len(AA_ORDER)), dtype=np.float64)
    col = {aa: j for j, aa in enumerate(AA_ORDER)}

    for s in seqs:
        n = len(s)
        for i in range(n):
            j = i if terminus == "N" else n - 1 - i
            counts[i, col[s[j]]] += 1.0

    n_per_pos = counts.sum(axis=1)
    return counts, n_per_pos


def pwm_matrix(sequences, terminus: str = "N", background=None) -> pd.DataFrame:
    """PWM(i, a) = log2[ p̂ᵢ(a) / p_bg(a) ].

    ``background`` : dict AA -> probabilité ; par défaut, composition
    globale des séquences fournies (reprise de peplib.stats). Un AA
    absent du dict a un fond nul. Les entrées où p̂ᵢ(a) = 0 ou où
    p_bg(a) = 0 sont masquées (NaN), jamais -inf ni inf.
    """
    if background is None:
        background = compute_composition(sequences)
    bg = np.array([background.get(aa, 0.0) for aa in AA_ORDER], dtype=np.float64)

    counts, n = build_counts(sequences, terminus)
    p_hat = counts / n[:, None]

    with np.errstate(divide="ignore", invalid="ignore"):
        pwm = np.log2(p_hat / bg)

    masque = (p_hat > 0) & (bg > 0)  # masquage exact des zéros
    pwm = np.where(masque, pwm, np.nan)
    return pd.DataFrame(pwm, columns=list(AA_ORDER))


def position_entropy(sequences, terminus: str = "N") -> pd.Series:
    """H(i) = - Σ_a p̂ᵢ(a) · log2 p̂ᵢ(a), 0 · log 0 = 0."""
    counts, n = build_counts(sequences, terminus)
    p_hat = counts / n[:, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        terms = p_hat * np.log2(p_hat)
    terms = np.where(p_hat > 0, terms, 0.0)  # 0 · log 0 = 0
    h = -terms.sum(axis=1)
    return pd.Series(h, index=range(len(h)), name="H")


def position_information(sequences, terminus: str = "N") -> pd.Series:
    """I(i) = log2(20) - [ H(i) + e(nᵢ) ], correction Schneider & Stephens."""
    h = position_entropy(sequences, terminus)
    _, n = build_counts(sequences, terminus)
    info = LOG2_S - (h.to_numpy() + small_sample_correction(n))
    return pd.Series(info, index=h.index, name="I")


def run_pwm_pipeline(
    processed_dir: str = "data/processed", results_dir: str = "results"
) -> Dict[str, pd.DataFrame]:
    """Calcule les PWM et informations N et C sur le corpus filtré.

    Lit ``{uniprot,amp,bactibase}.csv`` et sauvegarde dans ``results/`` :
    ``pwm_{N,C}.csv``, ``information_{N,C}.csv``. Retourne les matrices.
    """
    import os

    os.makedirs(results_dir, exist_ok=True)
    sources = ["uniprot", "amp", "bactibase"]
    all_seqs: List[str] = []
    for s in sources:
        df = pd.read_csv(os.path.join(processed_dir, f"{s}.csv"))
        all_seqs.extend(df["sequence"].dropna().astype(str).tolist())

    out: Dict[str, pd.DataFrame] = {}
    for terminus in ("N", "C"):
        pwm = pwm_matrix(all_seqs, terminus=terminus)
        info = position_information(all_seqs, terminus=terminus)
        pwm.to_csv(os.path.join(results_dir, f"pwm_{terminus}.csv"))
        info.to_csv(
            os.path.join(results_dir, f"information_{terminus}.csv"),
            header=True,
        )
        out[f"pwm_{terminus}"] = pwm
        out[f"information_{terminus}"] = info
    return out


if __name__ == "__main__":
    for name, mat in run_pwm_pipeline().items():
        print(name, mat.shape)
