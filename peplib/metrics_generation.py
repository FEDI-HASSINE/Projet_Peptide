
"""T59 — Metriques statistiques des peptides generes.

Ce module ne predit pas l'activite biologique des sequences.
Il mesure uniquement leur validite et leurs proprietes statistiques.
"""

from collections import Counter
import math

import numpy as np
from scipy.spatial.distance import jensenshannon
from scipy.stats import ks_2samp

from peplib.constants import AA_LIST


AA_SET = set(AA_LIST)


def is_valid_peptide(sequence: str) -> bool:
    """Vrai si la sequence est non vide et contient seulement les 20 AA."""
    return (
        isinstance(sequence, str)
        and len(sequence) > 0
        and all(aa in AA_SET for aa in sequence)
    )


def validity_rate(sequences: list[str]) -> float:
    """Proportion de peptides valides parmi tous les peptides generes."""
    if not sequences:
        raise ValueError("La liste de sequences est vide")

    valid_count = sum(is_valid_peptide(s) for s in sequences)
    return valid_count / len(sequences)


def uniqueness_rate(sequences: list[str]) -> float:
    """Proportion de sequences distinctes parmi les sequences generees."""
    if not sequences:
        raise ValueError("La liste de sequences est vide")

    # On conserve la representation exacte, sans nettoyer silencieusement.
    if not all(isinstance(s, str) for s in sequences):
        raise ValueError("Toutes les sequences doivent etre des chaines")

    return len(set(sequences)) / len(sequences)


def sequence_length_stats(sequences: list[str]) -> dict:
    """Statistiques de longueur pour les sequences valides."""
    valid = [s for s in sequences if is_valid_peptide(s)]

    if not valid:
        raise ValueError("Aucun peptide valide")

    lengths = np.array([len(s) for s in valid], dtype=float)

    return {
        "mean": float(np.mean(lengths)),
        "std": float(np.std(lengths)),
        "min": int(np.min(lengths)),
        "max": int(np.max(lengths)),
        "median": float(np.median(lengths)),
    }


def amino_acid_distribution(sequences: list[str]) -> np.ndarray:
    """Distribution normalisee des 20 acides amines."""
    valid = [s for s in sequences if is_valid_peptide(s)]

    if not valid:
        raise ValueError("Aucun peptide valide")

    counts = Counter("".join(valid))
    frequencies = np.array(
        [counts.get(aa, 0) for aa in AA_LIST],
        dtype=float,
    )

    return frequencies / frequencies.sum()


def js_divergence(
    generated: list[str],
    reference: list[str],
) -> float:
    """Divergence Jensen-Shannon entre distributions de 20 AA.

    Utilise le logarithme naturel : valeur entre 0 et ln(2).
    """

    p = amino_acid_distribution(generated)
    q = amino_acid_distribution(reference)

    # scipy renvoie la distance JS (racine de la divergence).
    distance = jensenshannon(p, q, base=math.e)
    return float(distance ** 2)


def ks_length_statistic(
    generated: list[str],
    reference: list[str],
) -> float:
    """Statistique KS entre distributions de longueurs des peptides valides."""

    gen_lengths = [
        len(s) for s in generated if is_valid_peptide(s)
    ]
    ref_lengths = [
        len(s) for s in reference if is_valid_peptide(s)
    ]

    if not gen_lengths or not ref_lengths:
        raise ValueError("Chaque groupe doit contenir un peptide valide")

    result = ks_2samp(gen_lengths, ref_lengths)
    return float(result.statistic)


def generation_metrics(
    generated: list[str],
    reference: list[str],
) -> dict:
    """Calcule les metriques principales de generation."""

    if not generated:
        raise ValueError("Aucune sequence generee")

    stats = sequence_length_stats(generated)

    return {
        "n_generated": len(generated),
        "n_valid": sum(is_valid_peptide(s) for s in generated),
        "validity_rate": validity_rate(generated),
        "uniqueness_rate": uniqueness_rate(generated),
        "length_mean": stats["mean"],
        "length_std": stats["std"],
        "length_min": stats["min"],
        "length_max": stats["max"],
        "length_median": stats["median"],
        "js_divergence": js_divergence(generated, reference),
        "ks_length": ks_length_statistic(generated, reference),
    }
