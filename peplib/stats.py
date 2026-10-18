import glob
import os
from typing import Dict, List
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

matplotlib.use("Agg")

AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")


def compute_composition(sequences: List[str] | pd.Series) -> Dict[str, float]:
    """T16 — Calcule la composition f(a) = compte(a) / N_total sur les 20 acides aminés.

    Chaque séquence est normalisée en majuscules avant comptage ; les
    caractères hors alphabet sont ignorés (ils ne contribuent ni au
    numérateur ni au dénominateur).
    """
    counts = {aa: 0 for aa in AA_ORDER}
    total = 0
    for seq in sequences:
        for ch in str(seq).upper():
            if ch in counts:
                counts[ch] += 1
                total += 1
    if total == 0:
        raise ValueError("Aucune séquence valide pour calculer la composition.")
    comp = {aa: counts[aa] / total for aa in AA_ORDER}
    assert (
        abs(sum(comp.values()) - 1.0) < 1e-12
    ), "La somme des fréquences doit valoir 1."
    return comp


def kl_divergence(p: np.ndarray, q: np.ndarray) -> float:
    """Calcule KL(p || q) en nats en masquant exactement les entrées où p == 0.

    Retourne inf si q == 0 pour un ou plusieurs indices où p > 0.
    Accepte des distributions non normalisées (le calcul est alors formel).
    """
    mask = p > 0
    with np.errstate(divide="ignore", invalid="ignore"):
        return float(np.sum(p[mask] * np.log(p[mask] / q[mask])))


def jensen_shannon_divergence(
    p_dict: Dict[str, float], q_dict: Dict[str, float]
) -> float:
    """T17 — Calcule JS(p, q) = 0.5 * KL(p || m) + 0.5 * KL(q || m) en nats (bornée par ln 2)."""
    p = np.array([p_dict[aa] for aa in AA_ORDER], dtype=np.float64)
    q = np.array([q_dict[aa] for aa in AA_ORDER], dtype=np.float64)
    m = 0.5 * (p + q)
    js = 0.5 * kl_divergence(p, m) + 0.5 * kl_divergence(q, m)
    return max(0.0, float(js))


def compute_js_matrix(
    compositions: Dict[str, Dict[str, float]],
) -> pd.DataFrame:
    """T17 — Construit la matrice N x N de divergence de Jensen-Shannon entre sources."""
    sources = list(compositions.keys())
    n = len(sources)
    mat = np.zeros((n, n), dtype=np.float64)
    for i, s1 in enumerate(sources):
        for j, s2 in enumerate(sources):
            mat[i, j] = jensen_shannon_divergence(
                compositions[s1], compositions[s2]
            )
    return pd.DataFrame(mat, index=sources, columns=sources)


def compute_entropy(comp_dict: Dict[str, float]) -> float:
    """T17 — Entropie de composition H = - sum f(a) log f(a) en nats."""
    p = np.array([comp_dict[aa] for aa in AA_ORDER], dtype=np.float64)
    mask = p > 0
    return float(-np.sum(p[mask] * np.log(p[mask])))


def plot_length_histograms(
    dfs: Dict[str, pd.DataFrame],
    out_path: str = "results/length_histograms_t16.png",
) -> str:
    """T16 — Génère et sauvegarde les histogrammes de longueurs par source."""
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    fig, axes = plt.subplots(1, len(dfs), figsize=(5 * len(dfs), 4), sharey=False)
    if len(dfs) == 1:
        axes = [axes]
    for ax, (source, df) in zip(axes, dfs.items()):
        lengths = df["sequence"].dropna().astype(str).str.len()
        ax.hist(lengths, bins=range(10, 52), edgecolor="black", alpha=0.75)
        ax.set_title(f"Longueurs - {source} (n={len(df)})")
        ax.set_xlabel("Longueur (AA)")
        ax.set_ylabel("Nombre de séquences")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


def run_stats_pipeline(
    processed_dir: str = "data/processed", results_dir: str = "results"
) -> pd.DataFrame:
    os.makedirs(results_dir, exist_ok=True)
    sources = ["uniprot", "amp", "bactibase"]
    dfs = {
        s: pd.read_csv(os.path.join(processed_dir, f"{s}.csv")) for s in sources
    }
    compositions = {
        s: compute_composition(df["sequence"]) for s, df in dfs.items()
    }
    pd.DataFrame(compositions).to_csv(
        os.path.join(processed_dir, "composition_t16.csv")
    )
    js_df = compute_js_matrix(compositions)
    js_df.to_csv(os.path.join(processed_dir, "js_divergence_t17.csv"))
    plot_length_histograms(
        dfs, out_path=os.path.join(results_dir, "length_histograms_t16.png")
    )
    return js_df


if __name__ == "__main__":
    print(run_stats_pipeline())