"""
Figures T21 : distribution des descripteurs, heatmap de la PWM, courbe
d'information par position.

Sorties dans ``results/`` :
- descriptors_distribution_t21.png : boxplots des 6 descripteurs par source ;
- pwm_heatmap_t21.png : heatmap de la PWM N et C (échelle commune) ;
- information_curve_t21.png : courbes d'information N et C.

Backend matplotlib ``Agg`` (aucune fenêtre graphique) ; ``plt.show()``
jamais utilisé.
"""

import os
from typing import Dict

import matplotlib
matplotlib.use("Agg")
import matplotlib.figure  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from peplib.descriptors import (  # noqa: E402
    calculate_boman_index,
    calculate_charge,
    calculate_gravy,
    calculate_hydrophobic_moment,
    calculate_molecular_mass,
    calculate_pi,
)
import matplotlib.colors  # noqa: E402

SOURCES = ["uniprot", "amp", "bactibase"]

DESCRIPTORS = [
    ("charge", "Charge nette (pH 7)"),
    ("pi", "Point isoélectrique pI"),
    ("gravy", "GRAVY (Kyte-Doolittle)"),
    ("mass", "Masse monoisotopique (Da)"),
    ("hydrophobic_moment", "Moment hydrophobe"),
    ("boman", "Indice de Boman (kcal/mol)"),
]

PWM_CSV = {"N": "pwm_N.csv", "C": "pwm_C.csv"}
INFO_CSV = {"N": "information_N.csv", "C": "information_C.csv"}


def load_corpus(processed_dir: str = "data/processed") -> Dict[str, pd.DataFrame]:
    """Lit les CSV traités par source. Returns: dict source -> DataFrame
    avec les colonnes sequence, source, id."""
    out = {}
    for s in SOURCES:
        out[s] = pd.read_csv(os.path.join(processed_dir, f"{s}.csv"))
    return out


def compute_descriptors(df: pd.DataFrame) -> pd.DataFrame:
    """Calcule les 6 descripteurs par séquence. Returns: DataFrame avec
    charge, pi, gravy, mass, hydrophobic_moment, boman. Le piège NaN
    (astype(str) avant dropna) est évité : dropna() d'abord."""
    seqs = df["sequence"].dropna().astype(str).tolist()
    rows = []
    for seq in seqs:
        rows.append(
            {
                "charge": calculate_charge(seq),
                "pi": calculate_pi(seq),
                "gravy": calculate_gravy(seq),
                "mass": calculate_molecular_mass(seq),
                "hydrophobic_moment": calculate_hydrophobic_moment(seq),
                "boman": calculate_boman_index(seq),
            }
        )
    return pd.DataFrame(rows)


def plot_descriptor_distributions(
    descriptors_per_source: Dict[str, pd.DataFrame], out_path: str
) -> matplotlib.figure.Figure:
    """Figure 1 : 6 boxplots (un par descripteur), une boîte par source.

    Returns: la Figure sauvegardée dans ``out_path``.
    """
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    for ax, (col, label) in zip(axes.flat, DESCRIPTORS):
        data = [d[col].to_numpy() for d in descriptors_per_source.values()]
        ax.boxplot(data)
        ax.set_xticks(range(1, len(descriptors_per_source) + 1))
        ax.set_xticklabels(list(descriptors_per_source.keys()))
        ax.set_ylabel(label)
        ax.grid(axis="y", alpha=0.3)
    fig.suptitle("Distribution des descripteurs par source (T21)")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return fig


def plot_pwm_heatmap(pwm_n: pd.DataFrame, pwm_c: pd.DataFrame, out_path: str) -> matplotlib.figure.Figure:
    """Figure 2 : heatmap de la PWM N-terminal et C-terminal, échelle de
    couleurs commune.

    Returns: la Figure sauvegardée dans ``out_path``.
    """
    a = pwm_n.to_numpy(dtype=float)
    b = pwm_c.to_numpy(dtype=float)
    finite = np.concatenate([a[np.isfinite(a)], b[np.isfinite(b)]])
    vmin, vmax = np.min(finite), np.max(finite)
    cmap = matplotlib.colormaps["viridis"].with_extremes(bad="lightgrey")

    fig, axes = plt.subplots(1, 2, figsize=(14, 7), sharey=True)
    for ax, mat, titre in ((axes[0], a, "PWM N-terminal"), (axes[1], b, "PWM C-terminal")):
        im = ax.imshow(np.ma.masked_invalid(mat), aspect="auto",
                       vmin=vmin, vmax=vmax, cmap=cmap)
        ax.set_xticks(range(20))
        ax.set_xticklabels(list(pwm_n.columns), fontsize=7)
        ax.grid(False)
    axes[0].set_ylabel("Position depuis le terminus")
    axes[0].set_title("PWM N-terminal")
    axes[1].set_title("PWM C-terminal")
    fig.colorbar(im, ax=axes, label="log2[ p(i,a) / p_bg(a) ]")
    fig.suptitle("Heatmap de la matrice PWM (T21)", y=0.98)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return fig


def plot_information_curve(
    info_n: pd.Series, info_c: pd.Series, out_path: str
) -> matplotlib.figure.Figure:
    """Figure 3 : courbe d'information par position, terminus N et C.

    Returns: la Figure sauvegardée dans ``out_path``.
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(np.arange(len(info_n)), info_n.to_numpy(), label="N-terminal")
    ax.plot(np.arange(len(info_c)), info_c.to_numpy(), label="C-terminal")
    ax.set_xlabel("Position depuis le terminus")
    ax.set_ylabel("Information (bits)")
    ax.set_title("Information par position (T21)")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return fig


def run_figures_pipeline(
    processed_dir: str = "data/processed", results_dir: str = "results"
) -> Dict[str, str]:
    """Orchestre T21 : corpus -> descripteurs -> PWM (peplib.pwm) ->
    3 figures dans results/. Returns: dict nom -> chemin de sortie."""
    from peplib.pwm import run_pwm_pipeline

    os.makedirs(results_dir, exist_ok=True)
    pwm_mat = run_pwm_pipeline(processed_dir=processed_dir, results_dir=results_dir)

    corpus = load_corpus(processed_dir)
    descriptors_per_source = {s: compute_descriptors(df) for s, df in corpus.items()}

    paths = {
        "descriptors_distribution": os.path.join(
            results_dir, "descriptors_distribution_t21.png"
        ),
        "pwm_heatmap": os.path.join(results_dir, "pwm_heatmap_t21.png"),
        "information_curve": os.path.join(results_dir, "information_curve_t21.png"),
    }
    plot_descriptor_distributions(descriptors_per_source, paths["descriptors_distribution"])
    plot_pwm_heatmap(pwm_mat["pwm_N"], pwm_mat["pwm_C"], paths["pwm_heatmap"])
    plot_information_curve(pwm_mat["information_N"], pwm_mat["information_C"], paths["information_curve"])
    return paths


if __name__ == "__main__":
    for name, path in run_figures_pipeline().items():
        print(f"{name}: {path}")
