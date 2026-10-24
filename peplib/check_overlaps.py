"""T38 : chevauchements exacts entre les sources du corpus.

Trois atomes (UniProt, APD3, AMPChallenge) : trois comparaisons deux a
deux, plus le triple recouvrement. Les quasi-doublons relevent de T42
(clusters) et ne sont pas traites ici.

Sorties :
- `data/processed/overlap_t38.csv` : 3 paires ;
- `data/processed/overlap_t38_triple.csv` : le triple recouvrement ;
- `results/overlap_t38.png` : figure en barres empilees.
"""

from itertools import combinations
from pathlib import Path

import pandas as pd

from peplib.build_corpus_v02 import ATOMES, SOURCES


def _ensembles(processed_dir="data/processed"):
    """Construit un ensemble de sequences par atome de source."""
    processed_dir = Path(processed_dir)
    ensembles = {atome: set() for atome in ATOMES}

    for filename, atome in SOURCES.items():
        df = pd.read_csv(processed_dir / filename, dtype=str)
        sequences = df["sequence"].dropna().str.strip().str.upper()
        ensembles[atome] |= set(sequences)

    return ensembles


def compute_overlaps(processed_dir="data/processed"):
    """Trois comparaisons deux a deux, ecrites dans overlap_t38.csv."""
    ensembles = _ensembles(processed_dir)
    lignes = []

    for atome_a, atome_b in combinations(ATOMES, 2):
        a = ensembles[atome_a]
        b = ensembles[atome_b]
        intersection = len(a & b)
        union = len(a | b)
        lignes.append(
            {
                "source_a": atome_a,
                "source_b": atome_b,
                "unique_a": len(a),
                "unique_b": len(b),
                "overlap": intersection,
                "jaccard": round(intersection / union, 6) if union else 0.0,
            }
        )

    resultat = pd.DataFrame(lignes)
    resultat.to_csv(Path(processed_dir) / "overlap_t38.csv", index=False)
    return resultat


def compute_triple(processed_dir="data/processed"):
    """Triple recouvrement, ecrit dans overlap_t38_triple.csv."""
    ensembles = _ensembles(processed_dir)
    intersection = set.intersection(*ensembles.values())
    union = set.union(*ensembles.values())
    ligne = {
        "sources": "+".join(ATOMES),
        "overlap": len(intersection),
        "union": len(union),
        "jaccard": round(len(intersection) / len(union), 6) if union else 0.0,
    }
    pd.DataFrame([ligne]).to_csv(
        Path(processed_dir) / "overlap_t38_triple.csv", index=False
    )
    return ligne


def parts_par_atome(processed_dir="data/processed"):
    """Pour chaque atome : uniques, partages par deux, presents dans les trois."""
    ensembles = _ensembles(processed_dir)
    tous = set.union(*ensembles.values())
    appartenance = {
        sequence: sum(sequence in ens for ens in ensembles.values())
        for sequence in tous
    }

    parts = {}
    for atome, ens in ensembles.items():
        parts[atome] = {
            "unique": sum(1 for s in ens if appartenance[s] == 1),
            "deux": sum(1 for s in ens if appartenance[s] == 2),
            "trois": sum(1 for s in ens if appartenance[s] == 3),
        }
    return parts


def plot_overlaps(processed_dir="data/processed", sortie="results/overlap_t38.png"):
    """Barres empilees : uniques / partages par deux / presents dans les trois."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    parts = parts_par_atome(processed_dir)
    atomes = list(parts)
    uniques = np.array([parts[a]["unique"] for a in atomes])
    deux = np.array([parts[a]["deux"] for a in atomes])
    trois = np.array([parts[a]["trois"] for a in atomes])

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.bar(atomes, uniques, label="uniques a la source")
    ax.bar(atomes, deux, bottom=uniques, label="partages par deux sources")
    ax.bar(atomes, trois, bottom=uniques + deux, label="presents dans les trois")
    ax.set_ylabel("nombre de sequences")
    ax.set_title("T38 - recouvrement exact des trois sources")
    ax.legend()
    fig.tight_layout()

    chemin = Path(sortie)
    chemin.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(chemin, dpi=150)
    plt.close(fig)
    return chemin


if __name__ == "__main__":
    print(compute_overlaps().to_string(index=False))
    print()
    print(compute_triple())
    print()
    print(plot_overlaps())
