
"""T38 : chevauchements exacts entre les sources du corpus."""

from pathlib import Path
from itertools import combinations

import pandas as pd

from peplib.build_corpus_v02 import SOURCES


def compute_overlaps(processed_dir="data/processed"):
    processed_dir = Path(processed_dir)
    sequences = {}

    for name, filename in SOURCES.items():
        df = pd.read_csv(processed_dir / filename, dtype=str)
        sequences[name] = set(
            df["sequence"].dropna().str.strip().str.upper()
        )

    rows = []

    for source_a, source_b in combinations(sequences, 2):
        a = sequences[source_a]
        b = sequences[source_b]

        intersection = len(a & b)
        union = len(a | b)

        rows.append({
            "source_a": source_a,
            "source_b": source_b,
            "unique_a": len(a),
            "unique_b": len(b),
            "overlap": intersection,
            "jaccard": intersection / union if union else 0.0,
        })

    result = pd.DataFrame(rows)
    result.to_csv(processed_dir / "overlap_t38.csv", index=False)
    return result


if __name__ == "__main__":
    result = compute_overlaps()
    print(result.to_string(index=False))
