#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Regenerer la figure 4_pipeline_s1_s2_s3.png avec texte non superpose."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_PNG = REPO_ROOT / "results" / "week3" / "visuals" / "4_pipeline_s1_s2_s3.png"

# Couleurs par etape
COLORS = {
    "S1": "#4c72b0",      # bleu
    "S2": "#dd8452",      # orange
    "S3": "#55a868",      # vert
    "arrow": "#555555",
    "text": "#222222",
}

fig, ax = plt.subplots(figsize=(14, 5))
ax.set_xlim(0, 14)
ax.set_ylim(0, 5)
ax.axis("off")

# ---- Etape S1 : Corpus + peplib v0.1 ----
x1 = 1.0
ax.add_patch(plt.Rectangle((x1, 1.5), 3.0, 2.0, facecolor=COLORS["S1"], alpha=0.15, edgecolor=COLORS["S1"], linewidth=2))
ax.text(x1 + 1.5, 3.2, "SEMAINE 1", ha="center", va="center", fontsize=14, fontweight="bold", color=COLORS["S1"])
ax.text(x1 + 1.5, 2.7, "Corpus + peplib v0.1", ha="center", va="center", fontsize=11, color=COLORS["text"])

s1_items = [
    "T1-T10  : Specs, repo, CI",
    "T11-T20 : Constants, filters, stats",
    "T21-T29 : PWM, descriptors, dedup",
    "Livrable: peplib/ + tests verts",
]
for i, item in enumerate(s1_items):
    ax.text(x1 + 0.15, 2.2 - i * 0.35, item, ha="left", va="center", fontsize=9, color=COLORS["text"])

# Fleche S1 -> S2
ax.annotate("", xy=(5.3, 2.5), xytext=(4.2, 2.5),
            arrowprops=dict(arrowstyle="->", lw=2, color=COLORS["arrow"]))

# ---- Etape S2 : Align + Splits v0.2 ----
x2 = 5.5
ax.add_patch(plt.Rectangle((x2, 1.5), 3.5, 2.0, facecolor=COLORS["S2"], alpha=0.15, edgecolor=COLORS["S2"], linewidth=2))
ax.text(x2 + 1.75, 3.2, "SEMAINE 2", ha="center", va="center", fontsize=14, fontweight="bold", color=COLORS["S2"])
ax.text(x2 + 1.75, 2.7, "Align + Splits v0.2", ha="center", va="center", fontsize=11, color=COLORS["text"])

s2_items = [
    "A1-A5  : NW/SW, BLOSUM62",
    "A6-A8  : Union-find 90%, splits 90/5/5",
    "A9-A10 : Fuite F(0.90)=0 verifiee",
    "Livrable: data/splits/*_v02.csv verrouilles",
]
for i, item in enumerate(s2_items):
    ax.text(x2 + 0.15, 2.2 - i * 0.35, item, ha="left", va="center", fontsize=9, color=COLORS["text"])

# Fleche S2 -> S3
ax.annotate("", xy=(10.3, 2.5), xytext=(9.2, 2.5),
            arrowprops=dict(arrowstyle="->", lw=2, color=COLORS["arrow"]))

# ---- Etape S3 : Modeles + Generation ----
x3 = 10.5
ax.add_patch(plt.Rectangle((x3, 1.5), 3.0, 2.0, facecolor=COLORS["S3"], alpha=0.15, edgecolor=COLORS["S3"], linewidth=2))
ax.text(x3 + 1.5, 3.2, "SEMAINE 3", ha="center", va="center", fontsize=14, fontweight="bold", color=COLORS["S3"])
ax.text(x3 + 1.5, 2.7, "Modeles + Generation", ha="center", va="center", fontsize=11, color=COLORS["text"])

s3_items = [
    "T52-T55 : LSTM / Transformer causal",
    "T56-T57 : Entrainement 6 runs (3 seeds x 2)",
    "T58-T60 : Generation + prox. SW + eval",
    "Livrable: val_pp < 18.65 (gate OK)",
]
for i, item in enumerate(s3_items):
    ax.text(x3 + 0.15, 2.2 - i * 0.35, item, ha="left", va="center", fontsize=9, color=COLORS["text"])

# Titre general
ax.text(7.0, 4.5, "Pipeline Semaines 1 \u2192 2 \u2192 3", ha="center", va="center",
        fontsize=16, fontweight="bold", color=COLORS["text"])

# Legende en bas
ax.text(7.0, 0.5,
        "Chaque etape produit des artefacts verrouilles (code, donnees, config) -> reproductibilite totale",
        ha="center", va="center", fontsize=10, style="italic", color="#666666")

plt.tight_layout()
OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(OUT_PNG, dpi=200, bbox_inches="tight")
plt.close()
print("Figure regeneree :", OUT_PNG.relative_to(REPO_ROOT))