"""Tracabilite des runs (Semaine 3, correction P0-1).

Tout nombre publie doit pouvoir etre relie a :

- un commit git (SHA court + arbre propre ou non) ;
- un hash de configuration (SHA-256 du YAML canonicalise) ;
- un environnement d'execution (Python, PyTorch, device, GPU) ;
- un horodatage UTC ;
- un hash SHA-256 des checkpoints utilises.

Ce module ne contient aucune logique d'entrainement : il ne fait que
collecter, hasher et mettre en forme. Il est utilisable sans GPU.
"""

from __future__ import annotations

import hashlib
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent

WEEK3_START = "<!-- semaine3:start -->"
WEEK3_END = "<!-- semaine3:end -->"


# ----------------------------------------------------------------------
# Identite du code et des artefacts
# ----------------------------------------------------------------------

def git_commit(repo_dir: str | Path = REPO_ROOT, length: int = 12) -> str:
    """SHA du commit HEAD, tronque ; 'unknown' hors depot git."""
    try:
        out = subprocess.check_output(
            ["git", "rev-parse", f"--short={length}", "HEAD"],
            cwd=str(repo_dir),
            stderr=subprocess.DEVNULL,
            text=True,
        )
        return out.strip() or "unknown"
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def git_is_dirty(repo_dir: str | Path = REPO_ROOT) -> bool:
    """True si l'arbre de travail contient des modifications non commitees."""
    try:
        out = subprocess.check_output(
            ["git", "status", "--porcelain"],
            cwd=str(repo_dir),
            stderr=subprocess.DEVNULL,
            text=True,
        )
        return bool(out.strip())
    except (OSError, subprocess.CalledProcessError):
        return True


def sha256_file(path: str | Path, chunk_size: int = 1 << 20) -> str:
    """SHA-256 hexadecinal d'un fichier, lu par blocs."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def config_sha256(config: dict) -> str:
    """SHA-256 d'une config YAML, insensible a l'ordre des cles."""
    canonical = yaml.safe_dump(config, sort_keys=True, allow_unicode=True)
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def environment_info() -> dict:
    """Versions et device effectivement utilises."""
    import torch  # import local : les tests purs n'ont pas besoin de CUDA

    cuda = torch.cuda.is_available()
    return {
        "python_version": platform.python_version(),
        "torch_version": torch.__version__,
        "cuda_available": cuda,
        "gpu_name": torch.cuda.get_device_name(0) if cuda else None,
        "platform": platform.platform(),
    }


def utc_timestamp() -> str:
    """Horodatage UTC ISO-8601 a la seconde."""
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def run_metadata(
    config: dict,
    seed: int,
    device: str,
    repo_dir: str | Path = REPO_ROOT,
) -> dict:
    """Metadonnees communes a toutes les lignes d'un run.

    Ces colonnes sont ecrites dans results/week3/runs.csv a chaque epoch :
    une ligne isolee suffit donc a savoir d'ou elle vient.
    """
    env = environment_info()
    return {
        "timestamp_utc": utc_timestamp(),
        "git_sha": git_commit(repo_dir),
        "git_dirty": git_is_dirty(repo_dir),
        "config_sha256": config_sha256(config),
        "seed": seed,
        "device": device,
        "torch_version": env["torch_version"],
        "python_version": env["python_version"],
        "gpu_name": env["gpu_name"] or "",
    }


# ----------------------------------------------------------------------
# Checkpoints
# ----------------------------------------------------------------------

def checkpoint_manifest(
    checkpoint_root: str | Path = "checkpoints",
    pattern: str = "**/*.pt",
) -> pd.DataFrame:
    """Table (checkpoint, sha256, size_bytes) des checkpoints presents."""
    root = Path(checkpoint_root)
    rows = []
    for path in sorted(root.glob(pattern)):
        if not path.is_file():
            continue
        rows.append(
            {
                "checkpoint": path.relative_to(root).as_posix(),
                "sha256": sha256_file(path),
                "size_bytes": path.stat().st_size,
            }
        )
    return pd.DataFrame(rows, columns=["checkpoint", "sha256", "size_bytes"])


# ----------------------------------------------------------------------
# PROVENANCE.md
# ----------------------------------------------------------------------

def upsert_markdown_block(
    path: str | Path,
    block: str,
    start_marker: str = WEEK3_START,
    end_marker: str = WEEK3_END,
    title: str = "## Semaine 3 -- Artefacts reproductibles",
) -> None:
    """Insere ou remplace un bloc delimite par des marqueurs HTML.

    Idempotent : deux appels avec le meme bloc laissent le fichier
    identique. Le reste de PROVENANCE.md est preserve.
    """
    path = Path(path)
    existing = path.read_text(encoding="utf-8") if path.exists() else ""
    wrapped = f"{start_marker}\n{block.strip()}\n{end_marker}\n"

    if start_marker in existing and end_marker in existing:
        head, _, rest = existing.partition(start_marker)
        _, _, tail = rest.partition(end_marker)
        new_text = head + wrapped + tail.lstrip("\n")
    else:
        separator = "" if not existing else ("\n" if existing.endswith("\n") else "\n\n")
        new_text = existing + separator + title + "\n\n" + wrapped

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(new_text, encoding="utf-8")


def _format_table(frame: pd.DataFrame, columns: list[str]) -> str:
    """Table markdown minimale, sans dependance a tabulate."""
    if frame.empty:
        return "_aucune ligne_"
    header = "| " + " | ".join(columns) + " |"
    sep = "|" + "|".join("---" for _ in columns) + "|"
    lines = [header, sep]
    for _, row in frame[columns].iterrows():
        cells = ["" if pd.isna(row[c]) else str(row[c]) for c in columns]
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def best_run_table(runs: pd.DataFrame) -> pd.DataFrame:
    """Meilleure epoch de chaque (model, seed), triee."""
    if runs.empty:
        return pd.DataFrame(
            columns=["model", "seed", "epoch", "train_loss", "val_loss", "val_pp"]
        )
    best = runs.loc[runs.groupby(["model", "seed"])["val_loss"].idxmin()].copy()
    best = best.sort_values(["model", "seed"])
    return best.reset_index(drop=True)


def model_summary(runs: pd.DataFrame) -> pd.DataFrame:
    """Moyenne et ecart-type du meilleur val_pp par modele (3 seeds)."""
    best = best_run_table(runs)
    if best.empty:
        return pd.DataFrame(columns=["model", "n_seeds", "val_pp_mean", "val_pp_std"])
    summary = (
        best.groupby("model")["val_pp"]
        .agg(n_seeds="count", val_pp_mean="mean", val_pp_std="std")
        .reset_index()
    )
    summary["val_pp_mean"] = summary["val_pp_mean"].round(4)
    summary["val_pp_std"] = summary["val_pp_std"].round(4)
    return summary


def render_week3_provenance(
    runs: pd.DataFrame,
    generation: pd.DataFrame,
    baselines: dict[str, float],
    manifest_csv: str | Path,
    generated_at: str | None = None,
) -> str:
    """Bloc markdown decrivant les artefacts de la Semaine 3.

    Parameters
    ----------
    runs : contenu de results/week3/runs.csv (une ligne par epoch).
    generation : contenu de results/week3/generation_summary.csv.
    baselines : perplexites de reference, ex. {'unigramme': 19.4874, ...}.
    manifest_csv : chemin du manifeste SHA-256 des checkpoints.
    generated_at : horodatage a afficher (defaut : maintenant, UTC).
    """
    generated_at = generated_at or utc_timestamp()

    best = best_run_table(runs)
    summary = model_summary(runs)

    lines = [
        f"Genere le {generated_at} par `notebooks/NB03_colab_train.ipynb`.",
        "",
        "### Environnement et code",
        "",
    ]

    if not runs.empty:
        meta = runs.iloc[0]
        env_rows = pd.DataFrame(
            [
                {"champ": "commit git", "valeur": meta.get("git_sha", "unknown")},
                {
                    "champ": "arbre git modifie",
                    "valeur": str(bool(meta.get("git_dirty", True))),
                },
                {"champ": "device", "valeur": meta.get("device", "")},
                {"champ": "GPU", "valeur": meta.get("gpu_name", "") or "aucun"},
                {"champ": "PyTorch", "valeur": meta.get("torch_version", "")},
                {"champ": "Python", "valeur": meta.get("python_version", "")},
            ]
        )
        lines.append(_format_table(env_rows, ["champ", "valeur"]))
    else:
        lines.append("_runs.csv absent : aucun run a documenter._")

    lines += [
        "",
        "### Commandes",
        "",
        "```bash",
        "python -m peplib.train --config configs/lstm_base.yaml --seed 42",
        "python -m peplib.train --config configs/transformer_base.yaml --seed 42",
        "python scripts/generate_week3.py --config configs/lstm_base.yaml "
        "--checkpoint checkpoints/lstm/seed_42/best.pt --seed 42 --n 1000",
        "python scripts/fix_visuals.py",
        "```",
        "",
        "### Resultats d'entrainement (meilleure epoch par run)",
        "",
        _format_table(
            best,
            ["model", "seed", "epoch", "train_loss", "val_loss", "val_pp"],
        ),
        "",
        "### Agregation par modele",
        "",
        _format_table(summary, ["model", "n_seeds", "val_pp_mean", "val_pp_std"]),
        "",
    ]

    # --- Reference de perplexite et gate ---
    lines.append("### Reference de perplexite (validation v0.2)")
    lines.append("")
    if baselines:
        base_rows = pd.DataFrame(
            [{"baseline": name, "perplexite": round(float(value), 4)}
             for name, value in sorted(baselines.items())]
        )
        lines.append(_format_table(base_rows, ["baseline", "perplexite"]))
        gate = baselines.get("unigramme_positionnel")
        if gate is not None and not summary.empty:
            lines.append("")
            lines.append(
                f"Gate retenu : `val_pp < {gate:.4f}` (unigramme positionnel)."
            )
            for _, row in summary.iterrows():
                verdict = "OK" if row["val_pp_mean"] < gate else "ECHEC"
                lines.append(
                    f"- {row['model']} : {row['val_pp_mean']:.4f} -> {verdict}"
                )
    else:
        lines.append("_baselines non calculees._")

    # --- Generation ---
    lines += ["", "### Generation sur checkpoints entraines", ""]
    if generation is not None and not generation.empty:
        cols = [
            c
            for c in [
                "model",
                "seed",
                "n_attempts",
                "n_valid",
                "taux_validite_brut",
                "taux_unicite",
                "taux_copie_exacte",
                "fraction_proche_0.7",
                "fraction_proche_0.9",
            ]
            if c in generation.columns
        ]
        lines.append(_format_table(generation, cols))
    else:
        lines.append("_generation non executee._")

    lines += [
        "",
        "### Checkpoints",
        "",
        f"Hashes SHA-256 dans `{Path(manifest_csv).as_posix()}`.",
        "Les checkpoints eux-memes ne sont pas versionnes (volume) : ils sont "
        "conserves sur Google Drive et leur hash permet de verifier qu'un "
        "fichier retrouve est bien celui du run publie.",
        "",
        "### Limites declarees",
        "",
        "- Les scorers et metriques de generation sont des proxys : aucune "
        "activite biologique n'est mesuree.",
        "- La proximite au train est estimee par Smith-Waterman apres "
        "prefiltrage par 3-mers ; les valeurs sous le seuil peuvent etre "
        "sous-estimees (documente dans `peplib/metrics_proximity.py`).",
        "- Aucun predicteur de structure ni poids pre-entraine n'est utilise.",
    ]
    return "\n".join(lines)
