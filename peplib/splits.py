"""Splits 90/5/5 par familles (T39 / ajustement A6) et split aléatoire de contrôle.

Conventions reprises de A1 (proposition, à confirmer mardi) :
  - découpe 90/5/5 des SÉQUENCES, tolerance ±0,5 point par partie ;
  - familles (clusters) entères : aucune famille à cheval train/val/test ;
  - les familles de plus de 1 % du corpus sont systématiquement en train
    (sinon un unique cluster déborderait de sa cible) ;
  - stratifié par source d'ingestion (`UniProt` / `APD3` / `AMPChallenge`),
    atomes de la colonne `sources` ;
  - seed 42 partout.

Algorithme d'affectation (déterministe, documenté, à valider en PR) :
  1. familles de taille > 1 % du corpus -> train ;
  2. familles restantes triées par taille décroissante (égalité : numéro
     de famille croissant) ;
  3. chaque famille va dans la partie dont le RESTE À POURVOIR est le plus
     grand (capacité restante en séquences), à égalité : train, puis val ;
  4. contrôle final : chaque partie dans sa tolerance ±0,5 pt ; sinon,
     recherche locale : on tente de déplacer la plus grande famille
     déplaçable de la partie en excès vers la partie en déficit
     (1000 itérations max) ; si l'échec persiste, ValueError — et non
     un epsilon de tolerance improvisé.

Le split aléatoire de contrôle (T16) utilise les MÊMES tailles exactes de
parties, seedée 42 ; `cluster_id` y est conservé (par construction il
traverse les parties, ce qui est le point du contrôle).

Clause conservée (PR #19 W4) : `toxicity_status` accompagne chaque séquence ;
`unknown` est traité comme `flagged` par tout consommateur génératif — aucune
exclusion décidée ici (arbitrage d'équipe, §6).
"""

import hashlib
import json
from pathlib import Path

import pandas as pd

SEED = 42
TAUX = (0.90, 0.05, 0.05)
TOLERANCE = 0.005            # ±0,5 point par partie
SEUIL_GRANDE = 0.01          # famille de > 1 % du corpus -> train
PARTIES = ("train", "val", "test")
SPLIT_DIR = "data/splits"
MANIFEST_PATH = "splits_manifest.json"


def assigner_familles(
    comptes_familles: dict,
    sources_atomes: dict | None = None,
    total: int | None = None,
    taux: tuple = TAUX,
    tolerance: float = TOLERANCE,
    seuil_grande: float = SEUIL_GRANDE,
) -> dict:
    """Affecte chaque famille à une partie, de façon déterministe.

    comptes_familles : {cluster_id: effectif}. sources_atomes : optionnel,
    {cluster_id: {atome: effectif}} — atomes de la colonne `sources`. Retour
    {cluster_id: "train"|"val"|"test"}.

    Algorithme (documenté, à valider en PR) :

      1. familles de plus de ``seuil_grande`` × total → train sans exception ;
      2. familles restantes triées par taille décroissante (égalité :
         numéro de famille croissant) ;
      3. affectation gloutonne au moindre coût : pour chaque famille, le
         coût d'une partie vaut la somme des carrés des dépassements —
         dépassement de taille (effectif partie > cible taille) et
         dépassements par atome (effectif atome de la partie > cible
         atome = taux_part × effectif global de l'atome) ; égalité de
         coût : ordre train < val < test ;
      4. vérification finale ±tolérance par partie ; sinon ValueError avec
         le relevé des écarts (aucun epsilon improvisé).

    Ce coût par atome EST la stratification par source : chaque partie est
    poussée vers la mixture globale (UniProt / APD3 / AMPChallenge) selon
    les familles disponibles. Les mixtures obtenues sont mesurées et
    écrites dans le manifeste.
    """
    if total is None:
        total = sum(comptes_familles.values())
    cibles = {part: taux[i] * total for i, part in enumerate(PARTIES)}

    atomes = sorted(
        {a for vect in (sources_atomes or {}).values() for a in vect}
    )
    totaux_atomes = {
        atome: sum(
            vect.get(atome, 0) for vect in (sources_atomes or {}).values()
        )
        for atome in atomes
    }
    cibles_atomes = {
        part: {
            atome: taux[PARTIES.index(part)] * totaux_atomes[atome]
            for atome in atomes
        }
        for part in PARTIES
    }

    affectation = {}
    remplies = {part: 0 for part in PARTIES}
    rempli_atomes = {part: {a: 0 for a in atomes} for part in PARTIES}

    # 1. familles trop grandes : train, sans exception.
    for cid, n in comptes_familles.items():
        if n > seuil_grande * total:
            affectation[cid] = "train"
            remplies["train"] += n
            for atome, v in (sources_atomes or {}).get(cid, {}).items():
                rempli_atomes["train"][atome] += v

    restantes = sorted(
        (cid for cid in comptes_familles if cid not in affectation),
        key=lambda cid: (-comptes_familles[cid], cid),
    )

    # 2-3. affectation gloutonne au moindre coût (taille + atomes).
    for cid in restantes:
        n = comptes_familles[cid]
        vect = (sources_atomes or {}).get(cid, {})
        couts = []
        for indice_part, part in enumerate(PARTIES):
            depasse = max(0, remplies[part] + n - cibles[part])
            cout = depasse * depasse
            for atome, valeur_atomique in vect.items():
                depasse_atom = max(
                    0,
                    rempli_atomes[part][atome]
                    + valeur_atomique
                    - cibles_atomes[part].get(atome, 0),
                )
                cout += depasse_atom * depasse_atom
            couts.append((cout, indice_part, part))
        affectation[cid] = min(couts)[2]
        part = affectation[cid]
        remplies[part] += n
        for atome, v in vect.items():
            rempli_atomes[part][atome] += v

    # 4. vérification finale.
    hors = {
        part: int(round(remplies[part] - cibles[part]))
        for part in PARTIES
        if abs(remplies[part] - cibles[part]) > tolerance * total
    }
    if hors:
        raise ValueError(f"tolérance ±0,5 pt violée : écarts {hors}")
    return affectation


def ecrire_parts(cumul, out_dir=SPLIT_DIR) -> dict:
    """Écrit les trois CSV LF sans BOM, retourne les empreintes {part: sha256}."""
    empreintes = {}
    repertoire = Path(out_dir)
    repertoire.mkdir(parents=True, exist_ok=True)
    for part in PARTIES:
        part_path = repertoire / f"{part}_v02.csv"
        cumul[part].to_csv(part_path, index=False, lineterminator="\n")
        empreintes[part] = hashlib.sha256(part_path.read_bytes()).hexdigest()
    return empreintes


def split_aleatoire_control(rows: pd.DataFrame, cibles: dict, seed: int = SEED):
    """Split aléatoire de contrôle : mêmes tailles exactes, graine 42.

    cibles : {"train": n, "val": n, "test": n}. Retour {part: DataFrame}.
    """
    melange = rows.sample(frac=1, random_state=seed).reset_index(drop=True)
    n_train = cibles["train"]
    n_val = cibles["val"]
    return {
        "train": melange.iloc[:n_train],
        "val": melange.iloc[n_train : n_train + n_val],
        "test": melange.iloc[n_train + n_val :],
    }


def ecrire_manifeste(
    comptes_parties, comptes_source_part, empreintes, chemin=MANIFEST_PATH
):
    """Manifeste des splits + empreintes (test anti-drift à l'appui)."""
    manifeste = {
        "version": "v0.2",
        "seed": SEED,
        "methodee": "split 90/5/5 par familles (clusters v0.2 couverture 0,90), "
        "stratifié par source d'ingestion, familles > 1 % en train",
        "tolerance": "±0,5 point",
        "comptes_par_partie": comptes_parties,
        "comptes_partie_par_source": comptes_source_part,
        "sha256": empreintes,
    }
    Path(chemin).write_text(
        json.dumps(manifeste, sort_keys=True, indent=2) + "\n",
        encoding="utf-8",
    )
    return manifeste


SEPARATEUR_SOURCES = "+"


def run(
    clusters_path="data/processed/clusters_v02_90.csv",
    corpus_path="data/processed/corpus_v0.2.csv",
    toxin_path="data/processed/toxin_flags_v02.csv",
    out_dir=SPLIT_DIR,
    manifest=MANIFEST_PATH,
    control_dir="data/splits_control",
):
    """Construit les splits T39/A6 + le contrôle aléatoire + le manifeste.

    Retourne (manifeste, statistiques) ; écrit :
      data/splits/{train,val,test}_v02.csv
      data/splits_control/{train,val,test}_v02.csv
      splits_manifest.json
    """
    clusters = pd.read_csv(clusters_path, dtype=str)
    if not {"sequence", "cluster_id"}.issubset(clusters.columns):
        raise ValueError("Colonnes manquantes dans l'artefact de clusters")

    corpus = pd.read_csv(corpus_path, dtype=str)
    if not {"sequence", "sources", "ids"}.issubset(corpus.columns):
        raise ValueError("Colonnes manquantes dans le corpus v0.2")

    toxin = pd.read_csv(toxin_path, dtype=str)
    if not {"sequence", "toxicity_status"}.issubset(toxin.columns):
        raise ValueError("Colonnes manquantes dans toxin_flags_v02.csv")

    if corpus["sequence"].nunique() != len(corpus):
        raise ValueError("Le corpus contient des doublons exacts")
    clusters_map = dict(zip(clusters["sequence"], clusters["cluster_id"]))
    manquants = [s for s in corpus["sequence"] if s not in clusters_map]
    if manquants:
        raise ValueError(f"{len(manquants)} séquences sans cluster")

    toxin_map = dict(zip(toxin["sequence"], toxin["toxicity_status"]))
    sans_annotation = [
        s for s in corpus["sequence"] if s not in toxin_map
    ]

    rows = corpus.copy()
    rows["cluster_id"] = rows["sequence"].map(clusters_map)
    rows["toxicity_status"] = rows["sequence"].map(
        lambda s: toxin_map.get(s, "unknown")
    )
    rows = rows[
        ["sequence", "sources", "ids", "cluster_id", "toxicity_status"]
    ]

    # Affectation au niveau des familles, guidée par la mixture par source.
    comptes_familles = rows.groupby("cluster_id").size().to_dict()
    total = len(rows)
    sources_atomes = {}
    for cid, groupe in rows.groupby("cluster_id"):
        vect: dict = {}
        for chaine in groupe["sources"]:
            for atome in chaine.split(SEPARATEUR_SOURCES):
                vect[atome] = vect.get(atome, 0) + 1
        sources_atomes[cid] = vect
    affectation = assigner_familles(
        comptes_familles, sources_atomes=sources_atomes, total=total
    )

    parts = {
        part: rows[rows["cluster_id"].map(affectation) == part].reset_index(
            drop=True
        )
        for part in PARTIES
    }

    # Rapport de stratification par source (atomes décomposés).
    comptes_source_part = {}
    for part in PARTIES:
        atomes = {}
        for chaine in parts[part]["sources"]:
            for atome in chaine.split(SEPARATEUR_SOURCES):
                atomes[atome] = atomes.get(atome, 0) + 1
        comptes_source_part[part] = atomes

    comptes_parties = {part: len(parts[part]) for part in PARTIES}
    empreintes = ecrire_parts(parts, out_dir=out_dir)

    # Split aléatoire de contrôle : mêmes tailles exactes, graine 42.
    controle = split_aleatoire_control(
        rows, cibles=comptes_parties
    )
    empreintes_control = ecrire_parts(controle, out_dir=control_dir)
    controle_dir = Path(control_dir)

    manifeste = ecrire_manifeste(
        comptes_parties,
        comptes_source_part,
        empreintes,
        chemin=manifest,
    )
    manifeste["controle_aleatoire"] = {
        "repertoire": controle_dir.as_posix(),
        "methode": "tirage aleatoire des memes tailles, seed 42",
        "sha256": empreintes_control,
    }
    Path(manifest).write_text(
        json.dumps(manifeste, sort_keys=True, indent=2) + "\n",
        encoding="utf-8",
    )

    stats = {
        "total": total,
        "familles": len(comptes_familles),
        "famille_plus_grande": max(comptes_familles.values()),
        "sans_annotation_toxin": len(sans_annotation),
        "comptes_parties": comptes_parties,
    }
    if sans_annotation:
        print(
            f"ATTENTION : {len(sans_annotation)} séquences absentes de "
            "toxin_flags_v02.csv -> toxicity_status=\"unknown\" "
            "(traitées comme flagged en usage génératif, clause W4)."
        )
    print(json.dumps(stats, indent=2))
    return manifeste, stats
