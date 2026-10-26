"""Clusters du corpus v0.2 au critère de couverture (A5, AJUSTEMENTS §3/§5).

Souvent résumé ainsi : deux séquences sont liées quand il existe un alignement
local (SW, o = 11, e = 1, BLOSUM62) dont l'identité ET la couverture valent au
moins le seuil tau (0,90 ; variante 0,70 rapportée). Les familles sont les
composantes connexes (union-find, chaînage).

Première version conservée : le clustering Levenshtein v0.1 (PR #16) reste en
place comme mesure complémentaire (T14(b) de l'encadrant) et n'est pas
remplacé par ce module.

Préfiltrage (borne démontrée : ``Docs/A6_preuve_borne_couverture.md``) :
    U(a, b) = Σ_s min(c_a(s), c_b(s))  sur les 3-mers partagés (multiset)
est une borne supérieure des occurrences intactes alignées ; si
U(a, b) < ⌊τ·m·(1 − k(1−τ)) − k⌋ avec m = min(len(a), len(b)), la paire est
rejetée sans Smith-Waterman (rejet sonore, aucun faux négatif).

Comparaison au seuil : les ratios produits ont des dénominateurs W ≤ 200
(2 séquences de 50 avec gaps) : deux valeurs de ratio distinctes sont
distantes d'au moins 1/200 = 0,005 ≫ 1e-16, l'écart flottant ne peut donc pas
faire franchir le seuil 0,90 (ni 0,70) par erreur d'arrondi — les
comparaisons ``>=`` sur floats sont exactes ici, sans epsilon arbitraire.

Sortie officielle : ``data/processed/clusters_v02_90.csv``
(``sequence, sequence_id, cluster_id, sources``) — même contrat que le
verrou v0.1 (``clusters_v01_90.csv``), la colonne monolithique ``source`` et
``id`` de la Semaine 1 devenant ``sources``/``ids`` dans le corpus v0.2.
"""

import math
from pathlib import Path

import pandas as pd

from peplib.cluster import UnionFind
from peplib.identite import identite_couverture

TAU = 0.90
K_NOMMER = 3          # k des k-mers (convention DECISIONS : k = 3)
OUT_PATH = "data/processed/clusters_v02_90.csv"
T13_PATH = "results/week2/clusters_v02_t13.csv"


def borne_intacts(m: int, tau: float = TAU, k: int = K_NOMMER) -> int:
    """K(m) = ⌊τ·m·(1 − k(1−τ)) − k⌋, majoration sans epsilon (preuve A6 §2)."""
    if not 0 < tau <= 1:
        raise ValueError("tau doit être dans ]0, 1]")
    return max(0, math.floor(tau * m * (1 - k * (1 - tau)) - k))


def candidats_par_couverture(
    sequences: list, tau: float = TAU, k: int = K_NOMMER, stats: dict | None = None
) -> set:
    """Paires survivant à la borne sonore sur les occurrences alignées.

    Les occurrences (multiset) et non les suffisent types — cf. preuve A6 §4.
    Retourne un ensemble de couples (i, j) avec i < j (sous-ensemble qui
    contient toutes les paires susceptibles de satisfaire le critère).
    """
    occurrences = [None] * len(sequences)
    groupes = {}
    for idx, s in enumerate(sequences):
        compte = {}
        for pos in range(len(s) - k + 1):
            km = s[pos : pos + k]
            compte[km] = compte.get(km, 0) + 1
            membres = groupes.get(km)
            if membres is None:
                groupes[km] = [idx]
            elif membres[-1] != idx:      # présence par type : un seul idx par séquence
                membres.append(idx)
        occurrences[idx] = compte

    seuils = {}   # cache des seuils K(m) par longueur minimale
    partiel = {}
    retenues = set()

    for km, membres in groupes.items():
        n = len(membres)
        for a in range(n - 1):
            i = membres[a]
            ci = occurrences[i][km]
            len_i = len(sequences[i])
            for b in range(a + 1, n):
                j = membres[b]
                paire = (i, j)
                if paire in retenues:
                    continue
                u = partiel.get(paire, 0) + min(ci, occurrences[j][km])
                m = len_i
                lj = len(sequences[j])
                if lj < m:
                    m = lj
                seuil = seuils.get(m)
                if seuil is None:
                    seuil = borne_intacts(m, tau, k)
                    seuils[m] = seuil
                if u >= seuil:
                    retenues.add(paire)
                    partiel.pop(paire, None)
                else:
                    partiel[paire] = u

    if stats is not None:
        stats["couples_partages"] = sum(
            n * (n - 1) // 2 for n in map(len, groupes.values())
        )
        stats["candidats_couverture"] = len(retenues)

    return retenues


def cluster_v02(
    sequences: list,
    tau: float = TAU,
    gap_open: int = 11,
    gap_extend: int = 1,
    stats: dict | None = None,
    candidats: set | None = None,
) -> list:
    """Clustering union-find : paires liées si identité ∧ couverture ≥ tau.

    Le retour est une liste d'identifiants de cluster alignés sur l'ordre
    d'entrée de `sequences` (numérotation déterministe des racines triées).
    `candidats` permet de réutiliser un ensemble de paires déjà produit par
    ``candidats_par_couverture`` (variantes de gaps : la borne A6 ne dépend
    que de tau, pas de (o,e) — le même ensemble est sain pour toutes les
    conventions de gap de la famille BLASTP testées).
    """
    if not 0 < tau <= 1:
        raise ValueError("tau doit être dans ]0, 1]")

    unions = UnionFind(len(sequences))
    if candidats is None:
        candidats = candidats_par_couverture(sequences, tau=tau, stats=stats)

    lignees = 0
    comparees = 0
    for i, j in sorted(candidats):
        a = sequences[i]
        b = sequences[j]
        identite, couverture = identite_couverture(a, b, gap_open, gap_extend)
        comparees += 1
        if identite >= tau and couverture >= tau:
            if unions.union(i, j):
                lignees += 1
        if comparees % 20000 == 0:
            print(f"Comparaisons : {comparees} | Fusions : {lignees}", flush=True)

    racines = [unions.find(i) for i in range(len(sequences))]
    numerotation = {
        racine: num for num, racine in enumerate(sorted(set(racines)))
    }
    etiquettes = [numerotation[racine] for racine in racines]

    if stats is not None:
        stats["comparees_sw"] = comparees
        stats["fusions"] = lignees

    return etiquettes


def tableau_t13(id_liste: list[int], sortie=T13_PATH) -> pd.DataFrame:
    """T13 : nb familles, taille max, singletons, fraction du plus grand."""
    lignes = []
    tailles = {}
    for cid in id_liste:
        tailles[cid] = tailles.get(cid, 0) + 1
    grande = max(tailles.values())
    lignes.append(
        {
            "nb_familles": len(tailles),
            "taille_max": grande,
            "singletons": sum(1 for t in tailles.values() if t == 1),
            "fraction_plus_grande": round(grande / len(id_liste), 6),
        }
    )
    df = pd.DataFrame(lignes)
    if sortie is not None:
        Path(sortie).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(sortie, index=False, lineterminator="\n")
    return df


def run(
    input_path="data/processed/corpus_v0.2.csv",
    output_path=OUT_PATH,
    figure_path="results/week2/cluster_sizes_v02.png",
    tau: float = TAU,
):
    """Clustering v0.2 de bout en bout : artefact + tableau T13 + figure."""
    df = pd.read_csv(input_path, dtype=str)
    attendu = {"sequence", "sources", "ids"}
    if not attendu.issubset(df.columns):
        raise ValueError("Colonnes manquantes dans le corpus v0.2")
    if df["sequence"].isna().any().any():
        raise ValueError("Séquences manquantes dans le corpus v0.2")

    sequences = sorted(df["sequence"].unique().tolist())
    print(f"Séquences uniques : {len(sequences)}")

    stats = {}
    etiquettes = cluster_v02(sequences, tau=tau, stats=stats)
    print(stats)

    mapping = dict(zip(sequences, etiquettes))
    meta = df.drop_duplicates("sequence", keep="first")

    resultat = pd.DataFrame(
        {
            "sequence": meta["sequence"].to_numpy(),
            "sequence_id": meta["ids"].to_numpy(),
            "cluster_id": meta["sequence"].map(mapping).to_numpy(),
            "sources": meta["sources"].to_numpy(),
        }
    )
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    resultat.to_csv(output_path, index=False, lineterminator="\n")

    tableau_t13(etiquettes)

    import matplotlib.pyplot as plt

    tailles = pd.Series(etiquettes).value_counts()
    Path(figure_path).parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(9, 5))
    plt.hist(tailles, bins=30, edgecolor="black")
    plt.xlabel("Nombre de séquences uniques par famille")
    plt.ylabel("Nombre de familles")
    plt.title(f"Distribution des clusters v0.2 à {tau:.0%} (couverture)")
    plt.tight_layout()
    plt.savefig(figure_path, dpi=150)
    plt.close()

    print(f"CSV : {output_path} | Figure : {figure_path}")
    return resultat


def representants(
    clusters_path="data/processed/clusters_v02_90.csv",
    output_path="data/processed/representants_v02.csv",
    separateur="+",
):
    """T42 : un représentant par famille + mapping famille -> sources.

    Critère du représentant (écrit, pas improvisé) : la longueur maximale ;
    à égalité, l'ordre alphabétique croissant de la séquence. Les sources
    fusionnent les ATOMES de tous les membres de la famille (dédoublonnés,
    séparés par `separateur`, ordre canonique UniProt < APD3 < AMPChallenge).

    Colonnes de sortie : `sequence, sequence_id, cluster_id, sources,
    taille_famille` (une ligne par famille) — consommables pour éviter de
    générer deux peptides de la même famille.
    """
    clusters = pd.read_csv(clusters_path, dtype=str)
    if not {"sequence", "sequence_id", "cluster_id", "sources"}.issubset(
        clusters.columns
    ):
        raise ValueError("Colonnes manquantes dans l'artefact de clusters")

    ordre_atomes = {"UniProt": 0, "APD3": 1, "AMPChallenge": 2}

    lignes = []
    for cid, groupe in clusters.groupby("cluster_id"):
        atomes = set()
        for chaine in groupe["sources"]:
            atomes.update(chaine.split(separateur))
        atomes_tries = sorted(
            atomes, key=lambda a: ordre_atomes.get(a, 99)
        )
        meilleur = sorted(
            groupe.itertuples(index=False),
            key=lambda t: (-len(t.sequence), t.sequence),
        )[0]
        lignes.append(
            {
                "sequence": meilleur.sequence,
                "sequence_id": meilleur.sequence_id,
                "cluster_id": cid,
                "sources": separateur.join(atomes_tries),
                "taille_famille": len(groupe),
            }
        )

    resultat = pd.DataFrame(lignes).sort_values("cluster_id").reset_index(drop=True)
    resultat.to_csv(output_path, index=False, lineterminator="\n")
    print(f"Représentants : {len(resultat)} familles -> {output_path}")
    return resultat
