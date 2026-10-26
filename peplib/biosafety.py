"""
Détection de séquences toxiques à exclure du corpus génératif (biosécurité).

Le cahier des charges impose de signaler ou d'exclure les toxines (par exemple
les conotoxines présentes dans UniProt) et de ne jamais diffuser de séquences
générées de type toxine.
"""

import os

import pandas as pd

# Motif de détection appliqué aux annotations textuelles ("Protein names"
# et "Keywords"). "Toxin" couvre déjà "Conotoxin" en casse ignorée ; les
# deux termes sont conservés pour rester explicites. "Venom" couvre les
# protéines de venin (dont "envenomation") et "Toxic" couvre toxic /
# toxicity. Le motif est volontairement SUR-inclusif (une description
# "non-toxic" est aussi signalée) : le mode d'échec doit être sûr — le
# doute est toujours signalé (relecture PR #19, point B2).
TOXIN_PATTERN = "Toxin|Conotoxin|Venom|Toxic"

# Colonnes candidates pour l'identifiant, par ordre de préférence.
ID_COLUMNS = ("Entry", "Accession")


def load_raw_tsv(filepath):
    """Charge un TSV brut UniProt en remplaçant les valeurs manquantes par "".

    ``fillna`` intervient *avant* toute conversion en chaîne : sinon les valeurs
    manquantes deviendraient la chaîne "nan" et fausseraient la détection.
    """
    return pd.read_csv(filepath, sep="\t").fillna("")


def detect_id_column(df):
    """Retourne la colonne identifiant : ``Entry``, sinon ``Accession``, sinon la première."""
    for name in ID_COLUMNS:
        if name in df.columns:
            return name
    return df.columns[0]


def flag_toxins(df, pattern=TOXIN_PATTERN):
    """Retourne une copie de ``df`` avec une colonne booléenne ``is_toxin``.

    Toutes les colonnes sont concaténées en une chaîne unique par entrée, puis
    comparées au motif sans tenir compte de la casse.
    """
    text = df.fillna("").astype(str).agg(lambda row: " ".join(row), axis=1)
    flagged = df.copy()
    flagged["is_toxin"] = text.str.contains(pattern, case=False, na=False)
    return flagged


def write_toxin_flags(df, filepath, id_column=None):
    """Écrit le rapport de biosécurité ``(identifiant, is_toxin)`` au format JSON.

    ``df`` doit déjà contenir la colonne ``is_toxin`` (voir :func:`flag_toxins`).
    """
    if id_column is None:
        id_column = detect_id_column(df)
    directory = os.path.dirname(filepath)
    if directory:
        os.makedirs(directory, exist_ok=True)
    df[[id_column, "is_toxin"]].to_json(filepath, orient="records", indent=2)
    return filepath


def run(raw_file, out_file):
    """Pipeline complet : charge le TSV brut, flague les toxines, écrit le JSON.

    Retourne le nombre d'entrées signalées comme toxiques.
    """
    df = load_raw_tsv(raw_file)
    flagged = flag_toxins(df)
    toxin_count = int(flagged["is_toxin"].sum())
    write_toxin_flags(flagged, out_file)
    return toxin_count
