"""
Acquisition des corpus bruts depuis l'API REST UniProt (T7, T8, T9).

Toute la logique de téléchargement, de pagination et d'empreinte vit ici ;
les scripts à la racine (``download_uniprot.py``, ``download_t8_t9.py``) ne
sont que des points d'entrée.
"""

import hashlib
import os
import re
from datetime import datetime

import requests

UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"

# T7 — corpus général, peptides de 10 à 50 résidus.
QUERY_T7 = "length:[10 TO 50] AND reviewed:true"
FIELDS_T7 = "accession,sequence,length,protein_name,keyword"

# T8 — peptides antimicrobiens (profil pharmaceutique).
QUERY_T8 = "(keyword:KW-0929 OR keyword:KW-0044) AND reviewed:true"
FIELDS_T8 = "accession,sequence,length,protein_name,organism_name"

# T9 — bactériocines (profil non pharmaceutique).
QUERY_T9 = "keyword:KW-0078 AND length:[5 TO 100]"
FIELDS_T9 = "accession,sequence,length,protein_name,organism_name"

PAGE_SIZE = 500
DEFAULT_TIMEOUT = 30
RAW_DIR = os.path.join("data", "raw")


def compute_sha256(filepath: str) -> str:
    """Calcule l'empreinte SHA-256 du fichier brut (pour T10)."""
    sha256_hash = hashlib.sha256()
    with open(filepath, "rb") as f:
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)
    return sha256_hash.hexdigest()


def count_lines(filepath: str) -> int:
    """Compte le nombre de lignes d'un fichier."""
    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        return sum(1 for _ in f)


def next_page_url(link_header):
    """Extrait l'URL de la page suivante depuis l'en-tête ``Link`` d'UniProt.

    Retourne ``None`` s'il n'y a pas de page suivante.
    """
    if not link_header or 'rel="next"' not in link_header:
        return None
    match = re.search(r'<([^>]+)>;\s*rel="next"', link_header)
    return match.group(1) if match else None


def page_body_lines(page_lines, is_first):
    """Retourne les lignes à écrire pour une page.

    La première page contient l'en-tête TSV ; les suivantes le répètent et il
    faut donc le retirer.
    """
    if is_first:
        return page_lines
    return page_lines[1:] if len(page_lines) > 1 else []


def download_paginated_tsv(
    query, filepath, fields, timeout=DEFAULT_TIMEOUT, session=None
):
    """Télécharge toutes les pages d'une requête UniProt dans un TSV.

    Retourne le nombre de séquences téléchargées (hors ligne d'en-tête).
    ``session`` permet d'injecter un client factice dans les tests.
    """
    session = session or requests.Session()
    params = {
        "query": query,
        "format": "tsv",
        "fields": fields,
        "size": PAGE_SIZE,
    }
    next_url = UNIPROT_SEARCH_URL
    is_first = True
    total_sequences = 0

    directory = os.path.dirname(filepath)
    if directory:
        os.makedirs(directory, exist_ok=True)

    with open(filepath, "w", encoding="utf-8", newline="\n") as f_out:
        while next_url:
            response = session.get(
                next_url, params=params if is_first else None, timeout=timeout
            )
            response.raise_for_status()

            lines = response.text.splitlines()
            if not lines:
                break

            body = page_body_lines(lines, is_first)
            if body:
                f_out.write("\n".join(body) + "\n")
            total_sequences += len(lines) - 1
            is_first = False

            next_url = next_page_url(response.headers.get("Link"))

    return total_sequences


def raw_path(prefix: str, date_str=None) -> str:
    """Construit le chemin brut daté : ``data/raw/<prefix>_<YYYYMMDD>.tsv``."""
    date_str = date_str or datetime.now().strftime("%Y%m%d")
    return os.path.join(RAW_DIR, f"{prefix}_{date_str}.tsv")


def download_t7(filepath):
    """T7 — corpus général (peptides 10–50, ``reviewed:true``)."""
    return download_paginated_tsv(QUERY_T7, filepath, FIELDS_T7)


def download_t8(filepath):
    """T8 — peptides antimicrobiens."""
    return download_paginated_tsv(QUERY_T8, filepath, FIELDS_T8)


def download_t9(filepath):
    """T9 — bactériocines."""
    return download_paginated_tsv(QUERY_T9, filepath, FIELDS_T9)
