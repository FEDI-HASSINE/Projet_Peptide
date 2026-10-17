import hashlib
import os
import re
from datetime import datetime
import requests
import urllib3

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

RAW_DIR = os.path.join("data", "raw")
os.makedirs(RAW_DIR, exist_ok=True)


def compute_sha256(filepath: str) -> str:
    sha256_hash = hashlib.sha256()
    with open(filepath, "rb") as f:
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)
    return sha256_hash.hexdigest()


def count_lines(filepath: str) -> int:
    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        return sum(1 for _ in f)


def download_paginated_uniprot(query: str, filepath: str) -> int:
    """Télécharge toutes les pages d'une requête REST."""
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": query,
        "format": "tsv",
        "fields": "accession,sequence,length,protein_name,organism_name",
        "size": 500,
    }
    session = requests.Session()
    next_url = base_url
    first_page = True
    total_seq = 0

    with open(filepath, "w", encoding="utf-8", newline="\n") as f_out:
        while next_url:
            r = session.get(
                next_url, params=params if first_page else None, timeout=30
            )
            r.raise_for_status()
            lines = r.text.splitlines()
            if not lines:
                break
            if first_page:
                f_out.write("\n".join(lines) + "\n")
                total_seq += len(lines) - 1
                first_page = False
            else:
                if len(lines) > 1:
                    f_out.write("\n".join(lines[1:]) + "\n")
                    total_seq += len(lines) - 1

            link = r.headers.get("Link")
            if link and 'rel="next"' in link:
                m = re.search(r'<([^>]+)>;\s*rel="next"', link)
                next_url = m.group(1) if m else None
            else:
                next_url = None
    return total_seq


def download_t8_amp():
    """T8 — Télécharge la base de peptides antimicrobiens (AMP)."""
    date_str = datetime.now().strftime("%Y%m%d")
    filepath = os.path.join(RAW_DIR, f"amp_database_{date_str}.tsv")

    print("1. Téléchargement de T8 (Base AMP complète)...")
    # On tente d'abord le fichier officiel DRAMP s'il est joignable
    dramp_urls = [
        "http://dramp.cpu-bioinfor.org/downloads/download.php?filename=download_data/DRAMP_General_amps.txt",
        "https://raw.githubusercontent.com/CPU-DRAMP/DRAMP_2.0/master/general_amps.txt",
    ]
    for url in dramp_urls:
        try:
            r = requests.get(url, timeout=15, verify=False)
            if (
                r.status_code == 200
                and len(r.text) > 5000
                and "<html" not in r.text[:300].lower()
            ):
                with open(filepath, "w", encoding="utf-8", newline="\n") as f:
                    f.write(r.text)
                lines = count_lines(filepath)
                sha = compute_sha256(filepath)
                print(f"✅ T8 (DRAMP) terminé ! Fichier : {filepath}")
                print(f"   URL utilisée : {url}")
                print(f"   Lignes brutes : {lines}")
                print(f"   SHA-256 : {sha}\n")
                return
        except Exception:
            pass

    # Si DRAMP est bloqué par le pare-feu, on récupère l'intégralité des AMP avec pagination complète
    query = "(keyword:KW-0929 OR keyword:KW-0044) AND reviewed:true"
    total_seq = download_paginated_uniprot(query, filepath)
    lines = count_lines(filepath)
    sha = compute_sha256(filepath)
    print(f"✅ T8 terminé ! Fichier : {filepath}")
    print(f"   Requête : {query}")
    print(f"   Lignes brutes : {lines} (soit {total_seq} séquences)")
    print(f"   SHA-256 : {sha}\n")


def download_t9_bacteriocins():
    """T9 — Télécharge la source non-pharmaceutique (Bactériocines de bioconservation alimentaire)."""
    date_str = datetime.now().strftime("%Y%m%d")
    filepath = os.path.join(RAW_DIR, f"bactibase_bacteriocins_{date_str}.tsv")

    print(
        "2. Téléchargement de T9 (Non-Pharma : Bactériocines de bioconservation)..."
    )
    # KW-0078 = Bacteriocin (syntaxe officielle valide sur UniProt REST)
    query = "keyword:KW-0078 AND length:[5 TO 100]"
    total_seq = download_paginated_uniprot(query, filepath)

    raw_lines = count_lines(filepath)
    sha = compute_sha256(filepath)
    print(f"✅ T9 terminé ! Fichier : {filepath}")
    print(f"   Requête : {query}")
    print(f"   Lignes brutes : {raw_lines} (soit {total_seq} séquences)")
    print(f"   SHA-256 : {sha}\n")


if __name__ == "__main__":
    download_t8_amp()
    download_t9_bacteriocins()