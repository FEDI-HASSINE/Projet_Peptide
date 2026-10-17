import hashlib
import os
import re
from datetime import datetime
import requests

# 1. Créer le dossier data/raw s'il n'existe pas encore
RAW_DIR = os.path.join("data", "raw")
os.makedirs(RAW_DIR, exist_ok=True)


def compute_sha256(filepath: str) -> str:
    """Calcule l'empreinte SHA-256 du fichier brut (pour T10)."""
    sha256_hash = hashlib.sha256()
    with open(filepath, "rb") as f:
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)
    return sha256_hash.hexdigest()


def download_uniprot():
    date_str = datetime.now().strftime("%Y%m%d")
    filepath = os.path.join(RAW_DIR, f"uniprot_reviewed_{date_str}.tsv")

    # Requête demandée dans T7 : longueur 10 à 50, reviewed:true
    query = "length:[10 TO 50] AND reviewed:true"
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": query,
        "format": "tsv",
        "fields": "accession,sequence,length,protein_name,keyword",
        "size": 500,
    }

    session = requests.Session()
    url = base_url
    first_page = True
    total_sequences = 0

    print("Début du téléchargement depuis UniProt...")

    with open(filepath, "w", encoding="utf-8", newline="\n") as f_out:
        while url:
            response = session.get(url, params=params if first_page else None)
            response.raise_for_status()

            lines = response.text.splitlines()
            if not lines:
                break

            if first_page:
                f_out.write("\n".join(lines) + "\n")
                total_sequences += len(lines) - 1
                first_page = False
            else:
                if len(lines) > 1:
                    f_out.write("\n".join(lines[1:]) + "\n")
                    total_sequences += len(lines) - 1

            print(f"  -> {total_sequences} séquences téléchargées...")

            # Gestion de la pagination (page suivante via l'en-tête Link)
            link_header = response.headers.get("Link")
            if link_header and 'rel="next"' in link_header:
                match = re.search(r'<([^>]+)>;\s*rel="next"', link_header)
                url = match.group(1) if match else None
            else:
                url = None

    # Calcul du nombre total de lignes brutes et du SHA-256 pour PROVENANCE.md (T10)
    with open(filepath, "r", encoding="utf-8") as f:
        raw_lines = sum(1 for _ in f)
    sha256 = compute_sha256(filepath)

    print("\n" + "=" * 50)
    print("✅ TÉLÉCHARGEMENT UNIPROT TERMINÉ (T7) !")
    print(f"📁 Fichier sauvegardé : {filepath}")
    print(
        f"📊 Nombre de lignes brutes : {raw_lines} (soit {total_sequences} séquences + 1 ligne d'en-tête)"
    )
    print(f"🔒 Empreinte SHA-256 : {sha256}")
    print("=" * 50)


if __name__ == "__main__":
    download_uniprot()