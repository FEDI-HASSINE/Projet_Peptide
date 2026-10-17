import os
import pandas as pd

RAW_FILE = "data/raw/uniprot_reviewed_20261017.tsv"
os.makedirs("data/processed", exist_ok=True)
OUT_FLAGS = "data/processed/toxin_flags.json"

# Charger le TSV UniProt en remplaçant les vides par des chaînes vides
df = pd.read_csv(RAW_FILE, sep="\t").fillna("")

# Identifier dynamiquement la colonne d'identifiant (ex. 'Entry', 'Accession' ou première colonne)
id_col = (
    "Entry"
    if "Entry" in df.columns
    else ("Accession" if "Accession" in df.columns else df.columns[0])
)

# Concaténer toutes les colonnes textuelles de manière sécurisée en convertissant en string
text_data = df.astype(str).agg(lambda x: " ".join(x), axis=1)

# Détecter les toxines (mots-clés 'Toxin' ou 'Conotoxin')
df["is_toxin"] = text_data.str.contains("Toxin|Conotoxin", case=False, na=False)

toxin_count = int(df["is_toxin"].sum())
print(f"Nombre de toxines détectées et flaggées : {toxin_count}")

# Exporter le rapport de flags au format JSON dans data/processed/
df[[id_col, "is_toxin"]].to_json(OUT_FLAGS, orient="records", indent=2)
print(f"Rapport de biosécurité enregistré sous {OUT_FLAGS}")