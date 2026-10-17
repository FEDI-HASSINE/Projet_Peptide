# Provenance des données (Corpus Peptides — Semaine 1)

Ce document retrace l'origine, les requêtes, les licences et les empreintes cryptographiques (SHA-256) des fichiers bruts stockés dans `data/raw/`. Conformément au cahier des charges, les fichiers dans `data/raw/` ne sont jamais modifiés.

---

## Source 1 : UniProtKB (Corpus général — T7)

- **Nom de la source** : UniProtKB (`reviewed:true` / Swiss-Prot)
- **URL de requête** : `https://rest.uniprot.org/uniprotkb/search?query=length:[10 TO 50] AND reviewed:true&format=tsv&fields=accession,sequence,length,protein_name,keyword&size=500`
- **Date de téléchargement** : 17/10/2026
- **Version de la base** : UniProt Release 2026 (API REST)
- **Conditions d'usage / licence** : Creative Commons Attribution (CC BY 4.0) — UniProt Terms of Use (vérifié le 17/10/2026, redistribution et usage académique autorisés avec citation).
- **Nombre de lignes brutes** : 12 424 (12 423 séquences + 1 ligne d'en-tête)
- **Empreinte SHA-256 du fichier brut** : `397959390c4ca2e6f55c5356e281f849398de0461328082f7549e2e1ad42f285`
- **Fichier brut sauvegardé sous** : `data/raw/uniprot_reviewed_20261017.tsv`
- **Téléchargement manuel ?** : Non (script `download_uniprot.py` avec pagination via en-tête `Link`).
- **Notes (Biosécurité)** : Le champ `keyword` et `protein_name` sont récupérés afin d'identifier et de signaler/exclure les toxines (ex. conotoxines, `KW-0800`) lors du prétraitement génératif.

---

## Source 2 : Base de Peptides Antimicrobiens — Profil Pharmaceutique (T8)

- **Nom de la source** : Référentiel AMP curé (Peptides antimicrobiens / Antibiotiques `KW-0929` & `KW-0044`, équivalent DRAMP/APD3)
- **URL de requête** : `https://rest.uniprot.org/uniprotkb/search?query=(keyword:KW-0929 OR keyword:KW-0044) AND reviewed:true&format=tsv&fields=accession,sequence,length,protein_name,organism_name&size=500`
- **Date de téléchargement** : 17/10/2026
- **Version de la base** : Release 2026
- **Conditions d'usage / licence** : Creative Commons Attribution (CC BY 4.0) — Usage académique et redistribution autorisés.
- **Nombre de lignes brutes** : 3 609 (3 608 séquences + 1 ligne d'en-tête)
- **Empreinte SHA-256 du fichier brut** : `ef21ecd472244cbd68d101651405ea90efd8f8e445fbe08b0769530dffe00134`
- **Fichier brut sauvegardé sous** : `data/raw/amp_database_20261017.tsv`
- **Téléchargement manuel ?** : Non (script `download_t8_t9.py` avec pagination).
- **Notes** : Contient des séquences de toutes longueurs ; le filtrage sur la fenêtre de longueur [10, 50] et sur l'alphabet canonique des 20 acides aminés est appliqué à l'étape C (T12–T14).

---

## Source 3 : Bactériocines de bioconservation — Profil Non-Pharmaceutique (T9)

- **Nom de la source** : Référentiel Bactériocines (type BACTIBASE / `KW-0078` Bacteriocin — bioconservation alimentaire)
- **URL de requête** : `https://rest.uniprot.org/uniprotkb/search?query=keyword:KW-0078 AND length:[5 TO 100]&format=tsv&fields=accession,sequence,length,protein_name,organism_name&size=500`
- **Date de téléchargement** : 17/10/2026
- **Version de la base** : Release 2026
- **Conditions d'usage / licence** : Creative Commons Attribution (CC BY 4.0) — Usage académique et redistribution autorisés.
- **Nombre de lignes brutes** : 787 (786 séquences + 1 ligne d'en-tête)
- **Empreinte SHA-256 du fichier brut** : `3832363fd7db973a4742a8dbea1335fccac19f65de4bccaed4f39d24e1ffceb4`
- **Fichier brut sauvegardé sous** : `data/raw/bactibase_bacteriocins_20261017.tsv`
- **Téléchargement manuel ?** : Non (script `download_t8_t9.py` avec pagination).
- **Notes** : Cible les peptides de bioconservation alimentaire (bactériocines). Le filtrage strict sur la longueur 10–50 et les caractères non standards est effectué dans `peplib/filters.py`.