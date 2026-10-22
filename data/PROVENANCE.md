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
- **Téléchargement manuel ?** : Non (module `python -m peplib.download_uniprot` avec pagination via en-tête `Link`).
- **Notes (Biosécurité)** : Le champ `keyword` et `protein_name` sont récupérés afin d'identifier et de signaler/exclure les toxines (ex. conotoxines, `KW-0800`) lors du prétraitement génératif.

---

## Source 2 : Peptides antimicrobiens extraits d'UniProtKB — Profil Pharmaceutique (T8)

- **Nom de la source** : UniProtKB (API REST), filtré par les mots-clés `KW-0929` (Antimicrobial) et `KW-0044` (Antibiotic), `reviewed:true`
- **Base dédiée de référence (non utilisée)** : APD3 / DBAASP / DRAMP n'ont pas été téléchargées cette semaine. UniProt a été retenu pour disposer d'une API REST unique, reproductible et sous licence CC BY 4.0. **Écart assumé au cahier des charges T8**, à substituer par une base dédiée si le coordinateur l'exige.
- **URL de requête** : `https://rest.uniprot.org/uniprotkb/search?query=(keyword:KW-0929 OR keyword:KW-0044) AND reviewed:true&format=tsv&fields=accession,sequence,length,protein_name,organism_name&size=500`
- **Date de téléchargement** : 17/10/2026
- **Version de la base** : Release 2026
- **Conditions d'usage / licence** : Creative Commons Attribution (CC BY 4.0) — Usage académique et redistribution autorisés.
- **Nombre de lignes brutes** : 3 609 (3 608 séquences + 1 ligne d'en-tête)
- **Empreinte SHA-256 du fichier brut** : `ef21ecd472244cbd68d101651405ea90efd8f8e445fbe08b0769530dffe00134`
- **Fichier brut sauvegardé sous** : `data/raw/amp_database_20261017.tsv`
- **Téléchargement manuel ?** : Non (module `python -m peplib.download_t8_t9` avec pagination).
- **Notes** : Contient des séquences de toutes longueurs ; le filtrage sur la fenêtre de longueur [10, 50] et sur l'alphabet canonique des 20 acides aminés est appliqué à l'étape C (T12–T14).

---

## Source 3 : Bactériocines extraites d'UniProtKB — Profil Non-Pharmaceutique (T9)

- **Nom de la source** : UniProtKB (API REST), filtré par le mot-clé `KW-0078` (Bacteriocin)
- **Base dédiée de référence (non utilisée)** : BIOPEP-UWM / BACTIBASE n'ont pas été téléchargées cette semaine. UniProt a été retenu pour la même raison qu'en Source 2. **Écart assumé au cahier des charges T9**, à substituer par une base dédiée si le coordinateur l'exige.
- **URL de requête** : `https://rest.uniprot.org/uniprotkb/search?query=keyword:KW-0078 AND length:[5 TO 100]&format=tsv&fields=accession,sequence,length,protein_name,organism_name&size=500`
- **Date de téléchargement** : 17/10/2026
- **Version de la base** : Release 2026
- **Conditions d'usage / licence** : Creative Commons Attribution (CC BY 4.0) — Usage académique et redistribution autorisés.
- **Nombre de lignes brutes** : 787 (786 séquences + 1 ligne d'en-tête)
- **Empreinte SHA-256 du fichier brut** : `3832363fd7db973a4742a8dbea1335fccac19f65de4bccaed4f39d24e1ffceb4`
- **Fichier brut sauvegardé sous** : `data/raw/bactibase_bacteriocins_20261017.tsv`
- **Téléchargement manuel ?** : Non (module `python -m peplib.download_t8_t9` avec pagination).
- **Notes** : Cible les peptides de bioconservation alimentaire (bactériocines). Le filtrage strict sur la longueur 10–50 et les caractères non standards est effectué dans `peplib/filters.py`.

---

## Biosécurité — signalement des toxines (T10)

UniProt contient des toxines (conotoxines notamment). Conformément au cahier des
charges, elles sont **signalées** et doivent être exclues du corpus génératif ;
aucune séquence de type toxine ne sera diffusée.

- **Script de génération** : `python -m peplib.check_toxins` (simple point d'entrée) → logique dans `peplib/biosafety.py`, tests dans `tests/test_biosafety.py`.
- **Fichier de sortie** : `data/processed/toxin_flags.json` — une entrée `{"Entry": ..., "is_toxin": bool}` par séquence.
- **Méthode** : concaténation de toutes les colonnes textuelles de chaque entrée, puis recherche insensible à la casse du motif `Toxin|Conotoxin`.
- **Empreinte SHA-256 de la sortie** : `9c1182834dc9d1dbedfd9248704c5bd630f33e43ab7e39aac5087d9e5685a626`
- **Résultat sur le corpus T7** : 1 921 entrées signalées sur 12 423, dont **340/340 conotoxines** (`protein_name` contenant « conotoxin »).

### Reproductibilité des empreintes

`core.autocrlf` convertit par défaut les fins de ligne en CRLF sous Windows, ce
qui rendait les SHA-256 ci-dessus invérifiables sur un poste Windows. Le fichier
`.gitattributes` neutralise cette conversion pour `data/raw/*.tsv` et
`data/processed/*.json` : les empreintes sont désormais identiques sur tous les
systèmes. Les valeurs ci-dessus sont celles des fichiers tels qu'ils sont
versionnés (fins de ligne LF).

## Source APD6 — Natural Antimicrobial Peptides (2024)

- **Nom de la source :** Antimicrobial Peptide Database (APD6), University of Nebraska Medical Center.
- **URL de téléchargement :** https://aps.unmc.edu/assets/sequences/naturalAMPs_APD2024a.fasta
- **Page des conditions :** https://aps.unmc.edu/downloads
- **Date de téléchargement :** 22/10/2026
- **Version du jeu de données :** APD 2024a (hébergé sur le site APD6).
- **Format brut :** FASTA.
- **Nombre d'entrées annoncé :** 3306 (comptage local à vérifier).
- **Nom du fichier :** naturalAMPs_APD2024a.fasta
- **Chemin local :** data/raw/naturalAMPs_APD2024a.fasta
- **SHA-256 :** B2D30BBBE57A383220D5FA68FDA5A73728D2345F88585A78AA5F9E9546CD3D66
- **Téléchargement manuel :** Oui.
- **Conditions d'usage :** utilisation pour la recherche et l'enseignement avec citation de la source, selon les indications du site APD.
- **Décision :** GO pour l'utilisation académique locale ; redistribution publique non autorisée explicitement, donc suspendue en attendant clarification.
- **Gestion Git :** fichier brut exclu via .gitignore ; ne pas publier le corpus dérivé contenant APD sans clarification des droits.
- **Résidus modifiés :** analyse et règles de traitement prévues dans T34.
- **Biosécurité :** réappliquer les contrôles existants avant toute utilisation dans le corpus génératif.
  