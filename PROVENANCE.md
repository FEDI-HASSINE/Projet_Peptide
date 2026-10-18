# PROVENANCE.md — Corpus de peptides (Semaine 1, peplib v0.1)

Ce document consigne la **provenance de chaque fichier brut** de `data/raw/` :
source, requête, date, licence, nombre de lignes et **empreinte SHA-256**.
Un fichier brut ne doit jamais être modifié ; toute évolution passe par
`data/processed/`, régénérable par le code de `peplib/`.

**Vérification d'une empreinte** (à relancer après tout transfert) :

```bash
python -c "import hashlib;print(hashlib.sha256(open('data/raw/uniprot_reviewed_20261017.tsv','rb').read()).hexdigest())"
```

**Reproductibilité** : les trois fichiers se régénèrent depuis l'API REST UniProt :

```bash
python -m peplib.download_uniprot   # T7 — corpus général
python -m peplib.download_t8_t9     # T8 — AMP, T9 — bactériocines
```

> Les noms de fichiers `amp_database_...` et `bactibase_bacteriocins_...`
> sont **hérités** : le contenu réel est un sous-ensemble **UniProtKB**
> obtenu par mots-clés (voir chaque bloc). Voir la note « Écart au spec ».

---

# Source 1 : UniProtKB (corpus général)

- **Nom de la source** : UniProtKB (`reviewed:true`)
- **URL de requête** : `https://rest.uniprot.org/uniprotkb/search?query=length:[10 TO 50] AND reviewed:true&format=tsv&fields=accession,sequence,length,protein_name,keyword&size=500`
- **Date de téléchargement** : 17/10/2026
- **Version de la base** : non spécifiée (UniProt REST, interrogation du 17/10/2026)
- **Conditions d'usage / licence** : UniProt Terms of Use — données **CC BY 4.0**
- **Nombre de lignes brutes** : 12 424 (12 423 séquences + 1 ligne d'en-tête)
- **Empreinte SHA-256 du fichier brut** : `397959390c4ca2e6f55c5356e281f849398de0461328082f7549e2e1ad42f285`
- **Fichier brut sauvegardé sous** : `data/raw/uniprot_reviewed_20261017.tsv`
- **Téléchargement manuel ?** : Non (automatisé, `peplib.download_uniprot`)
- **Notes** : peptides de 10 à 50 résidus, uniquement entrées **reviewed** ;
  pagination par en-tête `Link`. **1 921 entrées signalées comme toxines**
  (conotoxines et autres) dans `data/processed/toxin_flags.json` — à exclure
  du corpus génératif et à ne jamais diffuser.

---

# Source 2 : peptides antimicrobiens (AMP)

- **Nom de la source** : UniProtKB, sous-ensemble par mots-clés
  (`KW-0929` Antimicrobial / `KW-0044` Antibiotic)
- **URL de requête** : `https://rest.uniprot.org/uniprotkb/search?query=(keyword:KW-0929 OR keyword:KW-0044) AND reviewed:true&format=tsv&fields=accession,sequence,length,protein_name,organism_name&size=500`
- **Date de téléchargement** : 17/10/2026
- **Version de la base** : non spécifiée
- **Conditions d'usage / licence** : UniProt Terms of Use — données **CC BY 4.0**
- **Nombre de lignes brutes** : 3 609 (3 608 séquences + 1 ligne d'en-tête)
- **Empreinte SHA-256 du fichier brut** : `ef21ecd472244cbd68d101651405ea90efd8f8e445fbe08b0769530dffe00134`
- **Fichier brut sauvegardé sous** : `data/raw/amp_database_20261017.tsv`
- **Téléchargement manuel ?** : Non (automatisé, `peplib.download_t8_t9`)
- **Notes** : profil pharmaceutique. Le nom de fichier `amp_database` est
  **hérité et trompeur** : il ne s'agit **pas** d'APD3 (source recommandée par
  le spec), mais d'un sous-ensemble UniProtKB par mots-clés. À confirmer ou
  corriger par le coordinateur.

---

# Source 3 : bactériocines (non-pharmaceutique)

- **Nom de la source** : UniProtKB, sous-ensemble par mot-clé (`KW-0078` Bacteriocin)
- **URL de requête** : `https://rest.uniprot.org/uniprotkb/search?query=keyword:KW-0078 AND length:[5 TO 100]&format=tsv&fields=accession,sequence,length,protein_name,organism_name&size=500`
- **Date de téléchargement** : 17/10/2026
- **Version de la base** : non spécifiée
- **Conditions d'usage / licence** : UniProt Terms of Use — données **CC BY 4.0**
- **Nombre de lignes brutes** : 787 (786 séquences + 1 ligne d'en-tête)
- **Empreinte SHA-256 du fichier brut** : `3832363fd7db973a4742a8dbea1335fccac19f65de4bccaed4f39d24e1ffceb4`
- **Fichier brut sauvegardé sous** : `data/raw/bactibase_bacteriocins_20261017.tsv`
- **Téléchargement manuel ?** : Non (automatisé, `peplib.download_t8_t9`)
- **Notes** : profil non-pharmaceutique (biopréservation). Le nom de fichier
  `bactibase_bacteriocins` est **hérité et trompeur** : il ne s'agit **pas** de
  BACTIBASE, mais d'un sous-ensemble UniProtKB par mot-clé.

---

# Écart au spec (à trancher)

Le spec (Semaine 1, §B) recommandait une base AMP dédiée (**APD3**, repli
DBAASP/DRAMP) et une base non-pharma (**BIOPEP-UWM** ou **BACTIBASE**).
Les trois corpus livrés sont en réalité des **sous-ensembles UniProtKB**
(T7/T8/T9). Conséquences :

- **Avantage** : licence unique et claire (CC BY 4.0), acquisitions
  reproductibles et déjà implémentées.
- **Inconvénient** : les trois sources partagent le même biais de base de
  données ; la « différence de source » mesurée (composition, JS) reflète les
  **filtres de requête**, pas des bases indépendantes. À assumer dans
  l'analyse T27 et à mentionner au séminaire T28.

Décision à prendre par le coordinateur : conserver UniProt (documenté ici)
ou ajouter APD3/BACTIBASE en sources additionnelles.

---

# Licence, biosécurité et partage

- **Licence** : l'ensemble du corpus brut est sous **UniProt Terms of Use /
  CC BY 4.0**. Vérifier les conditions de redistribution avant tout partage
  public du corpus ou de séquences dérivées.
- **Biosécurité** : UniProt contient des toxines (conotoxines). **1 921 entrées**
  sont signalées dans `data/processed/toxin_flags.json` ; elles sont **exclues
  du corpus génératif**. **Ne jamais diffuser de séquences générées de type
  toxine.**
- **Immuabilité** : `data/raw/` n'est jamais modifié ; les empreintes ci-dessus
  permettent de détecter toute altération.
