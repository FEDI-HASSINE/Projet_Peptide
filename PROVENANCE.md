# PROVENANCE.md — Corpus de peptides (Semaine 1 + extension Semaine 2)

**Fichier CANONIQUE** — c'est le seul document de provenance du dépôt.
Il consigne la **provenance de chaque fichier brut** de `data/raw/` :
source, requête, date, licence, nombre de lignes et **empreinte SHA-256**.
Un fichier brut ne doit jamais être modifié ; toute évolution passe par
`data/processed/`, régénérable par le code de `peplib/`.

(Le fichier `data/PROVENANCE.md` est réduit à un pointeur : consolidation
de l'emplacement autorisée lors de la relecture de la PR #10.)

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

**Reproductibilité des empreintes (CRLF)** : `core.autocrlf` convertit par
défaut les fins de ligne en CRLF sous Windows, ce qui rendait les SHA-256
invérifiables sur un poste Windows. Le fichier `.gitattributes` neutralise
cette conversion pour `data/raw/*.tsv` et `data/processed/*.json` : les
empreintes sont identiques sur tous les systèmes. Les valeurs ci-dessous sont
celles des fichiers tels qu'ils sont versionnés (fins de ligne LF).

---

## Source 1 : UniProtKB (corpus général — T7)

- **Nom de la source** : UniProtKB (`reviewed:true` / Swiss-Prot)
- **URL de requête** : `https://rest.uniprot.org/uniprotkb/search?query=length:[10 TO 50] AND reviewed:true&format=tsv&fields=accession,sequence,length,protein_name,keyword&size=500`
- **Date de téléchargement** : 17/10/2026
- **Version de la base** : UniProt Release 2026 (API REST)
- **Conditions d'usage / licence** : Creative Commons Attribution (CC BY 4.0) — UniProt Terms of Use (consulté le 17/10/2026 ; usage académique et redistribution autorisés avec citation).
- **Nombre de lignes brutes** : 12 424 (12 423 séquences + 1 ligne d'en-tête)
- **Empreinte SHA-256 du fichier brut** : `397959390c4ca2e6f55c5356e281f849398de0461328082f7549e2e1ad42f285`
- **Fichier brut sauvegardé sous** : `data/raw/uniprot_reviewed_20261017.tsv`
- **Téléchargement manuel ?** : Non (module `python -m peplib.download_uniprot` avec pagination via en-tête `Link`).
- **Notes** : peptides de 10 à 50 résidus, uniquement entrées reviewed.
  Les colonnes `protein_name` et `keyword` sont récupérées pour identifier
  et signaler/exclure les toxines (ex. conotoxines, `KW-0800`) du corpus
  génératif — jamais de séquences toxines diffusées.

---

## Source 2 : Peptides antimicrobiens extraits d'UniProtKB — profil pharmaceutique (T8)

- **Nom de la source** :
  UniProtKB (API REST), filtré par les mots-clés `KW-0929` (Antimicrobial) et `KW-0044` (Antibiotic), `reviewed:true`
- **Base dédiée de référence (non utilisée en Semaine 1)** : APD3 / DBAASP / DRAMP
- **URL de requête** : `https://rest.uniprot.org/uniprotkb/search?query=(keyword:KW-0929 OR keyword:KW-0044) AND reviewed:true&format=tsv&fields=accession,sequence,length,protein_name,organism_name&size=500`
- **Date de téléchargement** : 17/10/2026
- **Version de la base** : Release 2026
- **Conditions d'usage / licence** : Creative Commons Attribution (CC BY 4.0) — usage académique et redistribution autorisés.
- **Nombre de lignes brutes** : 3 609 (3 608 séquences + 1 ligne d'en-tête)
- **Empreinte SHA-256 du fichier brut** : `ef21ecd472244cbd68d101651405ea90efd8f8e445fbe08b0769530dffe00134`
- **Fichier brut sauvegardé sous** : `data/raw/amp_database_20261017.tsv`
- **Téléchargement manuel ?** : Non (module `python -m peplib.download_t8_t9` avec pagination).
- **Notes** : contient des séquences de toutes longueurs ; le filtrage sur la
  fenêtre de longueur [10, 50] et sur l'alphabet canonique des 20 acides
  aminés est appliqué à l'étape C (T12–T14). Le nom de fichier
  `amp_database_...` est **hérité et trompeur** : il ne s'agit **pas** d'APD3,
  mais d'un sous-ensemble UniProtKB. **Écart assumé au cahier des charges
  T8**, en cours de résolution en Semaine 2 (blocs APD6 et AMP Challenge
  ci-dessous).

---

## Source 3 : Bactériocines extraites d'UniProtKB — profil non-pharmaceutique (T9)

- **Nom de la source** : UniProtKB (API REST), filtré par le mot-clé `KW-0078` (Bacteriocin)
- **Base dédiée de référence (non utilisée en Semaine 1)** : BIOPEP-UWM / BACTIBASE
- **URL de requête** : `https://rest.uniprot.org/uniprotkb/search?query=keyword:KW-0078 AND length:[5 TO 100]&format=tsv&fields=accession,sequence,length,protein_name,organism_name&size=500`
- **Date de téléchargement** : 17/10/2026
- **Version de la base** : Release 2026
- **Conditions d'usage / licence** : Creative Commons Attribution (CC BY 4.0) — usage académique et redistribution autorisés.
- **Nombre de lignes brutes** : 787 (786 séquences + 1 ligne d'en-tête)
- **Empreinte SHA-256 du fichier brut** : `3832363fd7db973a4742a8dbea1335fccac19f65de4bccaed4f39d24e1ffceb4`
- **Fichier brut sauvegardé sous** : `data/raw/bactibase_bacteriocins_20261017.tsv`
- **Téléchargement manuel ?** : Non (module `python -m peplib.download_t8_t9` avec pagination).
- **Notes** : cible les peptides de bioconservation alimentaire
  (bactériocines). Le nom de fichier `bactibase_bacteriocins_...` est
  **hérité et trompeur** : il ne s'agit **pas** de BACTIBASE, mais d'un
  sous-ensemble UniProtKB.

---

## Biosécurité — signalement des toxines (T10)

UniProt contient des toxines (conotoxines notamment). Conformément au
cahier des charges, elles sont **signalées** et doivent être exclues du
corpus génératif ; aucune séquence de type toxine ne sera diffusée.

- **Script de génération** : `python -m peplib.check_toxins` (point d'entrée) → logique dans `peplib/biosafety.py`, tests dans `tests/test_biosafety.py`.
- **Fichier de sortie** : `data/processed/toxin_flags.json` — une entrée `{"Entry": ..., "is_toxin": bool}` par séquence.
- **Méthode** : concaténation de toutes les colonnes textuelles de chaque entrée, puis recherche insensible à la casse du motif `Toxin|Conotoxin`.
- **Empreinte SHA-256 de la sortie** : `9c1182834dc9d1dbedfd9248704c5bd630f33e43ab7e39aac5087d9e5685a626`
- **Résultat sur le corpus T7** : 1 921 entrées signalées sur 12 423, dont **340/340 conotoxines** (`protein_name` contenant « conotoxin »).
- Le flag est **relancé sur le corpus v0.2** (avec APD6 et AMP Challenge) avant toute utilisation générative.

---

## Écart au spec (Semaine 1) et résolution (Semaine 2)

Le spec (Semaine 1, §B) recommandait une base AMP dédiée (**APD3**, repli
DBAASP/DRAMP) et une base non-pharma (**BIOPEP-UWM** ou **BACTIBASE**).
Les trois corpus livrés sont en réalité des **sous-ensembles UniProtKB**
(T7/T8/T9). Conséquences :

- **Avantage** : licence unique et claire (CC BY 4.0), acquisitions
  reproductibles et déjà implémentées.
- **Inconvénient** : les trois sources partagent le même biais de base de
  données ; la « différence de source » mesurée (composition, JS) reflète les
  **filtres de requête**, pas des bases indépendantes.

**Résolution (Semaine 2, T32)** : le coordinateur a accepté l'ajout de deux
vraies bases AMP, documentées dans les deux blocs ci-dessous : **APD6**
(listes 2024 de la base APD3) et la base **AMP Challenge 2027**
(agrégat de 13 bases AMP, incluant DBAASP et DRAMP — les replis recommandés).
La différence de composition UniProt ↔ APD6 ↔ AMPChallenge est devenue
mesurable (T38, recouvrement trilatéral).

---

## Source APD6 — Antimicrobial Peptide Database, listes 2024 (Semaine 2, T32)

- **Nom de la source** : Antimicrobial Peptide Database, University of Nebraska Medical Center (`aps.unmc.edu`)
- **Nomenclature APD3 / APD6** : APD3 est la base de référence (Wang et al., 2016) ; les listes 2024 hébergées sur le site sont appelées **APD6** ; l'ensemble est citable indifféremment sous les deux dénominations.
- **URL de téléchargement** : `https://aps.unmc.edu/assets/sequences/naturalAMPs_APD2024a.fasta`
- **Page des conditions** : `https://aps.unmc.edu/downloads`
- **Date de téléchargement** : 22/10/2026
- **Version du jeu de données** : APD 2024a (hébergée sur le site APD6)
- **Format brut** : FASTA
- **Nombre de lignes brutes** : 9 919 (comptage local du 22/10/2026)
- **Nombre de séquences (comptage local, mesuré)** : 3 306 séquences réelles,
  + 1 pseudo-entrée d'annonce en tête de fichier
  (`>Your search led to 3306 peptides` — à exclure au connecteur T34 ;
  0 séquence vide, 0 en-tête dupliqué, alphabet 100 % canonique)
- **Empreinte SHA-256 du fichier brut** : `b2d30bbbe57a383220d5fa68fda5a73728d2345f88585a78aa5f9e9546cd3d66`
- **Fichier brut sauvegardé sous** : `data/raw/naturalAMPs_APD2024a.fasta`
- **Téléchargement manuel ?** : Oui (le site refuse les clients sans en-tête navigateur, code 406)
- **Conditions d'usage / licence** — citation textuelle du site (consulté le 22/10/2026) : *« Term of use: If you have downloaded and used any APD list or info for your research and education, you can cite the APD3 or APD6. »* (page `https://aps.unmc.edu/downloads`)
- **Décision (GO / NO-GO)** : GO pour l'utilisation académique locale ;
  redistribution publique **suspendue** — les conditions du site exigent la
  citation mais n'autorisent pas explicitement la redistribution ; tout
  partage du corpus v0.2 contenant des séquences APD6 attend la clarification
  des droits.
- **Gestion Git** : fichier brut exclu via `.gitignore` (non commité).
- **Résidus modifiés** : analyse et règles de traitement prévues dans T34.
- **Biosécurité** : flag toxines réappliqué avant toute utilisation générative.

---

## Source AMP Challenge 2027 — antibacterial.fasta (Semaine 2, T32)

- **Nom de la source** : AMP Challenge 2027 (donnée fournie par l'encadrant)
- **Dépôt** : `szczurek-lab/amp-challenge-2027` (licence du dépôt : BSD-3-Clause)
- **URL de téléchargement** : `https://raw.githubusercontent.com/szczurek-lab/amp-challenge-2027/main/data/antibacterial.fasta`
- **Date de téléchargement** : 22/10/2026 (re-vérifié le même jour par le relecteur)
- **Version** : branche `main`, identifiée par son empreinte (le fichier peut évoluer ; le SHA-256 fait foi)
- **Format brut** : FASTA
- **Nombre de lignes brutes** : 78 896 (comptage local du 22/10/2026)
- **Nombre de séquences brutes** : 39 448 (0 doublon exact, alphabet 100 % canonique)
- **Fichier brut sauvegardé sous** : `data/raw/challenge_antibacterial_20261022.fasta`
- **Empreinte SHA-256 du fichier brut** : `cbbeac64ba95746d87961e8ad9dd0849ae8058d15a300b2e7f6990730ca521e9`
- **Téléchargement manuel ?** : Non (URL directe)
- **Licence et redistribution** :
  - licence du dépôt : BSD-3-Clause ;
  - le fichier agrège **13 bases amont** : dbAMP, DRAMP, DBAASP, CAMP,
    SATPdb, APD, AMPDB, DADP, InverPep, CancerPPD, BaAMPs, CyBase, ParaPep
    (fragments contenues dans le champ `dbs=` des en-têtes FASTA) ;
  - les licences des bases amont restent à vérifier individuellement **avant
    tout partage** du corpus v0.2 ;
  - **Décision** : GO provisoire pour l'analyse académique locale ; NO-GO
    provisoire pour la redistribution publique ; le fichier brut reste local,
    exclu du suivi Git via `.gitignore`.
- **Métadonnées d'en-têtes** : `len`, `charge`, `disulfide`, `activity`, `dbs` —
  parsées au connecteur T34 ; le champ `dbs=` est conservé en colonne
  `dbs_amont` mais ne pilote aucun calcul en Semaine 2.
- **Filtrage attendu (à vérifier en T34)** : 34 474 séquences dans la fenêtre 10–50.
- **Biosécurité** : contrôle des toxines réappliqué avant utilisation générative.

---

## Licence, biosécurité et partage (synthèse)

- **UniProt (3 fichiers)** : CC BY 4.0 ; vérifier les conditions de
  redistribution avant tout partage public du corpus ou de séquences dérivées.
- **APD6** : usage recherche + enseignement **avec citation** (conditions
  citées ci-dessus) ; **redistribution publique suspendue**.
- **AMP Challenge** : dépôt BSD-3-Clause ; licences des 13 bases amont à
  vérifier case à case ; **NO-GO redistribution provisoire**.
- **Biosécurité** : le flag toxines (`toxin_flags.json`) est relancé sur le
  corpus v0.2 fusionné. **Ne jamais diffuser de séquences générées de type toxine.**
- **Immuabilité** : `data/raw/` n'est jamais modifié ; les empreintes
  ci-dessus permettent de détecter toute altération.
- **Partage du corpus v0.2** : ne distribuer qu'après recalcul des empreintes
  et vérification de chaque licence amont ; ce document est mis à jour
  **avant** le partage.

## T34 - Normalisation APD et AMP Challenge

- APD : naturalAMPs_APD2024a.fasta (site APD6).
- AMP Challenge : challenge_antibacterial_20261022.fasta.
- Format : FASTA, lecture sans modification des fichiers bruts.
- Filtrage : longueur 10 a 50 et alphabet des 20 acides amines canoniques.
- Residus X/B/Z/J, U/O et autres caracteres non canoniques rejetes sans conversion arbitraire.
- Les modifications chimiques ne sont pas deduites des sequences FASTA.
- APD : 3306 sequences brutes, 2796 conservees.
- AMP Challenge : 39448 sequences brutes, 34474 conservees.
- Fichiers bruts exclus de Git.

