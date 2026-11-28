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

## Corpus v0.2 — data/processed/corpus_v0.2.csv (Semaine 2, T37)

- **Version** : v0.2
- **Date de génération** : 24/10/2026
- **Commande** : `python -m peplib.build_corpus_v02`
- **Entrées** : les cinq CSV de `data/processed/` (`uniprot.csv`, `amp.csv`,
  `bactibase.csv`, `amp_apd3.csv`, `challenge.csv`), eux-mêmes régénérables
  depuis `data/raw/` par les connecteurs T7/T8/T9/T34.
- **Atomes de source (T37)** : `UniProt`, `APD3`, `AMPChallenge`. Les trois
  fichiers `uniprot.csv`, `amp.csv` et `bactibase.csv` sont trois requêtes
  UniProtKB (cf. sources 1 à 3) : ils alimentent l'atome `UniProt` et ne sont
  pas comptés comme des sources indépendantes.
- **Comptes** : 50 352 entrées, 42 575 séquences uniques, 7 777 occurrences
  redondantes (15,45 %), dont 2 984 doublons internes à un fichier d'entrée et
  4 793 doublons inter-sources. 3 244 séquences sont présentes dans au moins
  deux atomes (2 772 dans deux, 472 dans les trois). La somme des
  `n_occurrences` vaut 50 352.
- **Recouvrement (T38)** : `data/processed/overlap_t38.csv` (trois paires) et
  `data/processed/overlap_t38_triple.csv` (triple recouvrement, 472 séquences,
  Jaccard 0,011086) ; figure `results/overlap_t38.png`. Les quasi-doublons
  relèvent de T42 et ne sont pas traités ici.
- **Métadonnées disponibles par source** :
  - UniProt (`uniprot.csv`, `amp.csv`, `bactibase.csv`) : `sequence`, `id` (accession) ;
  - APD3 (`amp_apd3.csv`) : `sequence`, `id` (APD) ; la cible et l'activité
    APD3 ne sont pas extraites par le connecteur T34 — perte héritée, non
    restituable par le corpus, à corriger dans le connecteur ;
  - AMPChallenge (`challenge.csv`) : `charge`, `disulfide`, `activity`,
    `dbs_amont` conservées ; `longueur` est recalculée pour toutes les lignes.
- **Empreinte SHA-256** : `f01e248c93179b70773b1775991619f3c4cb0c0975420ba23bdddacc71ae88da`
  (fichier en LF, sans BOM ; `data/processed/*.csv` est marqué `-text` dans
  `.gitattributes` pour que l'empreinte soit identique sur tout système)
- **Restriction de diffusion** : ce fichier contient des séquences APD6
  (redistribution publique suspendue) et AMP Challenge (NO-GO provisoire).
  Il reste **réservé au dépôt privé** : ne pas publier, ne pas joindre à une
  prépublication ni à un dépôt public sans vérification des licences amont.

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

## T34 — Normalisation APD et AMP Challenge

- Sources : `data/raw/naturalAMPs_APD2024a.fasta` (site APD6) et
  `data/raw/challenge_antibacterial_20261022.fasta` (AMP Challenge).
- Format : FASTA, **lecture seule** des fichiers bruts (`utf-8-sig`) ;
  aucune écriture dans `data/raw/`.
- Filtrage : `classify_sequence` de `peplib.filters`, déjà utilisé en
  semaine 1 (alphabet des 20 acides aminés canoniques, longueur 10–50,
  une seule cause de rejet par séquence, ordre des causes inchangé).
- **Résidus non canoniques : 0 occurrence observée** dans les deux bruts
  (alphabet vérifié caractère par caractère). La règle de rejet
  `X B Z J`, puis `U O`, puis autres caractères n'a donc coûté aucune
  séquence sur ces deux fichiers.
- **Amidation : 0 marqueur observé** dans le FASTA APD. Aucune colonne
  d'annotation n'est créée. Convention retenue pour une source ultérieure
  qui en contiendrait : annotation dans une colonne dédiée, symbole retiré
  de la séquence **avant** filtrage, convention revalidée à ce moment-là.
- **Métadonnées APD3 : dérogation actée le 24/10/2026.** Le FASTA ne
  contient que `id` et `sequence` ; les colonnes `nom, organisme, cible,
  activité, unité` exigées par la spécification T34 ne peuvent pas être
  remplies depuis ce fichier. Aucune colonne n'est fabriquée à vide ;
  l'ajout se fera depuis l'export tabulaire APD, au plus tard au mois 3,
  après accord du coordinateur.
- Métadonnées AMP Challenge : `len`, `charge`, `disulfide`, `activity` et
  `dbs_amont` extraites des en-têtes et renseignées à **100 %** des
  39 448 enregistrements bruts. `dbs_amont` conserve les 13 bases amont
  possibles ; c'est une colonne d'information, elle ne pilote ni filtre ni
  split (décision T39).
- Sorties : `data/processed/amp_apd3.csv` et
  `data/processed/challenge.csv`, colonnes `sequence, source, id`
  (+ métadonnées pour le Challenge), régénérables par
  `python -m peplib.download_apd3` et
  `python -m peplib.download_challenge`.
- Comptages mesurés : APD **3 306** séquences brutes → **2 796**
  conservées (82 trop courtes, 428 trop longues) ; AMP Challenge
  **39 448** brutes → **34 474** conservées (4 974 trop courtes).
  Invariants « conservées + rejetées = brutes » vérifiés pour les deux.
- Fichiers bruts exclus de Git (`.gitignore`).

## T35 — Clusters à 90 % d'identité (corpus v0.1)

- **Artefact** : `data/processed/clusters_v01_90.csv`
- **Commande** : `python -m peplib.cluster`
- **Entrée** : `data/processed/corpus_filtered_all.csv` (v0.1, 13 082 lignes,
  9 021 séquences uniques).
- **Contenu** : une ligne par séquence unique ; colonnes
  `sequence, sequence_id, cluster_id, source` (`sequence_id` et `source` =
  première occurrence dans le corpus). L'artefact est autosuffisant : aucune
  jointure et aucun ordre de lignes ne sont nécessaires.
- **Définition** : composantes connexes du graphe « identité ≥ 90 % »
  (Union-Find avec compression de chemin et union par rang) ; paires
  candidates issues de l'index inversé de 3-mers, puis distance de
  Levenshtein bornée. Le seuil est calculé exactement (`Fraction`), sans
  epsilon flottant.
- **Comptes** : 9 021 séquences → 6 855 clusters (2 166 fusions), plus grande
  famille 62, médiane 1 ; 6 422 239 paires candidates = 4 751 008 écartées par
  longueur + 1 671 231 comparées.
- **Chaînage** : un cluster est une composante connexe et peut donc contenir
  deux séquences d'identité inférieure à 90 %. Mesure sur v0.1 : 229 des 822
  clusters multi-membres (27,9 %) contiennent une paire sous 90 %, 12 (1,5 %)
  sous 70 %, minimum 45,5 % (cluster 1892). **T39 traite le cluster comme
  l'unité indivisible du split.**
- **Portée** : clusters calculés sur le corpus v0.1 uniquement. Le corpus v0.2
  (42 575 séquences uniques) en contient 33 554 (78,8 %) absentes :
  **recalculer les clusters sur v0.2 après T37 avant tout split**.
- **Empreinte SHA-256** : `6af381c2ac9e9ee2e4a9ccafb991009c4732ca6f6f4e455b35805993fed34037` (fichier LF, sans BOM ;
  `data/processed/*.csv` marqué `-text` dans `.gitattributes`).
- **Figure** : `results/cluster_sizes_v01.png`.


---

## T39 — Vues par source du corpus v0.2

- **Date du contrôle** : 25/10/2026
- **Commande** : `python -m peplib.source_views_v02`
- **Entrée** : `data/processed/corpus_v0.2.csv`
- **Méthode** : sélection des séquences selon les atomes de source `UniProt`, `APD3` et `AMPChallenge`. Les combinaisons de sources séparées par `+` ou `|` sont prises en charge.
- **Sorties** :
  - `data/processed/view_uniprot_v02.csv` (SHA-256 : `73d783156b45b511805f7454d40ddcc6304ee815d609caaaec86368a7475c507`)
  - `data/processed/view_apd3_v02.csv` (SHA-256 : `8703517a27d0417dfb6eb1d5107805f59ab3dae034091a9ad677ea34bde2fd26`)
  - `data/processed/view_ampchallenge_v02.csv` (SHA-256 : `62cc0ac3a45a640070ec62f282a445779e21b96fd5164a8782553b6c2ce061e1`)

### Résultats

| Source | Séquences |
|---|---:|
| UniProt | 9 021 |
| APD3 | 2 796 |
| AMPChallenge | 34 474 |
| **Total des appartenances** | **46 291** |
| **Séquences uniques du corpus** | **42 575** |

Les 42 575 séquences du corpus sont couvertes. Une séquence présente dans plusieurs sources apparaît dans chacune des vues correspondantes. La différence entre 46 291 appartenances et 42 575 séquences uniques provient des recouvrements.

- **Validation** : 5 tests spécifiques réussis dans `tests/test_source_views_v02.py`.
- **Diffusion** : les vues contenant des séquences APD6 ou AMP Challenge ne doivent pas être redistribuées publiquement avant clarification des licences.

---

## Biosécurité — Signalement des toxines sur le corpus v0.2

- **Date du contrôle** : 25/10/2026
- **Commande** : `python -m peplib.check_toxins_v02`
- **Entrée** : `data/processed/corpus_v0.2.csv`
- **Annotations utilisées** : `data/raw/uniprot_reviewed_20261017.tsv`
- **Sortie** : `data/processed/toxin_flags_v02.csv`
- **Méthode** : jointure par séquence exacte avec UniProt et recherche insensible à la casse de `Toxin|Conotoxin|Venom|Toxic` dans `Protein names` et `Keywords`. Le motif est volontairement sur-inclusif (une description « non-toxic » est aussi signalée) : le mode d'échec est sûr. Les annotations positives sont agrégées par séquence. Correction relecture PR #19 (B2) : « Venom » et « Toxic » ajoutés — 50 séquences du corpus portant des annotations de venin/toxicité étaient précédemment `not_flagged_by_uniprot`.

### Résultats

| Statut | Nombre |
|---|---:|
| Séquences totales | 42 575 |
| Signalées (`flagged`) | 1 786 |
| Non signalées par UniProt (`not_flagged_by_uniprot`) | 7 125 |
| Sans annotation UniProt (`unknown`) | 33 664 |

**Invariant** : 1 786 + 7 125 + 33 664 = 42 575.

### Limites

- Ce contrôle détecte des annotations textuelles, pas une toxicité biologique démontrée.
- `not_flagged_by_uniprot` ne signifie pas « non toxique ».
- `unknown` signifie qu'aucune annotation UniProt correspondante n'est disponible.
- Les annotations de toxicité d'APD3 et d'AMP Challenge ne sont pas couvertes par cette méthode.
- 110 séquences d'origine UniProt du corpus ne figurent pas dans le TSV d'annotations (accessions vérifiées absentes, ni par `Entry`, ni par séquence) : leur statut est forcé à `unknown`. À réconcilier (requête REST par accession ou alignement du corpus sur le raw) **avant le clustering v0.2** (ajustements A5/A6).
- L'ancien fichier `data/processed/toxin_flags.json` reste inchangé.

### Reproductibilité

- **SHA-256 de `toxin_flags_v02.csv`** : `171f150b7a41871c2d10b762756d8cf7ed1aae4fbbd69981557ec6056180d314` (fichier LF, sans BOM ;
  `data/processed/*.csv` marqué `-text` dans `.gitattributes`).
- **Tests spécifiques** : 11 réussis dans `tests/test_check_toxins_v02.py` (dont « Venom peptide » flague)
  ; le motif élargi est verrouillé par 3 tests ajoutés dans `tests/test_biosafety.py` ; les artefacts v0.2
  (comptes et empreinte) sont verrouillés par 4 tests dans `tests/test_artefacts_v02.py` (style B10).
- **Suite complète** : 369 tests réussis (Python 3.14 et 3.11).
- **Restriction** : ce contrôle ne constitue pas une validation de sécurité pour un usage thérapeutique.
- **Clause de consommation** : tout consommateur du corpus génératif v0.2 (splits A6, génération) doit joindre `toxin_flags_v02.csv` par `sequence` et traiter `unknown` comme `flagged` (mode d'échec sûr).


## Clusters v0.2 au critère de couverture 0,90 (A5, Semaine 2)

- **Date du calcul** : 26/10/2026
- **Commande** : `python -c "from peplib.cluster_v02 import run; run()"`
- **Entrée** : `data/processed/corpus_v0.2.csv` (42 575 séquences uniques,
  triées pour la numérotation).
- **Méthode** : composantes connexes (union-find, `peplib/cluster.py`) du
  graphe « il existe un alignement local SW (o = 11, e = 1, BLOSUM62) avec
  identité ≥ 0,90 ET couverture ≥ 0,90 ». Préfiltre multiset borné
  (occurrences intactes ≥ ⌊0,63·m − 3⌋, preuve écrite dans
  `Docs/A6_preuve_borne_couverture.md` ; la formule `max_distance` de la PR
  #16 n'est PAS réutilisée pour la couverture). Nombre de 3-mers intact
  exigé dépend de m = min(|a|, |b|) ; aucun epsilon ; comparaison au seuil
  exacte car le pas des ratios (dénominateurs ≤ 200) écrase l'erreur
  flottante (docstring `peplib/cluster_v02.py`).
- **Comptes mesurés** : 116 206 731 couples partageant un 3-mer →
  1 061 432 candidats après coupure sonore → 1 061 432 alignements SW →
  21 902 fusions. **T13 : 20 673 familles, taille max 4 541
  (10,6659 % du corpus), 17 592 singletons.**
- **Sorties** : `data/processed/clusters_v02_90.csv` (colonnes
  `sequence, sequence_id, cluster_id, sources`, LF sans BOM) ;
  `results/week2/clusters_v01_t13.csv` renommé `results/week2/clusters_v02_t13.csv` ;
  figure `results/week2/cluster_sizes_v02.png`.
- **Écart assumé v0.1 → v0.2** : la colonne `source` (atome unique) devient
  `sources` (fusion des atomes, séparateur `+`) : le corpus v0.2 est
  multi-sources ; le contrat de colonnes est identique par ailleurs.
- **Variante demandée non verrouillée** (seuil 0,70, identité globale,
  (o,e) alternatifs) : voir PR — budget de calcul arbitré le 26/10.
- **Empreintes SHA-256** :
  - `clusters_v02_90.csv` : `e9de2eb5f92fccb40d71b606b73280ade3f5a1461d88610301077bf5af3eed1f`
  - `clusters_v02_t13.csv` : `209ffc8d9ae5608e2e446ad47865140224c16b435ca01e0b9deeb7cb8fd57d7e`
- **Verrou de test** : `tests/test_artefacts_clusters_v02.py` (5 tests).


## Représentants par famille (T42, Semaine 2)

- **Date du calcul** : 26/10/2026
- **Commande** : `python -c "from peplib.cluster_v02 import representants; representants()"`
- **Entrée** : `data/processed/clusters_v02_90.csv` (artefact ci-dessus).
- **Critère du représentant (écrit, pas improvisé)** : la séquence la plus
  longue de la famille ; à égalité, l'ordre alphabétique croissant.
- **Sources** : fusion des atomes de TOUS les membres de la famille
  (dédoublonnés, ordre canonique UniProt < APD3 < AMPChallenge,
  séparateur `+`).
- **Comptes mesurés** : 20 673 lignes (1 par famille), 20 673 séquences
  uniques, somme des `taille_famille` = 42 575 (couverture complète du
  corpus). Colonnes : `sequence, sequence_id, cluster_id, sources,
  taille_famille`.
- **Empreinte SHA-256** :
  `931d8d7710183bd222630f7063861669baa69a4facb1729340c2b7db80998018`


## Splits v0.2 — 90/5/5 par familles et contrôle aléatoire (T39/A6, Semaine 2)

- **Date du calcul** : 26/10/2026
- **Commande** : `python -c "from peplib.splits import run; run()"`
- **Entrées** : `clusters_v02_90.csv` + `corpus_v0.2.csv` +
  `toxin_flags_v02.csv` (jointure par `sequence` ;
  `toxicity_status` suit la clause W4 de la PR #19 : `unknown` traité
  comme `flagged` pour tout usage génératif ; **0 séquence du corpus v0.2
  absente du TSV de toxines — la limite G4 de 110 UniProt ne touche pas
  les sorties de split ; elle est documentée dans la section corpus**).
- **Méthode** : affectation de familles ENTIÈRES (aucune à cheval),
  familles > 1 % du corpus systématiquement en train (la plus grande,
  4 541 > 425,75, y est), au moindre coût : dépassement de taille et
  dépassement par atome de source au carré (stratification par source,
  mixtures mesurées ci-dessous) ; tolérance ±0,5 pt par partie ; seed 42.
- **Comptes mesurés** (cible 38 317,5 / 2 128,75 / 2 128,75 ± 212,9) :
  train 38 336 (+18,5) ; val 2 095 (−33,75) ; test 2 144 (+15,25).
  20 665 familles mobilisées ; chevauchement de familles entre parts :
  **0 / 0 / 0**.
- **Mixtures par atome mesurées** (global : AMPChallenge 33 974,
  UniProt 9 021, APD3 2 796) :
  - train : AMPChallenge 30 996, UniProt 8 118, APD3 2 516 ;
  - val : 1 743 / 459 / 141 (attendu ≈ 1 699 / 451 / 140) ;
  - test : 1 735 / 444 / 139.
- **Split aléatoire de contrôle** : `data/splits_control/`, mêmes tailles
  exactes, seed 42, mêmes colonnes, `cluster_id` conservé (il traverse les
  parts par construction).
- **Fraction par source des statuts toxines (train)** : unknown 29 517,
  not_flagged_by_uniprot 7 017, flagged 1 785.
- **Empreintes SHA-256** (recouvrées du manifeste `splits_manifest.json`,
  SHA complet `670a896e625d65932a3f2ea7bdc56e503ac222272b8eceb2a2cff4b921b083cd`) :
  - `data/splits/train_v02.csv` : `f9016bc038ed473870e76e4a849962813ef010ac3079c31e227309aa4349e736`
  - `data/splits/val_v02.csv` : `f0a17569bf62be16ab931a1f8b47f8c62baad27f509bb62e816aad59f7dbb3a5`
  - `data/splits/test_v02.csv` : `b3c82a07673ab2ab534666103405186da5b1f4217693acf5f811afeb53e62cc5`
- **Verrou anti-dérive** : `tests/test_splits.py::test_run_sur_artefacts_synthetiques`
  (toy) + lecture du manifeste dans la PR ; `data/splits/*.csv` et
  `data/splits_control/*.csv` marqués `-text` dans `.gitattributes`.


## Courbe de fuite — val contre train (A6/D3, Semaine 2)

- **Date du calcul** : 26/10/2026
- **Commande** : `python -c "from peplib.fuite import run; run(tau_candidats=<t>)"` (t = 0,90 puis 0,70)
- **Méthode** : pour chaque séquence de val, meilleure paire (identité ,
  couverture) contre tout train, parmi les candidats du préfiltre multiset
  borné (A6) évalués uniquement en croisement val -> train ; SW complet
  pour chaque candidat retenu. F(θ) = proportion des séquences de val dont
  LA MÊME meilleure paire atteint identité ≥ θ ET couverture ≥ θ (le min
  de la meilleure paire, jamais un mélange de deux paires — fixé par un
  test). Un seul passage à tau = 0,70 couvre exactement F(0,90) ET F(0,70)
  (bornes K(m) croissantes en tau).
- **Pas de true biologique : les statuts toxines accompagnent les sorties,
  aucune décision d'exclusion n'est prise ici (§6 ouvert).**
- **Comptes mesurés (passe tau = 0,90)** : candidats val-train 17 916
  (split par familles) contre 95 816 (contrôle aléatoire) — le préfiltre
  mesure lui-même la séparation des familles. Résultats exacts :
  - split par familles : **F(0,90) = 0,000000** (0 val / 2 095 ;
    étanchéité vérifiée par la courbe, au-delà du test de familles) ;
  - contrôle aléatoire : F(0,90) = 0,557995, F(0,70) = 0,699761 —
    le contrôle fuit, comme attendu.
  - F(0,70) du split par familles mesuré sur le stock tau = 0,90 :
    0,333174 — **borne inférieure** ; remplacé ci-dessous par la mesure
    exacte.
- **Passe exacte tau = 0,70 (résultats complets, les deux splits)** :
  - candidats croisés évalués : 6 710 520 (split par familles),
    7 808 555 (contrôle aléatoire, exécution achevée le 27/10/2026) ;
  - **split par familles : F(0,90) = 0,000000 exact (récupéré avec le
    catalogue élargi) ET F(0,70) = 0,6071599** (0,607 = densité de cousins
    à ≥ 70 %/70 %, normalement élevée dans un corpus AMP ; le critère de
    fuite reste le seuil (0,90 ; 0,90), qui donne zéro) ;
  - **contrôle aléatoire (passe exacte achevée)** : F(0,90) = 0,5579952
    — identique au témoin de la passe tau = 0,90 (0,557995), cohérence
    confirmée — et **F(0,70) = 0,8252983 exact** (remplace la borne
    inférieure 0,699761). Verbatim du journal d'exécution :
    « Paires évaluées : 7800000/7808555 » puis
    « Split controle : {'F_70': 0.8252983293556085,
    'F_90': 0.5579952267303102} ». Le contrôle fuit, comme attendu ;
    la reprise de calcul n'a réexécuté QUE la partie contrôle (les faits
    du split familles, terminé le 26/10, sont réinjectés via
    ``FAITS_JOURNAL_CLUSTERS`` dans `peplib/fuite.py`).
- **Sorties** : `results/week2/fuite.csv` (détail par séquence, split
  contrôle de la passe exacte — le détail du split familles est couvert
  par sa passe complète du 26/10 consignée ci-dessus),
  `results/week2/fuite_fractions.csv` (4 lignes : F(0,70) et F(0,90)
  de chaque split), `results/week2/courbe_fuite.png`. Le journal
  d'exécution complet est conservé HORS DÉPÔT (règle `*.log` de
  `.gitignore`) ; les chiffres-témoins des deux journaux sont cités
  verbatim ci-dessus (« Paires évaluées : 6 710 520 » pour la partie
  familles, « Split controle : {'F_70': 0.8252983293556085,
  'F_90': 0.5579952267303102} » pour la partie contrôle achevée).
## A9 : mesure d'invariants sur paires reelles — results/align_invariants_a9.csv (PR #20)

- **Fichier** : `results/align_invariants_a9.csv` (500 lignes de donnees, entete
  `x,y,nw_score,sw_score,cov_xy,glob_xy,cov_yx,glob_yx,violations`),
  produit par `tests/test_align_exhaustive.py` (une seule execution, role
  « Youssef mesure » du plan A9 ; la verification d'OPTIMALITE par enumeration
  est faite separetement dans la mission Sofiene — roles complementaires).- **Methode** : tirage de 500 paires de sequences uniques du corpus v0.2
  (seed 42), application de NW, SW, `identite_couverture` et
  `identite_globale` a chaque paire, puis controle des invariants :
  score recalcule == score annonce (NW et SW), SW >= NW, SW >= 0,
  identites dans [0, 1], `identite_globale <= identite_couverture`,
  symetrie NW(x,y) = NW(y,x) et des deux identites.
- **Resultat** : les 500 lignes portent `violations = OK` (aucun invariant
  enfreint).
- **Limite** : controle d'INVARIANTS sur paires reelles, PAS une preuve
  d'optimalite par enumeration exhaustive (cette derniere verrouille
  `align.py` ailleurs ; ici on verifie la COHERENCE des sorties, pas le
  chemin d'alignement optimal).
- **Empreinte SHA-256** (blob referentiel, LF, attribut `-text`) :
  `ab743664e619cbf32a56fe6ae0567619bcd39cfa6432b3b53f42cca60ab4faa6`


## Evaluation finale des baselines n-gram (A7) — realisation apres merge de la PR #20

- **Date du calcul** : 26/10/2026 (execute dans la branche de mission Sofiene
  apres integration de `origin/main` = `ebe5734`, merge de la PR #20 — les
  splits A6 et le module `baseline_a7` s'y trouvent ensemble pour la
  premiere fois ; tache Cas 1 de la mission « completer PR #20 »).
- **Commande** : `python -c "from peplib.baseline_a7 import run; run()"`
- **Entrees** : `data/splits/{train,val}_v02.csv` (split par FAMILLES) et
  `data/splits_control/{train,val}_v02.csv` (split ALEATOIRE de controle) ;
  colonne `sequence` ; seed 42 ; grille d'alpha imposee `(0.01, 0.1, 1.0)`
  choisie sur val pour chaque ordre.
- **Mesure** : entropie croisee en nats/token sur 21 cibles (20 AA + EOS ;
  BOS/PAD jamais cibles), lissage add-alpha dans `peplib/ngram.py`.
- **Resultats** :
  | split | ordre | alpha | CE (nats/token) | perplexite |
  |---|---|---|---|---|
  | familles | 0 | 0.01 | 3.044522 | 21.000000 |
  | familles | 1 | 1.00 | 2.969770 | 19.487434 |
  | familles | 2 | 1.00 | 2.958478 | 19.268631 |
  | familles | 3 | 1.00 | 2.971710 | 19.525271 |
  | aleatoire | 0 | 0.01 | 3.044522 | 21.000000 |
  | aleatoire | 1 | 0.01 | 2.925784 | 18.648847 |
  | aleatoire | 2 | 0.01 | 2.887591 | 17.950008 |
  | aleatoire | 3 | 1.00 | 2.798959 | 16.427543 |
- **Ancre respectee** : ordre 0 = ln 21 exact sur les deux splits.
- **Lecture** : sur le split par familles, le trigramme NE progresse pas
  (2.9717 > 2.9585 du bigramme) — l'ordre 3 generalise mal quand val ne
  partage aucune famille avec train ; sur le split aleatoire la progression
  est continue (16.43). Les gains restent modestes : ces baselines definissent
  la reference a battre (gate mois-2 : etre meilleur que le unigramme 21
  tokens, PP < 19.49 sur familles).
- **Sorties** : `results/week2/baseline_ngram_A7.csv`,
  `results/week2/ecart_baselines.png`.
- **Empreintes SHA-256** :
  - `baseline_ngram_A7.csv` : `b1f23c74b623f1db5e804d48b9a497c025fe9fe759ff6f85f27f4e602529fe96`
  - `ecart_baselines.png` : `050dbc67fe066eb29b3c753275a2d80d8287bfe80b5a8a0dab9534ca2e5c5b69`


## Variantes du clustering v0.2 (solidite du critere — tableau de revue)

- **Date du calcul** : 26/10/2026 (nocturne) — execute dans la mission Sofiene
  apres merge de la PR #21 ; commandes : `python -c
  "from peplib.variantes_v02 import run; run()"` (candidats A6 recalcules
  UNE fois puis reutilises pour toutes les conventions de gap).
- **Resultats** (`results/week2/variantes_clusters_v02_t13.csv`) :

  | critere | tau | (o,e) | familles | taille_max | singletons | frac_max |
  |---|---|---|---|---|---|---|
  | couverture | 0.90 | 11,1 | 20673 | 4541 | 17592 | 0.106659 |
  | couverture | 0.90 | 10,1 | 20656 | 4544 | 17578 | 0.106729 |
  | couverture | 0.90 | 8,2 | 20596 | 4616 | 17526 | 0.108420 |
  | identite_globale | 0.90 | 11,1 | 26015 | 479 | 22558 | 0.011251 |

- **Lecture** : le critere de couverture 0,90 est **robuste a la convention de
  gap** (ecart max 77 familles sur 20 673 = 0,37 % ; taille max +75) ; la
  mesure identite globale NW (alignement sur la longueur ENTIERE des deux
  sequences) eclate les familles (26 015 familles, taille max 479,
  22 558 singletons) : elle regroupe beaucoup moins, car elle exige la
  ressemblance sur toute la longueur et ne capte pas les fusions par
  fragment — temoin complementaire (mise en garde A5/§5, ancre 20/50).
- **Variante differee** : couverture **0,70 sur le corpus entier** (K(m)
  tombe a 1 pour m = 10 → dizaines de millions de SW ; reporte, le seuil
  0,70 est mesure exactement par la courbe de fuite).

## Diffusion masquée · Contrôles de données S4

Contrôle local du 28 novembre 2026, horodatage de vérification des splits
`2026-11-28T10:23:16Z`. Base Git `b03e74b1f34c5caa1a66b6b553a4cbf1ac8d5398`.
Le code et le tableau mesuré sont identifiés par le commit
`9fad63517f52dbed1dcc6b130dd862d910642ed0`.
Windows AMD64, CPU AuthenticAMD Family 25 Model 80, Python 3.13.14,
PyTorch 2.14.1+cpu, CUDA indisponible.

| Artefact vérifié | SHA-256 des octets |
|---|---|
| `data/splits/train_v02.csv` | `f9016bc038ed473870e76e4a849962813ef010ac3079c31e227309aa4349e736` |
| `data/splits/val_v02.csv` | `f0a17569bf62be16ab931a1f8b47f8c62baad27f509bb62e816aad59f7dbb3a5` |
| `results/week4/masking_control.csv` | `0e07d7a0f8e1096846f78fed4413b72b064ed4b6fc13acf61400066cc32b5b7a` |
| `results/week3/baselines_metrics.csv` | `1e024ec69c3ce3ac9185b25b4ccfd614cd86b2376f600a95676442e55d02141b` |
| `configs/lstm_base.yaml` (S3) | `281dd7ddf3900eb4facf40babcde3c0f6918e3eb9ab4d7c5c28b53454d543715` |
| `configs/transformer_base.yaml` (S3) | `aea05eeb0a7a8700b8375eeb4ee017a6815a1cdebc5140ef99a46c0b83d3ab86` |

Le contrôle synthétique utilise T = 50, un schedule linéaire, la graine 42 et
100 000 positions par étape. Le CSV a été reproduit au même SHA-256 après
normalisation LF. Les effectifs et empreintes train/val concordent avec le
manifeste v0.2. Le split test est resté fermé.

Configuration d’entraînement S4, SHA-256 canonique du YAML, dates des runs,
matériel d’entraînement, métriques, manifeste et hashes des checkpoints S4
sont non disponibles. Les restrictions des sources décrites ci-dessus
s’appliquent. Les séquences générées et checkpoints restent privés.
Les données publiées ici se limitent aux comptages synthétiques et métadonnées.
Le [rapport S4](Docs/Semaine4_rapport.md) précise les exclusions de l’audit
et les validations expérimentales restantes.
