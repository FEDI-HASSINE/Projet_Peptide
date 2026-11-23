# Semaine 2 — Alignement, fusion UniProt + APD3, clusters et baselines (spécification adaptée 4 étudiants)

**Groupe : 4 étudiants (3 développeurs + 1 chef/coordinateur)**
**Encadrement : A. Ben Abdallah, ISIMS**
**Phase B du plan annuel (mois 2, édition semaine 2 ; tâches T30–T47), version peptides**

---

## 1. Objectif et livrables

La semaine 1 a produit un corpus propre (v0.1) et une bibliothèque testée. La semaine 2 prépare **l'entraînement sans triche** : on enrichit le corpus avec **deux bases AMP réelles** (APD3 et la base `antibacterial.fasta` fournie par l'AMP Challenge 2027), on mesure la ressemblance entre peptides (alignement), on regroupe les peptides en **familles** (clusters à 90 %), on découpe le corpus en train/validation/test **par familles**, on traduit les peptides en chiffres (tokenisation), et on installe des **baselines** à battre.

### Qu'est-ce qu'une famille (cluster) ?
Un cluster est un groupe de peptides dont les membres se ressemblent **à 90 % d'identité ou plus**, directement ou par chaînage (A ressemble à B, B ressemble à C ⇒ A, B, C dans le même groupe). Réalisable avec union-find, sans outil externe (BLAST/MMseqs2/cd-hit interdits).

### Qu'est-ce qu'un split étanche ?
Un split étanche est un découpage train/validation/test tel qu'**aucune famille n'apparaisse dans deux jeux à la fois**. Si le même peptide (ou son quasi-jumeau) est présent à l'entraînement et au test, le score du modèle est gonflé sans qu'il ait rien appris : c'est une fuite de données.

### Qu'est-ce que la version v0.2 du corpus ?
Le corpus v0.1 (13 082 peptides, sous-ensembles UniProtKB) enrichi d'APD3 **et de la base AMP Challenge**, fusionné **sans compter deux fois le même peptide** (colonne `sources`), re-filtré, re-clustérisé, avec la provenance, les licences et les empreintes SHA-256 mises à jour dans `PROVENANCE.md` et un `corpus_manifest.json`.

**Résumé en une phrase :** On ajoute deux bases AMP réelles → on mesure la ressemblance entre peptides → on groupe en familles → on découpe sans fuite → on traduit en chiffres → on installe trois baselines à battre.

**Livrables :**
- `peplib/align.py` (Needleman–Wunsch, Smith–Waterman, Gotoh, BLOSUM62) avec tests
- `peplib/cluster.py` (union-find, clusters à 90 %) + `clusters_90.csv` et distribution des tailles
- `peplib/split.py` (90/5/5 stratifié par source, par familles) + preuve d'étanchéité chiffrée
- `peplib/tokens.py` (V = 23, round-trip encode/decode) et `peplib/baselines.py` (aléatoire, unigramme, bigramme)
- `peplib/download_apd3.py` + `peplib/download_challenge.py` + corpus fusionné `corpus_v0.2.csv` (UniProt + APD3 + AMPChallenge) + rapport de recouvrement **trilatéral** (figure)
- `corpus_manifest.json` + `PROVENANCE.md` mis à jour + figures relancées (JS, descripteurs, PWM)
- Tag Git `v0.2`
- Demi-page T27 finalisée si elle ne l'est pas encore

**Critère de réussite :** les quatre étudiants relancent tout depuis `data/raw/` et obtiennent les mêmes comptages, **les mêmes familles** et les mêmes baselines (perplexités identiques).

---

## 2. Répartition des rôles (4 étudiants)

Chaque test et chaque fonction est relu par un autre étudiant que son auteur. Les rôles servent à répartir l'écriture, pas à isoler les connaissances : chacun doit pouvoir expliquer tout le pipeline.

| Rôle | Responsabilités | Tâches principales |
|------|----------------|--------------------|
| **Yassine Achour — Données** | Licences et téléchargements (APD3 + AMPChallenge), connecteurs, résidus modifiés, doublons inter-sources, corpus v0.2, recouvrement, `PROVENANCE.md` | T32, T34, T37, T38 |
| **Sofiene Ouanes — Modèles** | Union-find et clusters, split stratifié et étanche, dédoublonnage par clusters, contrôles qualité | T35, T39, T42, T44, T45 |
| **Youssef Jarraya* — Évaluation** | BLOSUM62, alignements NW/SW/Gotoh, tokenisation, baselines | T33, T36, T40, T43 |
| **Moi — Chef / Coordinateur (L)** | Clôture semaine 1 (tag `v0.1`, relance), manifeste de version, relance des analyses sur v0.2, relectures croisées, séminaire, fusion des branches, revue de tous les PR, merge vers `main` uniquement | T30, T31, T41, T46, T47 + revue de tous les PR |

---

## 3. Structure du dépôt et du module

```
data/raw/                # TOUJOURS immuable : uniprot_*.tsv, apd3_*,
                         # antibacterial.fasta (AMP Challenge 2027), SHA-256
data/processed/          # corpus_v0.2.csv, clusters_90.csv, split_*.csv,
                         # composition/js/redondance relancés, toxin_flags.json
peplib/constants.py      # + BLOSUM62 (saisie depuis la littérature, source citée)
peplib/align.py          # Needleman–Wunsch, Smith–Waterman, pénalités affines (Gotoh)
peplib/cluster.py        # union-find (compression de chemin + union par rang), clusters 90 %
peplib/split.py          # 90/5/5 stratifié par source d'ingestion, par familles
peplib/tokens.py         # Vocabulaire 23, encode/decode round-trip
peplib/baselines.py      # aléatoire uniforme, unigramme 21 tokens, bigramme 21×21
peplib/download_apd3.py  # point d'entrée APD3 (comme download_t8_t9.py)
peplib/download_challenge.py  # point d'entrée base AMP Challenge 2027
tests/                   # un fichier de tests par module
results/                 # figures relancées + recouvrement + distribution des familles
corpus_manifest.json     # version, comptages, hashes — régénérable
```

---

### Workflow Git et collaboration

**1. Chaque personne crée sa propre branche** à partir de `main` : `feat/yassine-apd3`, `feat/sofiene-clusters`, `feat/youssef-align-tokens`, etc.

**2. Chaque personne ouvre une PR** vers `main`. Relue par un autre étudiant que l'auteur.

**3. Toi seul merges vers `main`**, pytest vert (Python 3.11 **et** 3.14).

**Pourquoi ?** Mêmes raisons que la semaine 1 + nouveauté : **une seule nouveauté à la fois** — clusters testés avant fusion, fusion avant baselines. Une dépendance qui part de suite est plus dure à trouver qu'une dépendance intégrée pas à pas.

---

### Roadmap séquentielle des tâches T30–T47

1. **Tous — T30** : clôture S1 : demi-page T27, tag `v0.1`, relance complète
   → puis **Moi — T31** : vérification de l'inventaire semaines 1
   → puis **Yassine — T32** : licences APD3 et AMPChallenge **GO / NO-GO**, téléchargements bruts, PROVENANCE
   → puis **Youssef — T33** : BLOSUM62 saisie + tests
   → puis **Yassine — T34** : connecteurs APD3 + AMPChallenge, résidus modifiés (APD3), parse des en-têtes
   → puis **Sofiene — T35** : union-find, clusters à 90 %
   → puis **Youssef — T36** : Needleman–Wunsch (Gotoh) puis Smith–Waterman
   → puis **Yassine — T37** : doublons exacts inter-sources (3 sources), `corpus_v0.2.csv`
   → puis **Yassine — T38** : rapport de recouvrement trilatéral UniProt ↔ APD3 ↔ AMPChallenge
   → puis **Sofiene — T39** : split stratifié 90/5/5 par familles + preuve d'étanchéité
   → puis **Youssef — T40** : tokenisation V = 23 + round-trip
   → puis **Moi — T41** : `corpus_manifest.json` + test anti-drift
   → puis **Sofiene — T42** : dédoublonnage quasi-doublons inter-sources via clusters
   → puis **Youssef — T43** : baselines aléatoire, unigramme, bigramme
   → puis **Sofiene — T44** (si temps) : validation croisée des descripteurs vs APD3
   → puis **Sofiene — T45** (si temps) : témoins négatifs (séquences mélangées)
   → puis **Tous — T46** : relance complète, tag `v0.2` ; **Moi** merge vers `main`
   → puis **Tous — T47** : séminaire interne (30 min) et présentation à M. Abdallah.

---

## 4. Liste des tâches (T30–T47) avec détails techniques

### A. Clôture de la semaine 1 et nouvelle base

**T30 — Clôturer la semaine 1 (demi-page, tag, relance)**
Rédiger la demi-page T27 (1 paragraphe par étudiant, 6 questions) si elle ne l'est pas encore. Poser le tag `v0.1`. Relance complète : **régénérer `data/processed/` depuis les fichiers bruts existants**, PAS de re-téléchargement (un nouveau download créerait des fichiers datés `_YYYYMMDD` différents et les comptages dériveraient).

**À quoi ça sert en mots simples :** On ferme proprement la semaine précédente avant d'ouvrir la suivante. Le tag `v0.1` est le point de repère : si la v0.2 pose problème, on sait exactement ce qui a changé.

**Validation :** `pytest` vert des deux interprètes ; comptages v0.1 identiques chez les 4 ; tag poussé sur GitHub.

**T32 — Licences APD3 et AMPChallenge : GO / NO-GO, téléchargements et provenance**
Avant toute manipulation, deux **décisions écrites** :

*Source 1 — APD3* : lire les conditions d'usage et de **redistribution** (site de l'Université du Nebraska — `aps.unmc.edu`).
- **GO** : utilisation et redistribution libres ou académiques → traiter APD3 comme les autres sources ;
- **usage libre mais redistribution interdite** → APD3 utilisable localement, fichier brut **non commité** (noté dans `PROVENANCE.md` et `.gitignore`) ;
- **NO-GO** → basculer sur DBAASP ou DRAMP (le schéma du connecteur ne change pas) ; **filet doré** : ces deux bases figurent déjà dans la base AMP Challenge (cf. ci-dessous), donc le corpus AMP reste assuré même en cas de NO-GO APD3.

*Source 2 — AMP Challenge 2027* (`data/antibacterial.fasta` du dépôt `szczurek-lab/amp-challenge-2027`) :
- téléchargement direct : `https://raw.githubusercontent.com/szczurek-lab/amp-challenge-2027/main/data/antibacterial.fasta` ;
- licence du dépôt : **BSD-3-Clause** (permissive, favorable a priori) ; citer le dépôt et le challenge parleur ;
- **attention** : le fichier agrège **13 bases amont** (dbAMP, DRAMP, DBAASP, CAMP, SATPdb, APD, AMPDB, DADP, InverPep, CancerPPD, BaAMPs, CyBase, ParaPep — visible dans le champ `dbs=` des en-têtes). La licence BSD-3 du dépôt ne remplace pas les licences amont : avant tout **partage** du corpus v0.2, vérifier chaque licence amont et le consigner dans `PROVENANCE.md` ;
- **stats vérifiées lors de l'analyse du fichier** (à re-vérifier au téléchargement, le dépôt peut le mettre à jour) : 39 448 séquences, 8–50 résidus (médiane 18), **34 474 dans notre fenêtre 10–50**, 0 doublon exact, alphabet 100 % canonique, en-têtes `len, charge, disulfide, activity, dbs` ;
- empriente SHA-256 du fichier brut : `cbbeac64ba95746d87961e8ad9dd0849ae8058d15a300b2e7f6990730ca521e9`.

Télécharger les fichiers bruts dans `data/raw/`, calculer les SHA-256, consigner dans `PROVENANCE.md` (nombre de lignes, date, version, licence, téléchargement manuel Oui/Non) — **un bloc par source**.

**À quoi ça sert en mots simples :** C'est amener les bases AMP « réelles » (avec cibles et activités pour APD3, avec métadonnées en en-tête pour le challenge) dans notre dépôt, en respectant leurs règles. On ne regarde pas les licences après : on les lit **avant**, sinon tout le travail en aval peut être à refaire.

**Validation :** entrées `PROVENANCE.md` complètes (× 2) ; hashes recalculables à la main ; décisions de redistribution consignées.

**T33 — BLOSUM62 : saisie depuis la littérature**
Matrice 20 × 20 des scores de substitution, source citée (Henikoff & Henikoff, 1992). Convention : `SCORE(a, b)` entier, `SCORE(a, b) = SCORE(b, a)`. Stocker dans `peplib/constants.py` (ou `peplib/align.py`) sous forme de dict de dict ou matrice NumPy.

**À quoi ça sert en mots simples :** C'est la « table des prix » de l'alignement : remplacer L par I ne coûte presque rien, remplacer A par W coûte cher. Sans table justifiée, tout alignement est au hasard.

**Validation :** 400 valeurs saisies, symétrique, diagonale positive ; valeurs d'ancres rapprochées de la source citée (ex. A/A = 4, W/W = 11) ; **relire à quatre yeux** : une valeur fausse fausse tous les alignements.

### B. Connecteur APD3 et fusion

**T34 — Connecteurs APD3 et AMPChallenge, résidus modifiés (dépend de T32, T35 pour la fusion)**
Standardiser les deux sources vers le schéma commun : `sequence, source, id` (+ métadonnées persistées pour l'aval mois 3 : APD3 → `nom, organisme, cible, activité, unité` ; AMPChallenge → `charge, disulfide, activity, dbs` extraits des en-têtes). Appliquer **exactement** les mêmes filtres que la semaine 1 : `dropna()` → `strip()` → `upper()` → regex → 10–50.

**Connecteur AMPChallenge (le plus simple des deux)** : parser les en-têtes FASTA avec une regex documentée (`len=`, `charge=`, `disulfide=`, `activity=`, `dbs=`). Le fichier est déjà propre : **0 doublon exact, alphabet 100 % canonique, aucun résidu modifié**. La provenance multi-bases du champ `dbs=` est conservée en colonne `dbs_amont` (13 bases possibles) mais **ne pilote aucun calcul** cette semaine (cf. décision du split, T39).

**Résidus modifiés — décision écrite obligatoire (APD3 uniquement)** (APD3 contient des peptides amidés en C-terminal, des D-acides aminés, et des lettres non canoniques `B X Z J U O`) :
1. lettres ambigües/rares → **rejetées** comme en semaine 1 (première cause, logique inchangée). `classify_sequence` les marque déjà : **ne pas réinventer** ;
2. notation de l'amidation (souvent un symbole ou un signe en fin de séquence) → **annotée dans une colonne dédiée** puis retirée de la séquence brute avant filtrage ; **documenter la convention de parse** avec exemples réels ;
3. tout élément non couvert par ces deux règles → à lister, puis décision du coordinateur, notée dans `PROVENANCE.md`.

**À quoi ça sert en mots simples :** APD3 contient des peptides « modifiés » (sucres additionnels, extrémité bloquée). On explique ce qu'on en fait **par écrit** : ce qu'on rejette, ce qu'on garde, pour que le corpus reste propre et que personne ne découvre la pollution en semaine 4. Le connecteur du challenge, lui, est « tout-terrain » : les en-têtes donnent les métadonnées sans casse possible.

**Validation :** `classify_sequence` régénère les mêmes causes de rejets qu'en semaine 1 ; pour le challenge : 39 448 lignes parsées, champs d'en-tête renseignés à 100 % ; colonnes de métadonnées remplies pour ≥ 90 % des lignes APD3 conservées.

**T37 — Doublons exacts inter-sources et corpus v0.2 (dépend de T34)**
Intégrer APD3 et AMPChallenge au corpus. Pour chaque séquence existant dans plusieurs sources : **une seule ligne**, colonne `sources` = liste des atomes `{UniProt, APD3, AMPChallenge}` (ex. `"UniProt+APD3"`, `"UniProt+APD3+AMPChallenge"` — jusqu'à 7 combinaisons possible). Règle de fusion : on garde les **métadonnées les plus complètes** (APD3 : cible + activité ; AMPChallenge : charge/dbs) et l'accession UniProt en référence croisée (colonne ou champ `ids`).

**À quoi ça sert en mots simples :** Beaucoup de peptides APD3 sont déjà dans UniProt. Sans fusion intelligente, on compterait deux fois le même peptide : le corpus paraîtrait plus grand qu'il n'est, et toutes les statistiques seraient gonflées.

**Validation :** `conservés + rejetés = total brut` maintenu toutes sources confondues ; aucune `sequence` dupliquée dans `corpus_v0.2.csv` ; extraction par source toujours possible (vues `uniprot.csv`, `amp_apd3.csv`, `challenge.csv`, `corpus_v0.2.csv` — quatre vues issues du même pipeline).

**T38 — Rapport de recouvrement trilatéral (dépend de T37)**
Chiffres : pour chaque paire UniProt↔APD3, UniProt↔AMPChallenge, APD3↔AMPChallenge + le triple recouvrement : nb d'exactes communs, nb en quasi-doublons (via clusters, T42), nb uniques par source. Figure (barres empilées, matplotlib 2D).

**À quoi ça sert en mots simples :** C'est la preuve chiffrée que la fusion ne triche pas : on montre combien de peptides ont été « identifiés comme les mêmes » et on ne les cache pas. Avec trois sources, on voit aussi combien de peptides sont présents dans les trois à la fois.

**Validation :** la somme des parts = taille du corpus fusionné ; figure légendée, sauvegardée dans `results/`.

**T41 — Manifeste de version du corpus (dépend de T37)**
`corpus_manifest.json` : version (`v0.2`), date, total, **compte par source d'ingestion (3 entrées)**, nb de familles (T35), hashes SHA-256 des CSV majeurs **et des deux fichiers bruts AMP**, seed. Test : si un fichier du corpus change sans bump de version, le test échoue (hash dans manifeste ≠ hash fichier).

**À quoi ça sert en mots simples :** C'est le « sceau » du corpus : n'importe qui peut vérifier que ce qu'il a régénéré est bien la même version. Sans manifeste, deux versions cohabitent sans qu'on le sache.

**Validation :** régénération complète → manifeste identique octet pour octet.

### C. Alignement

**T36 — Needleman–Wunsch, Gotoh et Smith–Waterman (dépend de T33)**
- **Needleman–Wunsch (global)** avec **pénalités affines de Gotoh** : passages E (horizontal) et F (vertical) séparés de la diagonale H. **Convention documentée** : un gap de longueur `L` coûte `G_open + (L − 1) × G_extend`. Défauts proposés : `G_open = 11`, `G_extend = 1` (à justifier !). Retourner score + alignement (traceback).
- **Smith–Waterman (local)** : même moteur, avec `max(0, …)` et retour du meilleur score local + l'alignement correspondant.
- Identité via alignement : `identité = nb positions alignées identiques / longueur de l'alignement` — complémentaire de l'identité Levenshtein (T23) ; documenter la différence.

**À quoi ça sert en mots simples :** L'alignement compare deux peptides lettre à lettre en minimisant un coût : d'où viennent les substitutions, où il y a des trous. Le global aligne de bout en bout, le local ne garde que la portion qui colle. Gotoh fait que trois lettres d'un coup ne coûtent pas trois fois le prix d'un seul trou (un trou coûte l'ouverture, puis un petit supplément par allongement).

**Validation (ancres à la main, convention figée ci-dessus) :**
- `NW("AAAA", "AAAA") = 16` (4 × A/A = 4×4) ; `NW("GGGG", "GGGG") = 24` (4 × G/G = 4×6)
- `NW("AAAA", "AAAAA") = 16 − G_open = 5` (un trou de longueur 1)
- `NW("AAAA", "AAAAAAA") = 16 − (G_open + 2×G_extend) = 3` (un trou de longueur 3)
- `SW("AAAA", "AAAA") = 16` ; `SW ≥ 0` toujours ; le résultat de SW est une **sous-chaîne contiguë**.
- Si les ancres numériques ne correspondent pas, c'est la **convention** qui a bougé : corriger la documentation, pas la valeur du test.

### D. Clusters et splits

**T35 — Union-find et clusters à 90 % (dépend de T23, T33 pour l'identité alignée)**
Union-find avec compression de chemin et union par rang. Paires candidates depuis l'index de 3-mers existant (T23) ; mesure d'identité exacte sur candidats (Levenshtein normalisée, et alignement pour vérification sur échantillon) ; union des paires ≥ 90 %. **Chaînage transitif voulu** : si A~B et B~C mais A–C < 90 %, les trois sont dans le même groupe (c'est le but d'union-find : la transitivité garantit l'étanchéité du split). Livrables : `clusters_90.csv` (`sequence_id, cluster_id, source`), nb de familles, distribution des tailles (figure).

**À quoi ça sert en mots simples :** On transforme 13 000 peptides isolés en quelques centaines de **familles**. C'est l'unité anti-triche du split, et accessoirement la version « propre » du taux de quasi-doublons de la semaine 1.

**Validation :** transitive : {A,B} groupes puis {B,C} ⇒ même groupe ; complexité quasi-linéaire (compression de chemin) ; nb de familles cohérent avec les 5 675 paires de quasi-doublons de T25 (ordre de grandeur).

**T39 — Split stratifié 90/5/5 par familles + preuve d'étanchéité (dépend de T35)**
Découpage **par familles**, **stratifié par source d'ingestion** (au sein de chaque source : 90/5/5 des familles, tirage seed 42) pour garantir que chaque jeu contient des UniProt, des APD3 et des AMPChallenge. **Décision** : la stratification se fait par source d'ingestion (`UniProt` / `APD3` / `AMPChallenge`) — pas par l'une des 13 bases amont agrégées dans le champ `dbs=` du challenge, qui reste une colonne d'information. Écriture de `split_train.csv`, `split_val.csv`, `split_test.csv`.

**Rapport d'étanchéité (chiffres à montrer) :**
1. aucune famille dans deux jeux ;
2. proportions 90/5/5 respectées par source (± tolérance documentée) ;
3. **identité max entre un élément de train et un de test < 90 %** (calculée sur candidats 3-mers) — c'est LE chiffre de la semaine ;
4. taille des jeux et nb de familles par jeu.

**À quoi ça sert en mots simples :** C'est la règle « une famille ne se partage pas ». Sans stratification, on peut avoir tous les AMP dans le train et zéro dans le test : le split serait étanche mais inutilisable. Avec les deux, on a un découpage honnête et exploitable. L'ordre des étapes : **clusters fusionnés d'abord (une seule fois pour tout le corpus à 3 sources)**, ensuite seulement le split.

**Validation :** union des 3 = corpus ; intersections vides ; proportions 90/5/5 respectées par source ; `split` régénérable (seed) ; test anti-fuite dans `pytest`.

**T42 — Dédoublonnage quasi-doublons inter-sources via clusters (dépend de T35, T37)**
Au sein de chaque famille : choisir **un représentant** (règle : la ligne avec métadonnées les plus riches, préférer APD3 si cible/activité disponible). Garder `membres` ou un fichier de mapping `cluster_id → ids/sources`. Trois modes du corpus : brut (tous), fusionné (v0.2), dédupliqué (1 représentant/famille) — les 3 sont régénérables et documentés.

**À quoi ça sert en mots simples :** Le cluster sert deux fois : il protège le split **et** il dédoublonne entre bases. Un seul mécanisme, deux bénéfices — c'est pour ça qu'on attend la fin des clusters avant la fusion propre (T37).

**Validation :** aucune famille avec 2 représentants ; mapping complet (chaque peptide a un cluster) ; opération en temps quasi linéaire.

### E. Tokenisation et baselines

**T40 — Tokenisation, V = 23 (dépend de T37 pour le corpus)**
Vocabulaire : `0 = PAD`, `1 = BOS`, `2 = EOS`, puis les 20 acides aminés **dans l'ordre alphabétique** `ACDEFGHIKLMNPQRSTVWY` → ids 3..22. Encode : `BOS + ids + EOS` ; decode inverse. Round-trip sur 1 000 séquences: decode(encode(s)) == s pour 100 %.

**À quoi ça sert en mots simples :** Le réseau ne lit pas « ACDE » : il lit des chiffres. On traduit, et on vérifie que retraduire redonne le peptide d'origine, à 100 %. Une erreur ici contamine tout le mois 2.

**Validation :** `len(VOCAB) == 23` ; round-trip 100 % ; ids stables (pas d'ordre qui dérive) ; lettres hors alphabet impossibles (elles sont rejetées en amont).

**T43 — Baselines : aléatoire, unigramme, bigramme (dépend de T40, T39)**
Sur le split train uniquement :
- **Aléatoire uniforme** : tirage uniforme. Théorie : `−log p ≈ ln 20` nats/lettre sur les lettres.
- **Unigramme (21 tokens : 20 AA + EOS)** : `p(a) = compte(a) / total`, lissage α = 1 (pseudo-count justifié). Perplexité = exp(−(1/N) Σ −log p(x_i)) en nats/lettre.
- **Bigramme (21 × 21)** : `p(b|a)`, lissage identique ; init à BOS (token), conditionne et finit par EOS.

**À quoi ça sert en mots simples :** Trois joueurs d'entraînement au niveau croissant : le premier tape au hasard, le deuxième connaît la fréquence des lettres, le troisième connaît les paires de lettres qui collent bien. En semaine 3, le LSTM devra les battre dans cet ordre.

**Validation :** chaque ligne du tableau bigramme somme à 1 ; perplexités des 3 baselines rapportées **avec seed** ; bigramme ≤ unigramme ≤ aléatoire (ordre théorique attendu sur le même split).

### F. Contrôles qualité et clôture

**T44 — Validation croisée des descripteurs (si temps)**
Deux contrôles externes gratuits :
- **APD3** publie charge, pI, masse pour ses peptides → comparer nos valeurs : médiane et écart maximal des écarts ;
- **AMPChallenge** publie la **charge** dans les en-têtes de 39 448 séquences → comparer notre Henderson–Hasselbalch (pH 7) sur les lignes challenge du corpus fusionné.

Un écart systématique = signe d'un bug de constante (convention de protonation, masse moyenne vs monoisotopique) → à clarifier dans la doc, pas à masquer.

**À quoi ça sert en mots simples :** C'est un contrôle de qualité **externe** : les chiffres d'APD3 n'ont pas été calculés par notre code. Si les deux concordent, nos constantes sont bonnes ; sinon on a trouvé une erreur **avant** d'entraîner un modèle dessus.

**Validation :** rapport chiffré (médiane et écart maximal) ; discrepancies documentées et justifiées (convention de calcul), pas ignorées.

**T45 — Témoins négatifs : séquences mélangées (si temps)**
Prendre 500 AMP du corpus, **mélanger les lettres** à seed fixée, recalculer les descripteurs. Attendu : charge et moment hydrophobe distinguent AMP réels de mélangés (le moment hydrophobe dépend de l'**ordre**, pas seulement de la composition, il doit chuter après mélange ; la charge reste identique — c'est justement la preuve que chaque descripteur mesure une chose différente).

**À quoi ça sert en mots simples :** Si nos descripteurs ne voient pas la différence entre un peptide réel et son omelette mélangée, ils ne servent à rien. Deux figures suffisent : avant / après mélange.

**Validation :** 2 figures + comparaison chiffrée (ex. médiane moment hydrophobe AMP vs mélangés) ; seed écrite.

**T46 — Relance complète et tag `v0.2`** — comme T30 mais pour v0.2 : régénération depuis `data/raw/`, mêmes comptages (détail intégral dans le manifeste), **mêmes familles**, mêmes perplexités ; tag `v0.2`.

**T47 — Séminaire interne et préparation mardi** — ordre de passage : recouvrement et v0.2 (Yassine), familles et étanchéité (Sofiene), tokenisation et baselines (Youssef), synthèse coordonnée (Moi).

---

## 5. Formules à implémenter

| Quantité | Définition exacte |
|----------|-------------------|
| Score d'alignement | `score = Σ s(x_i, y_i) − Σ pénalités de gap`, `s` = entrée BLOSUM62 |
| Gap affine (Gotoh) | un gap de longueur `L` coûte `G_open + (L − 1) × G_extend` ; définir la convention en tête de module |
| NW global | récurrences H (diag), E (horizontal), F (vertical) ; `H(i,j) = max(H(i−1,j−1)+s, E(i,j), F(i,j))` ; `E(i,j) = max(H(i,j−1) − G_open, E(i,j−1) − G_extend)` ; F symétrique vertical |
| SW local | mêmes récurrences avec `max(0, …)` ; meilleur score local + traceback |
| Identité (alignement) | `identité = NB positions identiques / longueur d'alignement` |
| Identité (Levenshtein) | `1 − d_Lev(a,b) / max(|a|,|b|)` — la mesure du clustering 90 % |
| Union-find | `find` avec compression de chemin ; `union` par rang ; transitif |
| Tokenisation | `V = 23` : `PAD=0, BOS=1, EOS=2, AA 3..22` (ordre alphabétique) |
| Unigramme | `p(a) = (compte(a) + α) / (total + α·21)` avec `α = 1` (lissage additif documenté) |
| Bigramme | `p(b|a) = (compte(a,b) + α) / (compte(a) + α·21)` |
| Perplexité (nats/lettre) | `PP = exp( −(1/N) Σ −log p(x_i) )` — plus bas = mieux |
| Split | 90/5/5 des **familles**, stratifié par source, seed 42 |

---

## 6. Tests obligatoires

| Test | Résultat attendu |
|------|------------------|
| BLOSUM62 | 400 valeurs ; symétrique ; diagonale > 0 ; ancres rapprochées de la source citée (ex. A/A = 4, W/W = 11 ; **vérifier contre la source**) |
| Connecteurs | APD3 : causes de rejet régénérées. AMPChallenge : 39 448 séquences parsées ; en-têtes `len/charge/disulfide/activity/dbs` renseignés à 100 % ; 0 doublon exact ; alphabet 100 % canonique |
| NW identités | `NW("AAAA","AAAA") = 16` ; `NW("GGGG","GGGG") = 24` |
| NW affine | `NW("AAAA","AAAAA") = 5` ; `NW("AAAA","AAAAAAA") = 3` avec la convention `G_open=11, G_extend=1` |
| SW | `SW("AAAA","AAAA") = 16` ; `SW ≥ 0` ; résultat = sous-chaîne contiguë |
| Identité Levenshtein | `("kitten","sitting")` ⇒ 3 et 1 − 3/7 ≈ 0,5714 |
| Union-find | {A,B} puis {B,C} ⇒ même groupe (transitivité) ; find après union idempotent |
| Clusters | aucun membre dans 2 familles ; nb familles cohérent avec T25 |
| Round-trip | `decode(encode(s)) == s` sur 1 000 séquences ; `len(VOCAB) == 23` |
| Bigramme | chaque ligne conditionnelle somme à 1 à 10⁻⁹ près |
| Ordre baselines | perplexité bigramme ≤ unigramme ≤ aléatoire (sur le même split train) |
| Split étanche | union = corpus ; intersections vides ; aucune famille partagée ; proportions 90/5/5 par source d'ingestion (3) ± tolérance documentée |
| Split sans triche | **identité max train↔test < 0,90** (rapport chiffré) |
| Corpus v0.2 | conservés + rejetés = brut toutes sources ; aucune `sequence` dupliquée ; colonne `sources` ∈ atomes {UniProt, APD3, AMPChallenge} (combinaisons `+` documentées) ; vues enregistrables par source |
| Manifeste | régénération → manifeste identique (3 compteurs de source + 2 hashes AMP, seed) ; filet : test échoue si hash ≠ sans bump de version |

---

## 7. Pièges à éviter

- **Convention de gap** : `open + (L−1)×extend` vs `open + L×extend` — deux écoles. Documenter la nôtre en tête de module ; toute divergence d'ancre vient de là, pas du code.
- **BLOSUM62** : une valeur fausse fausse tous les alignements, donc tous les clusters. Les ancres « connues » doivent être re-vérifiées **dans la source citée**, pas de mémoire.
- **Transitivité des clusters** : c'est voulu (anti-fuite). Ne pas la « corriger » en retirant les chaînages longs.
- **Split stratifié** : ne pas stratifier **avant** d'avoir les familles ; ne pas mélanger famille et peptide comme unités de découpage.
- **Relance ≠ re-téléchargement** : `data/raw/` est immuable ; régénérer en aval, jamais télécharger de nouveau le corpus d'UniProt (fichiers datés → comptages qui dérivent).
- **PAD/BOS et EOS** : les ids envoient BOS et PAD avec des logits à −∞ quand on échantillonnera en semaine 3 ; dès la tokenisation, ne jamais laisser un modèle (même la baseline) tirer PAD ou BOS.
- **APD3 résidus modifiés** : amidation, D-acides aminés, lettres ambigües — règles de rejet **identiques** à la semaine 1 (première cause) ; toute exception doit être documentée dans `PROVENANCE.md`, pas codée en dur au fond d'une fonction.
- **Doublons inter-sources** : sans colonne `sources`, on compte deux fois le même peptide ; sans clusters, on ne voit pas les quasi-doublons entre bases. Avec trois sources, les 7 combinaisons de la colonne `sources` doivent être testées, pas seulement les paires.
- **Corpus plus gros (~45–50k séquences)** : l'index 3-mers tient, mais la **mesure d'identité exacte sur candidats** peut devenir longue → échantillonnage documenté pour la vérification d'alignement (le cluster lui-même reste exhaustif).
- **Parse des en-têtes FASTA du challenge** : regex figée et testée (`len=[0-9]+`, `charge=...`, ...) ; ne jamais dériver les métadonnées de l'ordre des en-têtes ni d'une virguleRelativeTo des espaces : chaque champ est parsé indépendamment.
- **Licences amont du challenge** : le BSD-3 couvre le dépôt (code), pas forcément chacune des 13 bases agrégées dans `dbs=` → avant tout partage, vérification ligne à ligne dans `PROVENANCE.md`.
- **Seeds** : 42 partout ; le split, les témoins négatifs, les baselines et le manifeste doivent être reproductibles.
- **Notebook contre module** : la logique vit dans `peplib` et ses tests ; le notebook appelle seulement (comme NB01).
- **Visualisation 3D interdite** (rappel) : matplotlib 2D uniquement.
- **UTF-8** : fichiers sans BOM, sans mojibake (vérifier chaque écriture comme en semaine 1).

---

## 8. Biosécurité, provenance et licences

- **APD3** (base du Nebraska, `aps.unmc.edu`) peut contenir des peptides hémolytiques ou toxiques ; le flag toxines de la semaine 1 (`toxin_flags.json`) est **relancé sur le corpus fusionné** ; jamais de séquences toxin-like dans le corpus génératif, jamais distribuées.
- **AMPChallenge** : le fichier est un pré-filtré « peptides antibactériens » sans résidus exotiques ; le flag toxines s'applique aussi. Le BSD-3 du dépôt ne remplace pas les licences des 13 bases amont agrégées (`dbs=`) : consignées dans `PROVENANCE.md` **avant** tout partage du corpus v0.2.
- **Redistribution** : la décision T32 consigne par écrit ce qui est autorisé et ce qui ne l'est pas ; si la redistribution d'APD3 est interdite, le brut reste local, non commité, et les produits dérivés à redistribuer ne l'incluent pas.
- Tout partage du corpus v0.2 passe par la mise à jour de `PROVENANCE.md` **avant** le partage.

---

## 9. Critères de réussite de la semaine

- Les tâches T30–T47 sont cochées, ou les restantes (T44, T45) reportées avec une raison écrite.
- Corpus v0.2 : **trois sources ingérées**, invariants de taille, de traçabilité (`sources`) et de filtre tenus ; aucune séquence dupliquée.
- Clusters : nb de familles rapporté ; distribution des tailles ; test de transitivité vert.
- Splits : étanches (test anti-fuite), stratifiés par source, **identité max train/test < 90 %** (chiffre rapporté).
- Tokenisation : round-trip 100 % sur 1 000 séquences.
- Baselines : 3 perplexités rapportées avec seeds ; ordre théorique respecté.
- Manifeste et `PROVENANCE.md` à jour ; les 4 étudiants régénèrent v0.2 avec les mêmes résultats ; tag `v0.2` poussé.

---

## 10. Ce que le superchef (moi) doit présenter à M. Abdallah

M. Abdallah doit voir des **visuels et documents**, pas du code — et cette semaine, des **chiffres de garde-fou**.

### Visuels à présenter (nouveaux cette semaine)
- **Rapport de recouvrement trilatéral UniProt ↔ APD3 ↔ AMPChallenge** (barres empilées : exacts / quasi-doublons / uniques, 3 paires + triple)
- **Distribution des tailles de familles** (combien de familles, comment réparties)
- **Preuve d'étanchéité des splits** (tableau : proportions par source, identité max train/test)
- Figures de la semaine 1 **relancées sur v0.2** (composition, JS, descripteurs, PWM) — montrer que « v0.2 régénère tout »
- **Échelle de baselines** (un tableau projectif : aléatoire / unigramme / bigramme, perplexités en nats/lettre)

### Tableaux et documents à présenter
- `PROVENANCE.md` mis à jour (licences APD3 + AMPChallenge, 13 bases amont, hashes, version v0.2)
- `corpus_manifest.json` (la « carte d'identité » de v0.2)
- Rapport des résidus modifiés d'APD3 (décision rejeter/annoter/convertir, chiffrée)
- Demi-page T27 finalisée

### Preuve de reproductibilité à montrer
- Tags `v0.1` **et** `v0.2` poussés
- Les 4 étudiants régénèrent v0.2 et obtiennent mêmes comptages, mêmes familles, mêmes perplexités

**On doit savoir expliquer chacun de ces éléments sans ouvrir le code**, surtout : pourquoi le cluster protège l'examen (le doublon fait tricher, pas le thème), pourquoi la stratification par source colle des AMP dans le test, et ce que signifie nats/lettre.

---

## 11. Commandes de démarrage rapide (à exécuter dans l'ordre)

Ouvrir un terminal **à la racine du projet** (nouvel emplacement : `C:\Users\DELL\Projects\Projet_Peptide`) :

```bash
# 1. (une fois) Créer l'environnement virtuel
python -m venv .venv

# 2. L'activer (Windows PowerShell)
.venv\Scripts\Activate.ps1

# 3. Installer les dépendances
pip install -r requirements.txt

# 4. Relancer les tests (suite complète, 2 interprètes)
pytest
python -m pytest

# 5. Régénérer le corpus traité depuis les bruts (NE PAS re-télécharger UniProt)
python -m peplib.filters          # ou l'entrée d'exécution utilisée en S1

# 6. Clusters, split et baselines (les modules de cette semaine)
#    python -m peplib.cluster ; python -m peplib.split ; python -m peplib.baselines
```

Résultat attendu : tous les tests passent, aucun échec ; les comptages et les perplexités sont identiques d'une relance à l'autre.

---

## 12. Template `PROVENANCE.md` (bloc APD3 à remplir)

```markdown
# Source : APD3

- **Nom de la source** : APD3 (Antimicrobial Peptide Database, Univ. du Nebraska — aps.unmc.edu)
- **URL de requête ou de téléchargement** : [page de téléchargement exacte, ex. https://aps.unmc.edu/download]
- **Date de téléchargement** : [JJ/MM/AAAA]
- **Version de la base** : [si indiquée, sinon "non spécifiée"]
- **Conditions d'usage / licence** : [copier la phrase exacte du site + verdict : redistribution autorisée / interdite / usage local uniquement]
- **Nombre de lignes brutes** : [nb de séquences brutes, + format (FASTA ?)]
- **Empreinte SHA-256 du fichier brut** : [hash]
- **Fichier brut sauvegardé sous** : `data/raw/apd3_YYYYMMDD.fasta` (non commité si redistribution interdite → noter ici)
- **Téléchargement manuel ?** : [Oui / Non]
- **Notes** : [résidus modifiés : amidation / D-acides aminés / lettres non canoniques — décision appliquée ; toxines ou hémolysines signalées ; recouvrement avec UniProt]
```

Exemple de rédaction de la ligne « Décision » (à adapter au verdict réel, Ne pas inventer) :

> « Usage académique autorisé sur consultation du site le [date] ; redistribution [autorisée / interdite] → le fichier brut est [commité dans `data/raw/` / ignoré par `.gitignore` et stocké localement]. Toute distribution du corpus v0.2 [inclut / n'inclut pas] les séquences APD3. »

Deuxième bloc à remplir — la base AMP Challenge (pré-rempli avec ce que l'analyse du 22/10 a déjà vérifié) :

```markdown
# Source : AMP Challenge 2027 — antibacterial.fasta

- **Nom de la source** : AMP Challenge 2027, dépôt `szczurek-lab/amp-challenge-2027` (donnée fournie par l'encadrant)
- **URL de téléchargement** : https://raw.githubusercontent.com/szczurek-lab/amp-challenge-2027/main/data/antibacterial.fasta
- **Date de téléchargement** : [JJ/MM/AAAA]
- **Version de la base** : branche `main` du dépôt (fichier susceptible d'être mis à jour → hash fait foi)
- **Conditions d'usage / licence** : dépôt sous **BSD-3-Clause** ; le fichier agrège 13 bases amont (dbAMP, DRAMP, DBAASP, CAMP, SATPdb, APD, AMPDB, DADP, InverPep, CancerPPD, BaAMPs, CyBase, ParaPep) dont les licences propres sont à vérifier **avant tout partage** du corpus
- **Nombre de lignes brutes** : 78 896 lignes ≈ 39 448 séquences (constaté le 22/10/2026 ; refaire le compte à chaque téléchargement)
- **Empreinte SHA-256 du fichier brut** : cbbeac64ba95746d87961e8ad9dd0849ae8058d15a300b2e7f6990730ca521e9
- **Fichier brut sauvegardé sous** : `data/raw/challenge_antibacterial_YYYYMMDD.fasta`
- **Téléchargement manuel ?** : Non (URL directe, sans clé ni inscription)
- **Notes** : séquences 8–50 (médiane 18), 34 474 dans la fenêtre 10–50, 0 doublon exact, alphabet 100 % canonique, en-têtes parsés `len/charge/disulfide/activity/dbs` ; aucun résidu modifié à traiter ; champ `dbs=` conservé en colonne `dbs_amont`, ne pilote aucun calcul en semaine 2
```
