# Semaine 1 — Corpus de peptides et bibliothèque `peplib` (spécification adaptée 4 étudiants)

**Groupe : 4 étudiants (3 développeurs + 1 chef/coordinateur)**  
**Encadrement : A. Ben Abdallah, ISIMS**  
**Phase A du plan annuel (tâches T1–T29), version peptides**

---

## 1. Objectif et livrables

Cette semaine, aucun modèle n'est entraîné. L'objectif est d'obtenir un **corpus propre, documenté et reproductible**, et une petite bibliothèque Python testée, `peplib`, qui servira pendant toute l'année.

### Qu'est-ce qu'un corpus ?
Un corpus est l'ensemble des séquences de peptides collectées (ici via UniProt et d'autres bases), avec leur provenance (source, date, requête, licence, nombre de lignes, empreinte SHA-256).

### Qu'est-ce qu'un corpus propre ?
Un corpus propre est ce même corpus après traitement :
- filtré par regex `^[ACDEFGHIKLMNPQRSTVWY]+$` et longueur 10–50
- `NaN` supprimés avant toute conversion en chaîne
- doublons exacts et quasi-doublons identifiés (index 3-mers + Levenshtein)
- toxines (conotoxines, etc.) signalées ou exclues
- exporté en CSV (`sequence`, `source`, `id`) dans `data/processed/`
- documenté dans `PROVENANCE.md` et reproductible depuis `data/raw/`

**Résumé en une phrase :** On téléporte des peptides du web → on jette les invalides → on enlève les doublons → on calcule leurs propriétés → on sauvegarde un fichier propre et reproductible.

**Livrables :**
- Module `peplib/` (constantes, filtrage, composition, PWM, descripteurs, redondance) avec tests `pytest`
- Notebook `NB01` version peptides (appelle `peplib` ; aucune logique métier copiée dans les cellules)
- `data/` avec fichiers bruts, fichiers traités et `PROVENANCE.md`
- Table des rejets par cause, histogrammes, composition, PWM et information par position
- Demi-page d'analyse rédigée
- Tag Git `v0.1`

**Critère de réussite :** les quatre étudiants relancent tout le pipeline de bout en bout, depuis les fichiers bruts, et obtiennent les mêmes comptages.

---

## 2. Répartition des rôles (4 étudiants)

Chaque test et chaque fonction est relu par un autre étudiant que son auteur. Les rôles servent à répartir l'écriture, pas à isoler les connaissances : chacun doit pouvoir expliquer tout le pipeline au séminaire de fin de semaine.

| Rôle | Responsabilités | Tâches principales |
|------|----------------|--------------------|
| **Yassine Achour — Données** | Acquisition, provenance, filtrage, statistiques de longueur et de composition, redondance exacte | T7–T14, T16, T17, T22, T25 |
| **Sofiene Ouanes — Modèles** | Constantes, descripteurs, PWM et information, quasi-doublons | T4, T18, T20, T21, T23 |
| **Youssef Jarraya* — Évaluation** | Seeds, tests unitaires de chaque module, entropie de composition, relance finale | T2, T3, T5, T15, T19, T24, T29 |
| **Moi — Chef / Coordinateur (L)** | Création du dépôt Git et de l'arborescence, coordination, relectures croisées, intégration, notebook, analyse, séminaire, fusion des branches, revue de tous les PR, merge vers `main` uniquement | T1, T6, T26, T27, T28 + revue de tous les PR |

---

## 3. Structure du dépôt et du module

```
README.md          # ordre d'exécution, versions, seeds, commandes de test
requirements.txt   # versions figées (Python, NumPy, SciPy, pandas, matplotlib, PyTorch, pytest)
data/raw/          # fichiers téléchargés, JAMAIS modifiés ; empreinte SHA-256 consignée
data/processed/    # CSV filtrés : colonnes sequence, source, id
peplib/constants.py   # 20 AA : masse résiduelle monoisotopique, pKa chaîne latérale, Kyte-Doolittle, Eisenberg ; pKa extrémités N et C ; référence citée pour chaque table
peplib/filters.py     # filtrage par alphabet et longueur ; renvoie séquences conservées + dictionnaire de rejets par cause
peplib/stats.py       # composition, histogramme de longueurs, KL et divergence de Jensen-Shannon, entropie
peplib/pwm.py         # PWM alignée depuis les extrémités N et C, entropie et information par position
peplib/descriptors.py # charge nette, GRAVY, moment hydrophobe, indice de Boman, point isoélectrique, masse
peplib/dedup.py       # doublons exacts, index inversé de k-mers, distance de Levenshtein, identité
tests/                # un fichier de tests par module
notebooks/NB01_peptides.ipynb   # appelle peplib ; aucune logique métier copiée dans le notebook
models/, results/
```

---

### Workflow Git et collaboration

**1. Chaque personne crée sa propre branche**  
Partir de `main`. Nommer la branche par rôle ou nom : `feat/yassine-donnees`, `feat/sofiene-modeles`, `feat/youssef-evaluation`, `feat/moi-git-arbo`.

**2. Chaque personne ouvre une Pull Request (PR) vers `main`**  
Quand une tâche est finie, on ouvre une PR. La PR doit être relue par **un autre étudiant que l'auteur** (règle du spec : chaque test et chaque fonction est relu par un autre étudiant).

**3. Toi seul merges vers `main`**  
Tu es le seul avec le droit de merger sur `main`. Tu vérifies que les tests passent (`pytest`) avant de merger. Tu gardes le contrôle sur ce qui entre dans la branche principale.

**Pourquoi ?** Cela empêche qu'un étudiant casse `main` pour les 3 autres, et ça force la relecture croisée exigée par le spec.

---

### Roadmap séquentielle des tâches T1–T29

1. **Moi — T1** : création du dépôt Git et de l'arborescence  
   → puis **Youssef — T2** : installation de l'environnement et `requirements.txt`  
   → puis **Youssef — T3** : écriture de `set_seed(42)`  
   → puis **Sofiene — T4** : saisie de `constants.py`  
   → puis **Youssef — T5** : premier test pytest (constantes)  
   → puis **Tous — T6** : lecture de fond (~2 h)  
   → puis **Yassine — T7** : requête REST UniProt (peptides 10–50 résidus)  
   → puis **Yassine — T8** : base secondaire (APD3 / DBAASP / DRAMP)  
   → puis **Yassine — T9** : source non-pharmaceutique (BIOPEP-UWM / BACTIBASE)  
   → puis **Yassine — T10** : rédaction de `PROVENANCE.md`  
   → puis **Yassine — T11** : production d'un CSV par source (`sequence`, `source`, `id`)  
   → puis **Yassine — T12** : chargement pandas, `dropna()` avant `astype(str)`  
   → puis **Yassine — T13** : filtre regex et fenêtre de longueur 10–50  
   → puis **Yassine — T14** : table des rejets par cause  
   → puis **Youssef — T15** : test du filtrage (NaN, vide, minuscules, X, U, longueurs 9 et 51)  
   → puis **Yassine — T16** : histogrammes et composition `f(a)`  
   → puis **Yassine — T17** : divergence de Jensen-Shannon entre sources  
   → puis **Sofiene — T18** : descripteurs écrits à la main (charge, pI, GRAVY, moment hydrophobe, Boman, masse)  
   → puis **Youssef — T19** : test des descripteurs sur peptides calculables à la main  
   → puis **Sofiene — T20** : PWM alignée depuis les extrémités N et C  
   → puis **Sofiene — T21** : figures (distribution, heatmap PWM, courbe d'information)  
   → puis **Yassine — T22** : doublons exacts et taux de redondance  
   → puis **Sofiene — T23** : quasi-doublons (index 3-mers + Levenshtein)  
   → puis **Youssef — T24** : test de la redondance + entropie de composition  
   → puis **Yassine — T25** : taux de quasi-doublons à 90 % d'identité  
   → puis **Tous — T26** : assemblage du notebook `NB01` version peptides  
   → puis **Tous — T27** : rédaction de la demi-page d'analyse  
   → puis **Tous — T28** : séminaire de fin de semaine (30 min)  
   → puis **Tous — T29** : tag `v0.1` et relance complète depuis `data/raw/` ; **Moi** merge vers `main` uniquement après que tous les tests passent.

---

## 4. Liste des tâches (T1–T29) avec détails techniques

### A. Environnement et constantes

**T1 — Créer le dépôt Git et l'arborescence**  
Créer `git init`, branches `main` et branches de feature par rôle. Arborescence : `data/`, `notebooks/`, `models/`, `results/`, `peplib/`, `tests/`. Fichier `.gitignore` excluant `__pycache__/`, `*.pyc`, `.venv/`, `data/raw/*.zip` non versionnés (ou versionnés avec LFS si > 50 Mo). `README.md` avec ordre d'exécution, versions, seeds, commande `pytest`.

**T2 — Installer l'environnement et figer `requirements.txt`**  
Python 3.10 (venv obligatoire). Paquets : `numpy`, `scipy`, `pandas`, `matplotlib`, `torch>=2.0`, `pytest`, `requests`. Gel exact avec `pip freeze > requirements.txt`. Détecter GPU : `torch.cuda.is_available()` ; noter dans le README (CPU uniquement attendu ce mois-ci, pas d'entraînement).

**T3 — Écrire `set_seed(42)` et documenter les seeds**  
Fonction dans `peplib/__init__.py` ou module dédié :
```python
import random, numpy as np, torch
def set_seed(seed: int = 42) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
```
**À quoi ça sert en mots simples :** `set_seed(42)` **verrouille le hasard**. Dans ce pipeline, l'ordinateur doit faire des choix aléatoires (mélanger les données, séparer en train/test/validation, etc.). Sans cette fonction, si vous relancez le code demain, vous obtiendrez des résultats différents parce que le hasard a choisi autre chose. Avec `set_seed(42)` : le "dé" est toujours lancé de la même manière ; les 4 étudiants obtiennent exactement les mêmes comptages ; le pipeline est reproductible : on peut relancer depuis `data/raw/` et retrouver les mêmes fichiers propres.

**Validation :** deux appels successifs redonnent les mêmes tirages aléatoires (test unitaire).

**T4 — Saisir `constants.py` depuis la littérature**  
Tables avec référence citée en commentaire pour chaque valeur :

- **20 acides aminés** : `ACDEFGHIKLMNPQRSTVWY`
- **Masse résiduelle monoisotopique** (Da) : poids de chaque acide aminé sans molécule d'eau.  
  A 71.03711, R 156.10111, N 114.04293, D 115.02694, C 103.00919, E 129.04259, Q 128.05858, G 57.02146, H 137.05891, I 113.08406, L 113.08406, K 128.09496, M 131.04049, F 147.06841, P 97.05276, S 87.03203, T 101.04768, W 186.07931, Y 163.06333, V 99.06841
- **pKa chaîne latérale** : acidité/basicité de la chaîne latérale.  
  D 3.9, E 4.07, C 8.18, Y 10.46, H 6.04, K 10.54, R 12.48
- **pKa extrémités** : N-terminale 9.69, C-terminale 2.34 (justifier le choix, ex. EMBOSS/ExPASy)
- **Échelle Kyte-Doolittle** : mesure de l'hydrophobie (affinité pour l'eau). Valeur positive = hydrophobe, négative = hydrophile.  
  A 1.8, R -4.5, N -3.5, D -3.5, C 2.5, Q -3.5, E -3.5, G -0.4, H -3.2, I 4.5, L 3.8, K -3.9, M 1.9, F 2.8, P -1.6, S -0.8, T -0.7, W -0.9, Y -1.3, V 4.2
- **Échelle Eisenberg (moment hydrophobe)** : valeurs normalisées par résidu pour calculer l'amphipathicité d'une hélice α.

**À quoi ça sert dans le pipeline :** `constants.py` est la table de référence du projet. Tous les autres modules ont besoin de ces nombres pour calculer quoi que ce soit :
- `descriptors.py` calcule la GRAVY, la charge, le pI, la masse → il puise dans `constants.py`
- `stats.py` calcule la composition → il utilise la liste des 20 acides aminés
- `pwm.py` utilise la composition de fond → elle vient de `constants.py`

Validation : test T5 vérifie 20 acides, masses positives, 20 valeurs par échelle.

**T5 — Premier test pytest : constantes**  
- `len(AA_LIST) == 20`
- tous les masses > 0
- `len(KYTE_DOOLITTLE) == 20`, `len(EISENBERG) == 20`
- `MASS["G"] == 57.02146`, `MASS["A"] == 71.03711`

**T6 — Lecture de fond (tous)**  
~2 h : acides aminés, liaison peptidique, peptide contre protéine, peptides antimicrobiens (cationiques, amphipathiques). Chacun explique en 2 minutes pourquoi les AMP sont cationiques et amphipathiques.  
Visionnage de fond optionnel : playlists sur la chimie des peptides, les acides aminés et la structure secondaire (matplotlib uniquement, pas de visualisation 3D).

**Ressources web pour T6 :**

* **Acides aminés et liaison peptidique :**
  - CultureSciences-Chimie (ENS) : *Les acides aminés et la synthèse peptidique* — https://culturesciences.chimie.ens.fr/thematiques/chimie-organique/methodes-et-outils/les-acides-amines-et-la-synthese-peptidique
  - OpenStax (Organic Chemistry) : *26.4 Peptides and Proteins* — https://openstax.org/books/organic-chemistry/pages/26-4-peptides-and-proteins
  - Prepas-SVT : *Peptide… ou protéine* (différence peptide / protéine) — https://www.prepas-svt.fr/wp-content/uploads/2023/01/2023-QR61-peptide%E2%80%A6-ou-proteine.pdf

* **Structure secondaire (hélice α, feuillet β) :**
  - EMBL-EBI : *Secondary structure: α-helices and β-sheets* — https://www.ebi.ac.uk/training/online/courses/foundations-protein-structure/principles-of-protein-folding-and-architecture/secondary-structure-%ce%b1-helices-and-%ce%b2-sheets/
  - Réussir SVT : *Leçon 3.1 Structures secondaires* — https://www.reussir-svt.com/cours/structure-des-proteines/lecon/structures-secondaires-helices-feuillets/
  - OpenStax : *26.9 Protein Structure* — https://openstax.org/books/organic-chemistry/pages/26-9-protein-structure

* **Peptides antimicrobiens (cationiques et amphipathiques) :**
  - ScienceDirect : *Peptides antimicrobiens cationiques : de l’étude de l’immunité innée à la production de médicaments* — https://www.sciencedirect.com/science/article/pii/S0399077X06003453
  - PMC (PubMed Central) : *Latest developments on the mechanism of action of membrane disrupting peptides* — https://pmc.ncbi.nlm.nih.gov/articles/PMC10244799/
  - Cambridge Core : *Understanding membrane-active antimicrobial peptides* — https://www.cambridge.org/core/journals/quarterly-reviews-of-biophysics/article/understanding-membraneactive-antimicrobial-peptides/D82F56B43F7165B1DC11B51F120875B7

---

### B. Acquisition des séquences et provenance

**T7 — UniProt : requête REST sur peptides 10–50 résidus**  
Endpoint : `https://rest.uniprot.org/uniprotkb/search?query=...&format=tsv&fields=accession,sequence,length&size=500` avec pagination via en-tête `Link` (cursor). Filtre `reviewed:true` d'abord, puis `reviewed:false` si nécessaire. Sauvegarder le fichier brut intact (`data/raw/uniprot_*.tsv`). Consigner le nombre de lignes.

**À quoi ça sert en mots simples :** C'est aller chercher des peptides sur Internet. On interroge UniProt (une énorme base de données publique de protéines et peptides) via une adresse web (l'API REST). On lui demande : « donne-moi uniquement les séquences de 10 à 50 acides aminés, avec leur identifiant, leur séquence et leur longueur ». Comme il y a trop de résultats pour un seul envoi, on les reçoit par paquets de 500 (pagination) jusqu'à tout récupérer. On privilégie d'abord les entrées `reviewed:true` (vérifiées à la main par des experts), puis `reviewed:false` si besoin (non vérifiées mais utiles). On sauvegarde le fichier tel quel dans `data/raw/uniprot_*.tsv` sans jamais le modifier. On note combien de lignes on a récupérées.

**À quoi ça sert dans le pipeline :** C'est la matière première. Sans T7, il n'y a pas de peptides à nettoyer, ni à analyser. Toutes les étapes suivantes (filtrage, descripteurs, statistiques) partent de ce fichier brut.

**T8 — Une base de données secondaire (DBAASP, APD3 ou DRAMP)**  
Vérifier conditions d'usage et de redistribution avant téléchargement. Sauvegarder brut. Consigner licence dans `PROVENANCE.md`.

**Différence entre les types de bases secondaires :**
- **APD3** (*Antimicrobial Peptide Database*) : peptides antimicrobiens avec séquence et activité biologique (profil pharma).
- **DBAASP** : *Database of Antimicrobial Activity and Structure of Peptides* — activité et structure de peptides antimicrobiens.
- **DRAMP** : *Database of Antimicrobial Peptides and Resistance* — peptides antimicrobiens et résistance.

Recommandation : prendre **APD3** comme source principale cette semaine. Si sa licence bloque le téléchargement ou la redistribution, basculer sur **DBAASP** ou **DRAMP**. Consigner la licence choisie dans `PROVENANCE.md`.

**T9 — Une source non pharmaceutique**  
Peptides alimentaires bioactifs ou bactériocines de biopreservation. Vérifier conditions d'usage. Sauvegarder brut.

**Différence entre les types de sources :**
- **BIOPEP-UWM** : peptides bioactifs d'origine alimentaire (Université de technologie de Varsovie). Source non-pharma alimentaire.
- **BACTIBASE** : bactériocines d'acide lactique (*LAB*), usage en biopreservation. Source non-pharma de conservation alimentaire.

Recommandation : **BIOPEP-UWM** pour les peptides alimentaires bioactifs ; **BACTIBASE** si le groupe se concentre sur les bactériocines de biopreservation. Vérifier les conditions d’usage et de redistribution avant téléchargement, et consigner la licence dans `PROVENANCE.md`.

**T10 — Rédiger `PROVENANCE.md`**  
Une entrée par source : nom, date, requête ou URL, version, licence, nombre de lignes, empreinte SHA-256 (`hashlib.sha256`). Téléchargement manuel noté s'il y en a un.

**T11 — Produire un CSV par source**  
Colonnes : `sequence`, `source`, `id`. Lignes comptées, identiques au brut (avant filtrage). `pandas.read_csv` / `to_csv`.

---

### C. Filtrage exact

**T12 — Charger avec pandas, `dropna()` avant tout `astype(str)`**  
`astype(str)` avant `dropna()` transforme `NaN` en `"NAN"` qui passe la regex. Ordre : `dropna()` → `strip()` → `upper()` → filtres. Aucune valeur manquante n'atteint la regex.

**T13 — Filtrer par regex et fenêtre de longueur**  
Regex : `^[ACDEFGHIKLMNPQRSTVWY]+$`. Longueur : `10 <= len(seq) <= 50`. 100 % des séquences conservées passent la regex.

**T14 — Table des rejets par cause**  
Ordre d'application documenté (chaque séquence comptée dans sa **première** cause) :
1. manquante (`NaN` après `dropna`)
2. vide (`""`)
3. lettres ambiguës (`X`, `B`, `Z`, `J`)
4. lettres rares (`U`, `O`)
5. autres caractères (chiffres, espaces, ponctuation)
6. trop courte (< 10)
7. trop longue (> 50)

Exporter les CSV traités dans `data/processed/`. Vérifier : `conservés + rejetés = total brut`.

**T15 — Tester le filtrage**  
Cas : `NaN`, chaîne vide, minuscules (`acdef...`), `X`, `U`, longueur 9, longueur 51. Rejets attendus par cause.

---

### D. Statistiques, descripteurs, PWM et information

**T16 — Histogrammes et composition**  
Histogrammes des longueurs par source (`matplotlib`). Composition `f(a) = compte(a) / N_total`. Test : `sum(f(a)) = 1` (tolérance `1e-12`) sur chaque source.

**T17 — Divergence de Jensen-Shannon entre sources**  
Matrice 3 × 3 (ou N × N selon nombre de sources).  
`m = (p + q) / 2` ; `JS(p, q) = 0.5 * KL(p || m) + 0.5 * KL(q || m)`.  
Propriétés : `JS(p, p) = 0` ; bornée par `ln 2` (nats).

**À quoi ça sert en mots simples :** T17 compare la **"recette"** (proportion de chaque acide aminé) de chaque source de peptides. Résultat 0 : les deux sources sont identiques. Résultat élevé (max ≈ 0,69) : elles sont très différentes. Sert à répondre à la question d'analyse — *"Les peptides pharma (AMP) et non-pharma sont-ils composés différemment ?"* Si JS est élevé entre UniProt et BIOPEP, les deux familles ne se ressemblent pas.

**T18 — Descripteurs écrits à la main** (dépend de T4)
- **Charge nette à pH 7** : résidus positifs `1 / (1 + 10^(pH - pKa))` + résidus négatifs `1 / (1 + 10^(pKa - pH))`, extrémités N et C comprises.
- **Point isoélectrique (pI)** : pH tel que `charge(pH) = 0`, trouvé par bissection (charge strictement décroissante en pH).
- **GRAVY** : moyenne des valeurs Kyte-Doolittle sur les résidus.
- **Moment hydrophobe** : `μ_H = sqrt( [Σ H_i cos(k·i)]^2 + [Σ H_i sin(k·i)]^2 )`, `k = 100°` pour une hélice α.
- **Indice de Boman** : somme des énergies de liaison aux protéines des résidus / nombre de résidus (table Boman 2003).
- **Masse monoisotopique** : Σ masses résiduelles + masse H₂O (18.01056 Da).

**À quoi ça sert en mots simples :** Pour chaque peptide, on calcule des nombres simples qui le décrivent : son poids, sa charge électrique, s'il aime ou déteste l'eau (GRAVY), sa forme en hélice (moment hydrophobe), sa capacité à se lier aux protéines (Boman), et son point de neutralité (pI). Rôle dans le pipeline : donner une carte d'identité physico-chimique à chaque peptide.

**T19 — Tester les descripteurs sur peptides calculables à la main**  
Valeurs attendues (section Tests obligatoires ci-dessous).

**À quoi ça sert en mots simples :** On vérifie les calculs sur des petits peptides dont on connaît le résultat à la main (ex: GRAVY de "A" = 1,8). Rôle : s'assurer que le code ne dit pas n'importe quoi avant de l'utiliser sur tout le corpus.

**T20 — PWM alignée depuis les extrémités N et C**  
Pas de troncature à la longueur minimale. Aligner depuis N-terminale (position 0 = premier résidu) et depuis C-terminale (position -1 = dernier résidu).  
`PWM(i, a) = log[ p̂ᵢ(a) / p_bg(a) ]`, `p_bg` = composition globale. Entrées nulles masquées ou pseudo-count justifié par écrit.  
Information par position : `I(i) = log 20 - [H(i) + e(n)]`, correction de petit échantillon `e(n) = (s - 1) / (2 * ln 2 * n)`, `s = 20`, `n` = nombre de séquences à la position `i` (Schneider & Stephens, 1990). Masquer `log(0)` exactement (`0 log 0 = 0`).

**À quoi ça sert en mots simples :** On regarde quelles lettres reviennent souvent au début et à la fin des peptides, position par position, sans tronquer les séquences. Rôle : repérer des motifs conservés aux bords des peptides.

**T21 — Figures**  
Distribution des descripteurs par source, heatmap de la PWM, courbe d'information par position. Axes légendés, lisibles.

**À quoi ça sert en mots simples :** On dessine des graphiques : histogrammes des descripteurs, heatmap de la PWM, courbe d'information par position. Rôle : voir visuellement si les sources de peptides se ressemblent ou diffèrent.

---

### E. Redondance, analyse et clôture

**T22 — Doublons exacts et taux de redondance**  
Ensemble Python (`set`). Taux par source et entre sources. `conservés uniques / conservés totaux`.

**À quoi ça sert en mots simples :** On enlève les copies identiques (même séquence lettre pour lettre) pour ne pas compter plusieurs fois le même peptide. On calcule le taux de redondance par source et entre sources. Sert à : savoir combien de peptides uniques on a vraiment.

**T23 — Quasi-doublons**  
Index inversé des 3-mers : dictionnaire `kmer -> ensemble d'IDs`. Deux séquences sont candidates si elles partagent au moins un 3-mer. Deux séquences sans 3-mer commun ne sont pas candidates. Distance de Levenshtein exacte sur les candidats. Identité : `1 - d_Lev(a, b) / max(|a|, |b|)`.  
Test : `Levenshtein("kitten", "sitting") == 3`.

**À quoi ça sert en mots simples :** On cherche les séquences très proches (pas identiques, mais très similaires) : on repère rapidement les candidats qui partagent un petit bout de 3 lettres (3-mers), puis on mesure leur distance exacte (Levenshtein). Sert à : détecter les peptides qui se ressemblent à 90 % ou plus.

**T24 — Tester la redondance et calculer l'entropie de composition**  
`H = - Σ f(a) * log f(a)`. Séquences identiques candidates ; séquences sans 3-mer commun non candidates.

**À quoi ça sert en mots simples :** On vérifie que la détection des doublons marche bien, et on calcule la diversité des acides aminés dans le corpus (entropie). Sert à : valider le code de détection et mesurer si le corpus est varié ou répétitif.

**T25 — Taux de quasi-doublons à 90 % d'identité**  
Rapporté + justification écrite pourquoi le clustering complet (union-find) est écrit au mois 2 (dépend de T23).

**À quoi ça sert en mots simples :** On rapporte combien de peptides sont quasi-identiques à 90 %, et on explique par écrit pourquoi le regroupement complet (clustering union-find) sera fait le mois 2, pas maintenant. Sert à : quantifier la similarité du corpus et justifier le report du clustering.

**T26 — Assembler `NB01` version peptides**  
Notebook Jupyter appelant `peplib`. Aucune logique métier copiée dans les cellules (imports + appels uniquement).

**À quoi ça sert en mots simples :** On crée un notebook Jupyter qui appelle simplement `peplib` (imports et appels uniquement). Aucun calcul n'est recopié dans les cellules du notebook. Sert à : montrer le pipeline en action, sans mélanger la logique métier et la présentation.

**T27 — Rédiger la demi-page d'analyse**  
Un paragraphe par étudiant. Questions traitées :
- Les distributions de longueurs sont-elles unimodales ? Les sources diffèrent-elles ?
- Quelles différences de composition entre peptides pharma (AMP) et non-pharma ?
- Que dit la divergence de Jensen-Shannon ?
- Quelles positions des extrémités N et C sont conservées ?
- Quelle fraction du corpus est redondante, exactement et à 90 % d'identité ?
- Quelles conséquences pour le split du mois 2 ?

**T28 — Séminaire de fin de semaine (30 min)**  
0–10 min : corpus (Yassine Achour)  
10–20 min : composition, descripteurs, PWM, information (Sofiene Ouanes)  
20–27 min : tests, redondance, taux de quasi-doublons (Youssef Jarraya*)  
27–30 min : discussion, ce qui change pour le split mois 2

**T29 — Poser le tag `v0.1` ; relance complète**  
Chacun des quatre relance tout depuis `data/raw/` ; mêmes comptages ; tous les tests passent.

---

## 5. Formules à implémenter

| Quantité | Définition exacte |
|----------|-------------------|
| Composition | `f(a) = compte(a) / N_total` ; `Σ f(a) = 1` |
| Charge nette à pH | résidus positifs `1/(1+10^(pH-pKa))` + résidus négatifs `1/(1+10^(pKa-pH))`, extrémités comprises |
| pI | pH tel que `charge(pH) = 0` par bissection |
| GRAVY | moyenne Kyte-Doolittle sur les résidus |
| Moment hydrophobe | `μ_H = sqrt( (Σ H_i cos(k·i))² + (Σ H_i sin(k·i))² )`, `k = 100°` |
| Indice de Boman | Σ énergies de liaison résidu / longueur |
| Masse monoisotopique | Σ masses résiduelles + 18.01056 |
| Entropie position | `H(i) = -Σ p̂ᵢ(a) log p̂ᵢ(a)`, `0 log 0 = 0` masqué exactement |
| Information position | `I(i) = log 20 - [H(i) + e(n)]`, `e(n) = (20-1)/(2 ln 2 * n)` |
| PWM | `log[p̂ᵢ(a) / p_bg(a)]`, entrées nulles masquées ou pseudo-count justifié |
| JS divergence | `0.5*KL(p||m) + 0.5*KL(q||m)`, `m = (p+q)/2`, bornée par `ln 2` |
| Identité | `1 - d_Lev(a,b) / max(|a|,|b|)` |
| Entropie de composition | `H = -Σ f(a) log f(a)` |

---

## 6. Tests obligatoires

| Test | Résultat attendu |
|------|------------------|
| `GRAVY("A")` | `1.8` |
| `GRAVY("GG")` | `-0.4` |
| `GRAVY("AG")` | `(1.8 + (-0.4)) / 2 = 0.7` |
| Masse `"G"` | `~75.0320` (57.02146 + 18.01056) |
| Masse `"AA"` | `~160.0848` (2 × 71.03711 + 18.01056) |
| Charge `"DE"` et `"KK"` à pH 7 | `DE` négative, `KK` positive |
| Charge `"G"` à pH 1 | valeur à la main, égale au code à 10⁻³ près |
| pI | `charge(pI) = 0` à 10⁻³ près ; charge strictement décroissante en pH |
| Moment hydrophobe, k = 0 | égal à `|H|` ; un résidu : `|H|` |
| Information, séquences identiques | `H = 0`, donc `I = log 20 ≈ 4.3219` avant correction ; `I = log 20 - e(n)` après |
| Information, 20 séquences d'un AA chacune | `H = log 20`, donc `I = 0` avant correction |
| JS | `JS(p,p) = 0` ; `JS(δ_A, δ_C) = ln 2` |
| Composition | `Σ f = 1` à 10⁻¹² près, sur chaque source |
| Filtrage | `NaN`, vide, minuscules, `X`, `U`, longueur 9 et 51 : rejets attendus ; conservés + rejetés = total brut |
| Levenshtein | `("kitten", "sitting") == 3` |
| Index k-mers | identiques = candidates ; sans 3-mer commun = non candidates |

---

## 7. Pièges à éviter

- `astype(str)` avant `dropna()` fait passer `"NAN"` pour une séquence valide → `dropna()` d'abord.
- PWM : ne **pas** tronquer à la longueur minimale ; aligner depuis les extrémités N et C.
- `log(0)` : masquer les entrées nulles ou justifier un pseudo-count ; ne pas laisser `-inf` dans une figure.
- Information sur petits échantillons : sans correction `e(n)`, `I(i)` est biaisée vers le haut quand `n` est petit.
- Constantes : une erreur de saisie dans une échelle fausse tous les résultats → relecture croisée et test sur valeurs connues.
- Notebook contre module : la logique vit dans `peplib` et ses tests, pas dans des cellules copiées.
- Fichiers bruts : ne jamais les modifier ; chaque fichier traité se reproduit par code.
- **Visualisation 3D interdite** : PyMOL, NGL Viewer, AlphaFold et autres prédicteurs de structure 3D sont interdits. Seuls matplotlib (2D) et le modèle HP sur réseau 2D sont autorisés.

---

## 8. Biosécurité, provenance et licences

UniProt contient des toxines (conotoxines, par exemple). Signaler ces séquences ou les exclure du corpus génératif, et ne pas diffuser de séquences générées de type toxine. Vérifier pour chaque base ses conditions d'usage et de redistribution avant tout partage du corpus ; consigner le résultat dans `PROVENANCE.md`.

---

## 9. Critères de réussite de la semaine

- Les 29 tâches sont cochées, ou les restantes reportées avec une raison écrite.
- 100 % des séquences conservées passent la regex ; conservés + rejetés = total brut.
- Tous les tests de la section 6 passent.
- Fichiers bruts avec empreintes et `PROVENANCE.md` complet, licences vérifiées ; toxines signalées ou exclues.
- Les quatre étudiants obtiennent les mêmes comptages en relançant depuis `data/raw/`.

---

## 10. Ce que le superchef (moi) doit présenter à M. Abdallah

À la fin de la semaine, M. Abdallah doit voir des **visuels et des documents**, pas du code. Le superchef doit être capable de les expliquer sans montrer le code source.

### Visuels à présenter
- **Histogrammes des longueurs** de peptides par source (UniProt, APD3, BIOPEP-UWM…)
- **Heatmap de la PWM** (quelles lettres reviennent souvent au début et à la fin des peptides, position par position)
- **Courbe d'information par position** (quelles positions des extrémités N et C sont conservées)
- **Graphiques de distribution des descripteurs** par source (GRAVY, charge nette, masse, pI, moment hydrophobe)

### Tableaux et documents à présenter
- **Table des rejets par cause** : combien de séquences rejetées et pourquoi (trop courte, lettre X/U, vide, NaN, etc.)
- **Matrice de divergence de Jensen-Shannon (JS)** : dans quelle mesure les sources de peptides se ressemblent ou diffèrent
- **Taux de redondance** : pourcentage de doublons exacts et de quasi-doublons à 90 % d'identité
- **Fichier `PROVENANCE.md`** : provenance de chaque source, date, licence, empreinte SHA-256
- **Demi-page d'analyse** : un paragraphe par étudiant répondant aux 6 questions d'analyse

### Preuve de reproductibilité à montrer
- **Tag Git `v0.1`** posé sur le dépôt
- Les 4 étudiants ont relancé tout depuis `data/raw/` et ont obtenu **les mêmes comptages** (même nombre de conservés, même nombre de rejetés par cause)

**On doit savoir expliquer au superchef chacun de ces éléments à M. Abdallah sans ouvrir le code.**

---

## 11. Commandes de démarrage rapide (à exécuter dans l'ordre)

Ouvrir un terminal à la racine du projet :

```bash
# 1. Créer l'environnement virtuel
python -m venv .venv

# 2. L'activer (Windows PowerShell)
.venv\Scripts\Activate.ps1

# 3. Installer les dépendances
pip install -r requirements.txt

# 4. Initialiser le dépôt Git (fait par le chef, T1)
git init
git branch -m main

# 5. Fixer la seed et lancer les tests
pytest
```

Résultat attendu : tous les tests passent, aucun échec.

---

## 12. Template `PROVENANCE.md` (à remplir pour chaque source)

Copier-coller ce bloc pour chaque base de données utilisée :

```markdown
# Source : [NOM DE LA BASE]

- **Nom de la source** : [ex. UniProtKB / APD3 / BIOPEP-UWM / BACTIBASE]
- **URL de requête ou de téléchargement** : [ex. https://rest.uniprot.org/uniprotkb/search?query=...]
- **Date de téléchargement** : [JJ/MM/AAAA]
- **Version de la base** : [si indiqué, sinon "non spécifiée"]
- **Conditions d'usage / licence** : [ex. CC-BY 4.0, licence propriétaire, interdit de redistribution...]
- **Nombre de lignes brutes** : [nombre]
- **Empreinte SHA-256 du fichier brut** : [ex. a3f5c8...]
- **Fichier brut sauvegardé sous** : `data/raw/[nom_fichier].tsv`
- **Téléchargement manuel ?** : Oui / Non
- **Notes** : [toute remarque sur la qualité, les filtres appliqués à la source, les toxines détectées, etc.]
```

Exemple rempli pour UniProt :

```markdown
# Source : UniProtKB

- **Nom de la source** : UniProtKB (reviewed:true)
- **URL de requête** : https://rest.uniprot.org/uniprotkb/search?query=length:[10 TO 50] AND reviewed:true&format=tsv&fields=accession,sequence,length&size=500
- **Date de téléchargement** : 15/01/2026
- **Version de la base** : non spécifiée
- **Conditions d'usage / licence** : UniProt Terms of Use (consulté le 15/01/2026)
- **Nombre de lignes brutes** : 12 450
- **Empreinte SHA-256 du fichier brut** : a1b2c3d4e5f6...
- **Fichier brut sauvegardé sous** : `data/raw/uniprot_reviewed_20260115.tsv`
- **Téléchargement manuel ?** : Non
- **Notes** : Pagination gérée via en-tête Link (cursor). Toxines (conotoxines) signalées et exclues du corpus génératif.
```