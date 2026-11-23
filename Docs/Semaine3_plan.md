# Semaine 3 — Modèles de base from scratch : LSTM et Transformer décodeur, génération et évaluation (spécification adaptée 4 étudiants)

**Groupe : 4 étudiants (3 développeurs + 1 chef/coordinateur)**
**Encadrement : A. Ben Abdallah, ISIMS**
**Phase C du plan annuel (mois 2, édition semaine 3 ; tâches T48–T68), version peptides**

> **Contrainte de calendrier.** Seulement **4 jours de travail (vendredi 30/10 → lundi 02/11)**, réunion avec M. Abdallah **mardi 03/11**. Yassine Achour est indisponible **vendredi dès 16:00** et **samedi dès 12:00**. Le plan réel est la **section 3 bis** ; le périmètre y est réduit et classé par priorité.

> **Convention de lecture.** Les valeurs entre crochets `[à vérifier]` sont reprises de documents de travail précédents et doivent être relues dans `corpus_manifest.json` et `results/` avant toute présentation. Les hyperparamètres proposés sont des **points de départ à justifier**, pas des vérités.

---

## 1. Objectif et livrables

Les semaines 1 et 2 ont produit un corpus v0.2 (UniProt + APD3 + AMPChallenge, fusionné sans doublon), des **familles** (clusters à 90 %), un **split étanche** train/val/test, une tokenisation V = 23 et trois **baselines** (aléatoire, unigramme, bigramme). La semaine 3 est la première où un **réseau de neurones apprend**. On entraîne deux modèles de base **sans aucune condition**, entièrement codés par nous, on génère des peptides, et on les évalue **honnêtement** (sans fuite, sans tricher sur le test).

### Contrainte fondatrice du projet (AGENTS.md)
**Aucun poids pré-entraîné, jamais.** Tout est codé et entraîné par le groupe, avec `nn.Module` et l'autograd de PyTorch uniquement. Cela exclut les options « classiques » d'un projet LLM :

| Option | Verdict | Raison à savoir défendre |
|--------|---------|--------------------------|
| Prompt engineering sur LLM / API | ❌ | poids pré-entraînés, aucune maîtrise de l'algorithme, non reproductible (seed 42 impossible) |
| RAG + base vectorielle | ❌ | hors contrainte **et** conceptuellement un anti-modèle : il recombine des voisins du train, c'est-à-dire exactement la fuite que le split par familles a supprimée |
| Fine-tuning LoRA/QLoRA d'un LLM ouvert | ❌ | poids pré-entraînés ; inadapté à des séquences de 10–50 acides aminés avec un vocabulaire de 23 tokens |
| **LSTM char-level from scratch** | ✅ | baseline neuronale robuste sur petit corpus |
| **Transformer décodeur from scratch** | ✅ | comparaison à budget égal avec le LSTM |
| Diffusion discrète, VAE | ⏳ | prévus plus tard (ordre du planning annuel figé) : ne pas ouvrir une nouvelle famille avant le jalon S3 |

Cette exclusion est une **limite déclarée** dans le rapport, pas un oubli.

### Qu'est-ce qu'un modèle autorégressif ?
Un modèle qui prédit **le prochain token** connaissant les précédents : `p(x₁…x_L) = Π p(x_t | x_<t)`. Pour générer : on part de `BOS`, on tire un token, on le réinjecte, et on s'arrête quand le modèle tire `EOS` (ou à 50 résidus). Entraîner = minimiser l'**entropie croisée** (la « surprise » moyenne du modèle devant les vrais peptides du train).

### Qu'est-ce que la perplexité (PP) ?
`PP = exp(entropie croisée moyenne en nats/token)`. Lecture intuitive : « le modèle hésite en moyenne entre PP tokens également probables ». Plus bas = mieux. Un tirage uniforme sur 21 tokens (20 acides aminés + EOS) donne PP = 21 (ln 21 = 3,045). La baseline unigramme de la semaine 2 donne **PP = 19,49 `[à vérifier]`** : c'est la **barrière à franchir**.

### Qu'est-ce qu'un test verrouillé ?
Le test set sert **une seule fois**, à la fin, pour rapporter un chiffre final. Tout réglage (architecture, taux d'apprentissage, température d'échantillonnage, early stopping) se fait sur **val**. Si on regarde le test pour choisir, il n'est plus un test : c'est une fuite, de la même famille que celle que les clusters ont supprimée.

**Résumé en une phrase :** On gèle le test → on charge les données → on code et entraîne un LSTM puis un Transformer (3 graines chacun) → on génère 1 000 peptides par modèle → on mesure perplexité, validité, unicité, distance au train et proxys physico-chimiques contre des lignes de référence → on ouvre le test **une seule fois** → on rédige.

**Livrables :**
- `peplib/config.py` + `configs/*.yaml` (configurations versionnées, une par modèle)
- `peplib/dataset.py` (`make_loader` existant, PR #20, étendu : padding, masques) `[chemin à confirmer]`
- `peplib/models/lstm.py` et `peplib/models/transformer.py` (`nn.Module`, from scratch)
- `peplib/train.py` (boucle d'entraînement, early stopping, checkpoints, logs CSV, reprise)
- `peplib/metrics_lm.py` (entropie croisée, perplexité, partagée avec les baselines)
- `peplib/sampling.py` (température, top-k, top-p, masques PAD/BOS)
- `peplib/eval_gen.py` (validité, unicité, longueurs, composition, descripteurs, distance au train)
- `data/generated/` : CSV des séquences générées (non distribuées), `results/week3/` : tableaux et figures
- `corpus_manifest.json` étendu (hashes des checkpoints, configs, seeds) + `PROVENANCE.md` complété
- `Docs/Semaine3_rapport.md` (rapport) + tag Git `v0.3`

**Critère de réussite :** les quatre étudiants relancent l'entraînement **CPU** d'un modèle (config + seed) et obtiennent la **même perte, bit à bit**, et les mêmes séquences générées pour un même seed d'échantillonnage.

---

## 2. Répartition des rôles (4 étudiants)

Chaque fonction et chaque test est relu par un autre étudiant que son auteur. Les rôles servent à répartir l'écriture, pas à isoler les connaissances : chacun doit pouvoir expliquer tout le pipeline.

| Rôle | Responsabilités | Tâches principales |
|------|----------------|--------------------|
| **Yassine Achour — Données** | Configurations, chargeur de données, script de lancement des graines (lancées par un autre), stockage et traçabilité des générations, drapeau toxines, manifeste. **Contraint : vendredi < 16:00, samedi < 12:00** | T49, T50, T58 |
| **Youssef Jarraya — Modèles** | LSTM, Transformer, boucle d'entraînement, early stopping, checkpoints, ablations | T52, T53, T54, T55, T56, T62 |
| **Sofiene Ouanes — Évaluation** | Perplexité partagée, échantillonnage, métriques de génération, distance au train, lignes de référence | T51, T57, T59, T60, T61, T63 |
| **Moi — Chef / Coordinateur (L)** | Gel du test, protocole, tableau comparatif, **seule ouverture du test**, rapport, relectures croisées, merge vers `main` | T48, T64, T65, T66, T67, T68 + revue de tous les PR |

---

## 3. Structure du dépôt et des modules

```
data/raw/                 # immuable (inchangé)
data/processed/           # corpus_v0.2.csv, clusters_90.csv, split_train/val/test.csv (inchangés)
data/generated/           # séquences générées (CSV) — NON distribuées, hash dans le manifeste
configs/                  # lstm_base.yaml, transformer_base.yaml (versionnés)
checkpoints/              # poids entraînés — .gitignore, hashes dans le manifeste
peplib/config.py          # chargement/validation des YAML (dataclasses), seeds
peplib/dataset.py         # make_loader (PR #20) + collate (padding, masque) 
peplib/models/lstm.py     # LSTM char-level from scratch
peplib/models/transformer.py  # Transformer décodeur causal from scratch
peplib/train.py           # boucle d'entraînement, early stopping, reprise
peplib/metrics_lm.py      # entropie croisée / perplexité (nats/token), partagée avec baselines.py
peplib/sampling.py        # génération BOS→EOS, température, top-k, top-p
peplib/eval_gen.py        # validité, unicité, longueurs, JS, KS, distance au train
tests/                    # un fichier de tests par module
results/week3/            # tableaux, figures, logs CSV d'entraînement, LOCK du test
corpus_manifest.json      # + configs, seeds, hashes checkpoints et générations
```

**Réutilisation explicite des semaines 1 et 2 :**

| Livrable S1/S2 | Usage en S3 |
|----------------|-------------|
| `corpus_v0.2.csv` filtré (alphabet canonique, 10–50) | corpus d'entraînement ; **re-filtrage des sorties** avec le même `filters.py` (validité) |
| Splits étanches par familles (`split_*.csv`) : 38 336 / 2 095 / 2 144 `[à vérifier]` | le modèle n'apprend que sur train ; val = early stopping et réglages ; test verrouillé |
| `tokens.py` (V = 23 : PAD=0, BOS=1, EOS=2, AA 3..22) | embeddings, décodage des générations (**V ne change pas en S3**) |
| `baselines.py` (aléatoire, unigramme, bigramme) | barrières à battre ; **générateurs de référence** pour valider la chaîne d'évaluation |
| `align.py`, `cluster.py` (Levenshtein 3-mers, identité) et `fuite.py`, `dedup.py` | distance des générations au train, diversité |
| `descriptors.py` (charge, pI, GRAVY, moment hydrophobe) | proxys physico-chimiques des générations |
| Témoins négatifs (séquences mélangées, T45) | ligne de référence « mélangé » |
| `toxin_flags.json` + règles de biosécurité | drapeau toxines sur les générations |
| `make_loader` seedé + test bit-à-bit (PR #20) | base du chargeur |
| `PROVENANCE.md`, `corpus_manifest.json`, Git/PR | traçabilité des résultats |

---

### Workflow Git et collaboration

**1. Chaque personne crée sa propre branche** à partir de `main` : `feat/yassine-config-loader`, `feat/sofiene-lstm`, `feat/sofiene-transformer`, `feat/youssef-sampling-eval`, etc.

**2. Chaque personne ouvre une PR** vers `main`, relue par un autre étudiant que l'auteur.

**3. Seul le coordinateur merge vers `main`**, pytest vert (Python 3.11 **et** 3.14).

**Pourquoi ?** Mêmes raisons qu'en semaines 1–2, plus une règle propre à S3 : **une seule nouveauté à la fois** — chargeur testé avant modèle, modèle testé avant entraînement, entraînement validé avant génération. Un bug dans le masque causal découvert après trois jours d'entraînement coûte trois jours.

---

### Roadmap séquentielle des tâches T48–T68

1. **Moi — T48** : clôture S2 (tag `v0.2` vérifié), gel du test, protocole d'évaluation figé
   → puis **Yassine — T49** : configs YAML + validation, seeds
   → puis **Yassine — T50** : chargeur (`make_loader` étendu : padding, masque) + tests
   → puis **Sofiene — T51** : perplexité partagée (tokens = AA + EOS), reproduction du 19,49
   → puis **Youssef — T52** : LSTM + tests
   → puis **Youssef — T53** : boucle d'entraînement, early stopping, checkpoints, reprise
   → *(en parallèle)* **Sofiene — T57** : échantillonneur, **testé sur les baselines unigramme/bigramme**
   → puis **Youssef (+ Sofiene pour les graines 43/44)  — T54** : entraînement LSTM, 3 graines
   → puis **Youssef — T55** : Transformer décodeur + test du masque causal
   → puis **Youssef (+ Sofiene pour les graines 43/44) — T56** : entraînement Transformer, 3 graines, budget égal
   → puis **Yassine — T58** : stockage des générations, re-filtrage, drapeau toxines
   → puis **Sofiene — T59** : métriques de génération (validité, unicité, longueurs, JS, KS)
   → puis **Sofiene — T60** : distance au train
   → puis **Sofiene — T61** : lignes de référence (baselines générées, val réelle, mélangés)
   → puis **Youssef — T62** (si temps) : ablations (profondeur, dimension, dropout)
   → puis **Sofiene — T63** (si temps) : analyse de mémorisation
   → puis **Moi — T64** : tableau comparatif et figures
   → puis **Moi — T65** : **évaluation test, une seule fois**
   → puis **Tous — T66** : inspection humaine + rapport
   → puis **Tous — T67** : relance complète, tag `v0.3`, **Moi** merge vers `main`
   → puis **Tous — T68** : séminaire interne (30 min) et présentation à M. Abdallah.

**Calendrier initial sur 7 jours (remplacé par la section 3 bis, 4 jours) :**

| Jour | Tâches | Risque principal | Repli |
|------|--------|------------------|-------|
| J1 | T48, T49, T50, T51 (+ chronométrer 1 epoch CPU) | padding/masque mal gérés | masques verrouillés par test bit-à-bit ; PP reproduit à 10⁻⁶ près la baseline |
| J2 | T52, T53, T57 (échantillonneur testé sur baselines) | écart entre PP du modèle et PP des baselines (tokens comptés différemment) | un seul calcul de PP partagé (T51) |
| J3 | T54 (3 graines LSTM), début T55 | overfitting (≈ 38 k séquences) | early stopping strict, dropout, weight decay ; **repli** : réduire la taille du modèle |
| J4 | T55, T56 (3 graines Transformer) | session Colab coupée, divergence | checkpoints + reprise ; petite architecture d'abord ; warmup + clipping |
| J5 | T58, T59, T60, T61 (grille d'échantillonnage sur val) | effondrement de modes (« tout Leu/Lys ») | top-k/top-p, la métrique d'unicité le détecte |
| J6 | T62/T63 (si temps), T64, **T65 (test, une fois)** | regarder le test trop tôt | script de test avec fichier LOCK ; choix signés avant |
| J7 | T66, T67, T68 | relance qui ne reproduit pas | CPU bit-à-bit ; sinon documenter la tolérance GPU |

---

## 3 bis. Plan réduit : 4 jours (vendredi → lundi), réunion mardi

**Contraintes de calendrier :**
- Travail effectif : **vendredi 30/10, samedi 31/10, dimanche 01/11, lundi 02/11** ; **réunion avec M. Abdallah mardi 03/11**.
- **Yassine Achour** : vendredi, plus rien à partir de **16:00** ; samedi, plus rien à partir de **12:00**. Dimanche et lundi : disponible `[à confirmer]`.
- Le plan initial de 7 jours est compressé en 4 : le **périmètre est réduit** (voir « Priorités ») et les tâches de Yassine sont **placées avant ses coupures**, sans qu'aucune tâche bloquante ne dépende de lui après ces heures.

### Règles de planification

1. **Tout ce qui bloque les autres est livré par Yassine avant sa coupure** : configs et chargeur (T49, T50) vendredi avant 14:00 ; stockage des générations, script de lancement des graines et extension du manifeste (T58) samedi avant 12:00, avec un mini-mode d'emploi dans la PR.
2. **Yassine ne lance aucune graine après sa coupure** : il prépare `run_all_seeds` (script unique pour 42/43/44) et Youssef ou Sofiene le lancent. Les entraînements longs tournent **sans surveillance** pendant la nuit et l'après-midi des autres.
3. **Les tâches de Sofiene ne dépendent pas des vrais modèles au début** : échantillonneur et métriques sont développés et validés sur les **baselines** et sur **val réelle/mélangée** ; les modèles s'y branchent ensuite.
4. **Les réglages sont figés dimanche soir** (choix signés par commit) ; le test est ouvert **lundi, une seule fois**.
5. **Décision CPU ou Colab T4** vendredi soir, après chronométrage du premier epoch : si le CPU rend un entraînement de 30 epochs infaisable avant dimanche midi, Colab dès samedi (checkpoints sur Drive pour la reprise). Le seuil exact est à fixer vendredi.

### Priorités (MoSCoW)

| Niveau | Contenu | Tâches |
|--------|---------|--------|
| **Indispensable pour mardi** | Socle (test gelé, configs, chargeur, PP partagée) ; LSTM et Transformer entraînés **au moins en graine 42** ; PP val < baseline unigramme ; échantillonneur ; validité/unicité/longueurs/composition ; lignes de référence ; test ouvert une fois ; rapport court | T48–T53, T55, T57–T59, T61, T64–T66 |
| **Important** | Graines 43 et 44 des deux modèles ; distance au train complète ; drapeau toxines ; relance reproductible d'une config ; tag `v0.3` | T54, T56, T60, T58, T67 |
| **Si temps (sinon reporté, avec raison écrite)** | Ablations ; analyse de mémorisation ; grille d'échantillonnage complète | T62, T63 |
| **Allégé** | Séminaire interne = répétition de 30 min lundi soir ; la **réunion de mardi** remplace la présentation | T68 |

**Grille d'échantillonnage réduite** (au lieu de 12 combinaisons) : température ∈ {0,8 ; 1,0 ; 1,2} × top-p ∈ {1,0 ; 0,9}, top-k désactivé (6 combinaisons), réglage sur val.

### Qui fait quoi, en gros

| Étudiant | Responsabilité globale de ces 4 jours | Livrables |
|----------|----------------------------------------|-----------|
| **Yassine Achour — Données** (**contraint : Ven < 16:00, Sam < 12:00**) | Fondations de données et traçabilité : tout ce qui doit exister **avant** les modèles, puis la traçabilité des générations | Configs YAML (T49), chargeur avec padding/masque (T50), stockage et re-filtrage des générations, drapeau toxines, script de lancement des 3 graines, extension du manifeste (T58), puis `PROVENANCE.md` |
| **Youssef Jarraya — Modèles** | Les deux modèles et leur entraînement : le chemin critique du projet | LSTM (T52), boucle d'entraînement (T53), Transformer (T55), entraînements 3 graines (T54, T56), courbes train/val |
| **Sofiene Ouanes — Évaluation** | Tout ce qui mesure : perplexité partagée, échantillonnage, métriques, références, distance au train | PP partagée (T51), échantillonneur (T57), métriques de génération (T59), références (T61), distance au train (T60) |
| **Moi — Chef / Coordinateur** | Protocole, gel du test, relectures, intégration, tableau final, **ouverture unique du test**, rapport, tag, préparation de la réunion | T48, T64, T65, T66, T67 + revue de **tous** les PR |

### Planning détaillé par jour

**Vendredi 30/10**

| Qui | Matin | Après-midi | Soir |
|-----|-------|------------|------|
| **Moi** | **T48 avant 11:00** : clôture S2, hash du test gelé, protocole écrit et commité ; ouvrir les branches | Relire et merger les PR de Yassine (T49, T50) | Relire T51/T52 ; vérifier la décision CPU/Colab |
| **Yassine** (jusqu'à 16:00) | T49 : configs et `set_seed` (jusqu'à ~11:30) | T50 : chargeur étendu, padding, masque, tests (**terminé et en PR avant 14:00**) ; puis début du script de lancement `run_all_seeds` et du squelette de T58 jusqu'à 16:00, avec une note de passage | **indisponible** |
| **Youssef** | T52 : LSTM + tests sur tenseurs synthétiques (sans attendre le chargeur) | Branchement au chargeur dès 14:00 ; début de T53 (boucle, early stopping, checkpoints) | T53 terminée ; **premier epoch chronométré** ; lancement du LSTM graine 42 |
| **Sofiene** | T51 : perplexité partagée ; reproduire le 19,49 de la baseline | T57 : échantillonneur, testé sur unigramme/bigramme | T57 terminé ; début T59 (métriques) sur val réelle |

*Jalons : test gelé (11:00) ; chargeur mergé (~15:00) ; LSTM 42 lancé (soir) ; décision CPU ou Colab.*

**Samedi 31/10**

| Qui | Matin | Après-midi | Soir |
|-----|-------|------------|------|
| **Moi** | Relire T53, T57, T59 ; préparer le script du tableau comparatif (T64) | Relire T58 ; vérifier les premières courbes LSTM | Rédiger le squelette du rapport |
| **Yassine** (jusqu'à 12:00) | T58 : stockage des générations, re-filtrage `filters.py`, drapeau toxines, extension du manifeste, **script des 3 graines** ; PR avec mode d'emploi **avant 12:00** | **indisponible** | **indisponible** |
| **Youssef** | Lecture des résultats LSTM 42 (PP val < 19,49 ?) ; **début T55 Transformer** | T55 terminé + test du masque causal ; lancement **Transformer 42** | Lancement LSTM 43/44 avec le script de Yassine |
| **Sofiene** | T59 terminé ; T61 : références (aléatoire, unigramme, bigramme, val réelle, val mélangée) | T60 : distance au train, d'abord val vs train (référence), puis première génération LSTM 42 | Premiers tableaux de métriques LSTM |

*Jalons : handoff de Yassine à 12:00 ; Transformer 42 lancé ; LSTM 43/44 en cours.*

**Dimanche 01/11**

| Qui | Matin | Après-midi | Soir |
|-----|-------|------------|------|
| **Moi** | Vérifier la reproductibilité bit-à-bit sur 1 config (CPU) | T64 : tableau comparatif et figures (échelle de PP) | **Signer les choix** (réglages d'échantillonnage, modèles) par commit/PR relue |
| **Yassine** (disponible `[à confirmer]`) | T58 sur les générations LSTM : re-filtrage, toxines, hashes | Idem pour le Transformer 42 ; `PROVENANCE.md` « Modèles et générations » | Relectures croisées (T60, T59) |
| **Youssef** | Transformer 42 : courbes ; lancement Transformer 43/44 | Suivi des graines ; ajustements **seulement si divergence** | T62 (si temps) |
| **Sofiene** | Grille d'échantillonnage réduite **sur val** pour le LSTM, puis le Transformer | T59/T60 sur les deux modèles ; sélection selon le critère écrit en T48 | Tableau des métriques par modèle ; choix proposés au coordinateur |

*Jalon : **réglages figés dimanche soir** ; aucun changement ensuite.*

**Lundi 02/11**

| Qui | Matin | Après-midi | Soir |
|-----|-------|------------|------|
| **Moi** | Dernières graines terminées ? Agrégation moyenne ± σ ; **T65 vers 12:00 : test ouvert une seule fois** | T66 : rapport, commentaire de l'écart val/test ; **tag `v0.3`** après T67 | **Répétition** de la réunion (30 min) |
| **Yassine** | Générations finales (réglage figé) pour toutes les graines : re-filtrage, toxines, hashes ; manifeste final | `PROVENANCE.md` final ; relance de contrôle depuis `data/raw/` | Répétition : sa partie (données, traçabilité) |
| **Youssef** | Vérifier les 6 runs (3 graines × 2 modèles), courbes finales | T66 : section « modèles » du rapport ; relance bit-à-bit d'**une** config | Répétition : sa partie (architectures, entraînement) |
| **Sofiene** | Métriques finales avec réglages figés pour toutes les graines (T59/T60) | T66 : inspection humaine de 20 séquences par modèle (hors toxines) ; section « évaluation » du rapport | Répétition : sa partie (échantillonnage, métriques, références) |

*Jalons : test ouvert (lundi ~12:00) ; tag `v0.3` ; support de réunion prêt avant minuit.*

**Mardi 03/11 — réunion avec M. Abdallah** : ordre de passage inchangé (socle et configs : Yassine ; architectures et entraînement : Youssef ; échantillonnage et évaluation : Sofiene ; synthèse, test verrouillé et limites : Moi).

### Plan de repli si on prend du retard

| Situation | Décision |
|-----------|----------|
| Le CPU est trop lent | Colab T4 dès samedi, checkpoints sur Drive ; **ne pas** réduire les règles (test, seeds, protocole) |
| Les graines 43/44 ne finissent pas avant lundi midi | Test ouvert avec la graine 42 + les graines terminées ; **dire clairement** « n = 1 » ou « n = 2 » dans le rapport et le tableau |
| Transformer instable (perte `NaN`) | Baisser le taux d'apprentissage, ajouter un warmup, réduire la taille ; documenter ; en dernier recours, ne présenter que le LSTM complet et le Transformer partiel |
| PP val ne dépasse pas la baseline unigramme | Résultat **à rapporter honnêtement** avec l'analyse (capacité, taux d'apprentissage, overfitting) ; ne pas bidouiller le test |
| Un livrable de Yassine n'est pas prêt à 16:00 vendredi ou 12:00 samedi | Sofiene ou Youssef reprennent la partie critique en priorité (T50 d'abord) ; T58 est reportée à dimanche/lundi, car elle ne bloque pas l'entraînement |
| Distance au train trop lente | Échantillonnage documenté (≥ 200 séquences par modèle, comme en S2), le calcul exhaustif reste en option |

### Ce qui sera réellement présentable mardi

Au minimum : l'échelle de perplexités (aléatoire / unigramme / bigramme / LSTM / Transformer), les courbes train/val, le tableau comparatif avec lignes de référence, le test ouvert une fois avec `TEST_LOCK.json`, la preuve de reproductibilité sur une config, et les limites. Les ablations, l'analyse de mémorisation et le conditionnement sont **annoncés comme suite**, pas comme acquis.

---

## 4. Liste des tâches (T48–T68) avec détails techniques

### A. Clôture de la semaine 2 et socle

**T48 — Clôturer S2, geler le test et figer le protocole (Moi)**
Vérifier que le tag `v0.2` est poussé et que la relance S2 reproduit les mêmes comptages. **Geler le test** : calculer le SHA-256 de `split_test.csv` et l'écrire dans `corpus_manifest.json` (champ `test_frozen_sha256`) avec la date. Figer par écrit le **protocole** : métriques, seuils, grille d'échantillonnage, règle « test ouvert une seule fois ». Ouvrir la branche d'intégration.

**À quoi ça sert en mots simples :** On décide des règles du jeu **avant** de voir les résultats. Un protocole écrit après coup se plie toujours en faveur du modèle.

**Validation :** tag `v0.2` sur GitHub ; hash du test dans le manifeste ; protocole commité et relu par un autre étudiant ; `pytest` vert (3.11 et 3.14).

**T49 — Configurations YAML versionnées et seeds (Yassine)**
Un fichier par modèle dans `configs/`, chargé par `peplib/config.py` dans une dataclass **validée** (clés manquantes ou inconnues = erreur). Contenu minimal :

```yaml
# configs/lstm_base.yaml  — points de départ, à justifier
model:
  type: lstm
  vocab_size: 23
  d_embed: 128
  d_hidden: 256
  n_layers: 2
  dropout: 0.2
  tie_weights: false
train:
  seed: 42            # graines de confirmation : 43, 44
  batch_size: 128
  max_epochs: 30      # budget égal Transformer
  optimizer: adamw
  lr: 1.0e-3
  weight_decay: 0.01
  grad_clip: 1.0
  early_stopping: {metric: val_ppl, patience: 3}
data:
  train: data/processed/split_train.csv
  val: data/processed/split_val.csv
  max_len: 50
sampling_grid:        # réglé sur val uniquement
  temperature: [0.8, 1.0, 1.2]
  top_k: [0, 10]      # 0 = désactivé
  top_p: [1.0, 0.9]
```

Une fonction `set_seed(seed)` unique : `random`, `numpy`, `torch.manual_seed`, `torch.cuda.manual_seed_all`, et `torch.use_deterministic_algorithms(True)` quand c'est possible.

**À quoi ça sert en mots simples :** Une expérience = une config + une seed. Si elles sont dans Git, n'importe qui relance exactement la même expérience, et on ne retrouve jamais « le bon réglage » de mémoire.

**Validation :** une config invalide échoue avec un message clair (test) ; la config est enregistrée **avec les logs** de chaque run ; deux appels de `set_seed(42)` donnent la même suite de tirages.

**T50 — Chargeur de données : padding, masque et reproductibilité (Yassine)**
Partir du `make_loader` existant (PR #20, DataLoader seedé, test bit-à-bit). Vérifier/étendre : lecture de `split_*.csv`, encodage `BOS + ids` (entrée) et `ids + EOS` (cible), **padding avec PAD=0** jusqu'à la longueur max du batch, retour de `lengths` et d'un masque booléen `pad_mask`. Option de **bucketing par longueur** (regrouper les séquences de longueur proche) pour accélérer le CPU, **désactivable** car il change l'ordre des batchs.

**À quoi ça sert en mots simples :** Les peptides n'ont pas la même longueur ; pour les empiler dans un tableau, on complète par des `PAD`. Il faut que ces `PAD` n'influencent **jamais** la perte, sinon on entraîne le modèle à prédire du vide.

**Validation :** entrée et cible décalées d'un cran (test sur 3 exemples écrits à la main) ; `PAD` absent des cibles utiles (masque) ; deux itérations avec la même seed → batchs identiques bit à bit ; `decode(encode(s)) == s` conservé ; **aucune séquence de val ou de test** n'est accessible depuis le chargeur d'entraînement.

**T51 — Entropie croisée et perplexité partagées (Sofiene)**
Une seule fonction `perplexity(logits, targets, mask)` dans `metrics_lm.py`, **utilisée aussi par `baselines.py`**. Convention : on compte les positions **résidus + EOS** (21 tokens possibles), on **exclut** `PAD` et `BOS` (jamais cibles). Retourne nats/token et PP = exp(nats/token).

**À quoi ça sert en mots simples :** Si le modèle et la baseline ne comptent pas les mêmes tokens, comparer leurs perplexités revient à comparer des mètres et des yards. Une seule fonction pour tout le monde supprime ce risque.

**Validation :** reproduit la perplexité de la baseline unigramme (`≈ 19,49` `[à vérifier]`) à 10⁻⁶ près en lui passant les mêmes probabilités ; un modèle uniforme sur 21 tokens donne exactement PP = 21 ; changer le nombre de `PAD` dans un batch ne change pas le résultat.

### B. Modèles

**T52 — LSTM char-level from scratch (Youssef)**
`nn.Module` : `Embedding(23, d_embed)` → `LSTM` (2 couches, dropout entre couches) → `Linear(d_hidden, 23)`. L'**interdiction des poids pré-entraînés** porte sur les poids, pas sur les briques : `nn.LSTM`, `nn.Embedding` et `nn.Linear` de PyTorch sont autorisés (**à confirmer avec l'encadrant ; sinon coder la cellule LSTM à la main**). Initialisation documentée, `forward(x, lengths)` → logits `[B, T, 23]`, méthode `step(token, state)` pour la génération incrémentale.

**À quoi ça sert en mots simples :** Le LSTM lit le peptide de gauche à droite en gardant une « mémoire » ; à chaque position il propose une probabilité pour chacun des 23 tokens.

**Validation :** formes de sortie correctes ; **causalité** : changer le token à la position t ne modifie pas les logits aux positions < t (test) ; **sur-apprentissage d'un mini-batch** : la perte tombe vers 0 en quelques centaines de pas (preuve que le gradient circule) ; nombre de paramètres rapporté.

**T53 — Boucle d'entraînement, early stopping, checkpoints (Youssef)**
`train.py` : perte = entropie croisée avec `ignore_index=PAD`, optimiseur `torch.optim.AdamW` `[à confirmer avec AGENTS.md : l'optimiseur est une brique, pas un poids]`, **gradient clipping** (1,0), évaluation sur val à chaque epoch, **early stopping** sur la perplexité val (patience 3), sauvegarde du **meilleur** checkpoint et du dernier (pour la reprise), log CSV par epoch (`epoch, train_loss, val_loss, val_ppl, lr, secondes`). Copie de la config et de la seed dans le dossier du run.

**À quoi ça sert en mots simples :** Entraîner trop longtemps « apprend par cœur » le train : la perte train baisse, la perte val remonte. L'early stopping s'arrête au meilleur moment, mesuré sur des peptides que le modèle n'a jamais vus.

**Validation :** relancer deux fois la même config sur CPU → **pertes identiques bit à bit** ; reprise après interruption → même trajectoire qu'un run non interrompu ; le test set n'est **pas lu** par ce script (test automatique : le chemin `split_test.csv` est refusé).

**T54 — Entraînement du LSTM, 3 graines (Youssef, Sofiene)**
Youssef lance la graine 42 ; Youssef ou Sofiene lancent 43 et 44 avec le script préparé par Yassine (en parallèle de T55). Chronométrer d'abord **un epoch CPU en J1** pour décider CPU ou Colab T4. Sauvegarder `checkpoints/lstm_seed{42,43,44}_best.pt` et hash dans le manifeste.

**À quoi ça sert en mots simples :** Un seul run peut être chanceux. Trois graines donnent une moyenne et un écart-type : on distingue « le modèle est meilleur » de « la graine a eu de la chance ».

**Validation :** les 3 runs convergent ; **PP val < baseline unigramme** (19,49 `[à vérifier]`) ; PP val rapportée en moyenne ± écart-type (ddof = 1) ; courbes train/val sauvegardées dans `results/week3/`.

**T55 — Transformer décodeur causal from scratch (Youssef)**
`nn.Module` : embedding des tokens + **encodage positionnel** (sinusoïdal ou appris, taille max 52 = BOS + 50 + EOS ; choix justifié), `N` blocs (**pre-LayerNorm** → auto-attention **masquée causalement** multi-têtes → MLP GELU), tête linéaire. Point de départ : `d_model = 128, 4 têtes, 4 couches, d_ff = 512, dropout 0,1–0,2`, **nombre de paramètres comparable au LSTM** (écart ≤ 20 %, à rapporter). L'attention peut être codée à la main (produit scalaire, masque, softmax) : recommandé pour la défense orale.

**À quoi ça sert en mots simples :** Au lieu d'une mémoire qui se transmet de lettre en lettre, chaque position « regarde » directement toutes les lettres précédentes et décide lesquelles comptent. Le masque causal empêche de regarder le futur (sinon la tâche serait triviale).

**Validation :** **test du masque causal** : modifier la lettre en position t ne change aucun logit aux positions < t ; sur-apprentissage d'un mini-batch ; le `PAD` n'influence pas les positions valides (test avec et sans padding supplémentaire) ; formes et nombre de paramètres.

**T56 — Entraînement du Transformer à budget égal (Youssef, Sofiene)**
Même données, même `max_epochs`, même `batch_size`, même critère d'arrêt que le LSTM ; seul l'optimiseur change si nécessaire (warmup court des premières centaines de pas, puis taux constant ou cosinus) et doit être **documenté**. 3 graines. **Budget égal** = mêmes epochs maximales + paramètres comparables + temps de calcul rapporté.

**À quoi ça sert en mots simples :** Comparer deux modèles n'a de sens que s'ils ont eu les mêmes chances. Sinon on compare des budgets, pas des architectures (le mois 8 exige ce benchmark à budget égal).

**Validation :** PP val moyenne ± écart-type sur 3 graines ; comparaison avec LSTM **dans un même tableau** ; tout écart de configuration entre les deux modèles est listé ; pas de divergence (perte `NaN`) : sinon baisse du taux d'apprentissage, documentée.

### C. Génération

**T57 — Échantillonneur autorégressif (Sofiene)**
`sampling.py` : partir de `BOS`, tirer jusqu'à `EOS` ou 50 résidus. À **chaque pas** : logits de `PAD` et `BOS` forcés à **−∞** (jamais tirables) ; division par la **température** ; **top-k** (garder les k plus probables) ; **top-p** (plus petit ensemble de probabilité cumulée ≥ p) ; tirage avec un générateur seedé. Règles de longueur **explicites** : 
- *mode libre* (référence) : le modèle choisit quand s'arrêter ; une séquence < 10 ou sans EOS à 50 est **comptée invalide** (pas de correction silencieuse) ;
- *mode contraint* (variante rapportée à part) : `EOS` masqué avant 10 résidus, forcé à 50.

L'échantillonneur est développé et **testé d'abord sur les baselines unigramme/bigramme** (qui sont aussi des « modèles » avec `p(token | contexte)`), avant l'arrivée des vrais modèles.

**À quoi ça sert en mots simples :** La température règle l'audace : basse = prudente et répétitive, haute = variée mais plus souvent invalide. Top-k/top-p coupent la queue des choix absurdes. On mesure l'effet de chaque réglage au lieu de le deviner.

**Validation :** `PAD`/`BOS` **jamais** générés sur 10 000 tirages ; température → 0 donne la séquence gloutonne ; même seed → mêmes séquences bit à bit ; top-k = 1 équivaut au glouton ; top-p = 1 et top-k = 0 reproduisent la distribution du modèle (test statistique sur la baseline unigramme : fréquences générées ≈ fréquences théoriques).

**T58 — Stockage des générations, re-filtrage et drapeau toxines (Yassine)**
Pour chaque (modèle, graine d'entraînement, réglage d'échantillonnage) : N = 1 000 séquences dans `data/generated/<run>.csv` (`id, sequence, length, valid, model, seed_train, seed_sampling, temperature, top_k, top_p`). **Re-filtrage** par `filters.py` (alphabet canonique, 10–50) → colonne `valid`. **Drapeau toxines** : relancer la logique de `toxin_flags.json` sur chaque génération, colonne `toxin_flag`. Hash SHA-256 de chaque CSV dans le manifeste.

**À quoi ça sert en mots simples :** Chaque séquence générée est un résultat scientifique : elle doit pouvoir être régénérée et auditée. Et un générateur de peptides peut, par accident, ressembler à une toxine : on le vérifie systématiquement.

**Validation :** régénération avec les mêmes seeds → CSV identique octet pour octet ; toute séquence signalée toxine est **exclue des figures, exemples et exports** ; compteurs `valides + invalides = N`.

### D. Évaluation

**T59 — Métriques de génération (Sofiene)**
Pour chaque jeu de N = 1 000 séquences :
- **Validité** : % passant `filters.py` (alphabet canonique, 10–50, EOS présent).
- **Unicité** : % de séquences distinctes ; **entropie des 3-mers** des sorties.
- **Longueurs** : histogramme vs val ; statistique **KS** (Kolmogorov–Smirnov, quelques lignes de numpy).
- **Composition** : divergence de **Jensen–Shannon** des fréquences d'acides aminés vs val (réutilise le code S1).
- **Descripteurs** (charge à pH 7, pI, GRAVY, moment hydrophobe) : médiane, IQR et KS vs val.

**À quoi ça sert en mots simples :** Une bonne perplexité ne garantit pas de bons peptides : un modèle peut produire toujours la même séquence plausible. Ces mesures répondent à : « valides ? variés ? qui ressemblent aux vrais en longueur, composition et propriétés ? ».

**Validation :** KS = 0 et JS = 0 quand on compare val avec elle-même ; les métriques calculées sur **séquences de val réelles** donnent validité 100 % ; formules testées sur de petits exemples à la main.

**T60 — Distance au train : pas de recopie (Sofiene)**
Pour chaque génération valide, trouver son **plus proche voisin du train** : candidats via l'index 3-mers (comme `cluster.py`), identité = Levenshtein normalisée (la **même mesure que le clustering à 90 %**), alignement Smith–Waterman sur un sous-échantillon (≥ 200 séquences) en vérification. Rapporter : distribution de l'identité max, **% de générations avec identité ≥ 0,90** (recopie probable), % de copies exactes. **Référence** : val vs train (< 0,90 par construction du split).

**À quoi ça sert en mots simples :** Un modèle peut obtenir une excellente perplexité en **mémorisant** le train. S'il recopie des peptides existants, il ne « génère » rien. On mesure à quel point les sorties sont proches des peptides d'entraînement, avec la même règle des 90 % que pour le split.

**Validation :** copie exacte d'une séquence du train → identité = 1,0 (test) ; val vs train **< 0,90** reproduit le chiffre S2 ; **seuil de lecture** : objectif 0 % de générations ≥ 0,90 ; tolérance ≤ 1 % `[à valider avec l'encadrant]` ; au-delà → signe de mémorisation, à analyser (T63) ; temps de calcul chronométré (échantillonnage documenté, comme en S2).

**T61 — Lignes de référence (Sofiene)**
Les chiffres d'un modèle n'ont de sens que **comparés**. Mesurer, avec **exactement les mêmes fonctions**, quatre références :
1. générateur **aléatoire uniforme** (borne basse) ;
2. générateur **unigramme** et **bigramme** (S2) ;
3. **val réelle** (borne haute : validité 100 %, distance au train < 0,90) ;
4. **val mélangée** (lettres permutées, seed 42 : même composition, ordre détruit).

**À quoi ça sert en mots simples :** Le moment hydrophobe dépend de l'**ordre** des lettres : le modèle doit s'approcher des vrais peptides et s'éloigner des mélangés. Sans ces références, « 93 % de validité » ne dit pas si c'est bon.

**Validation :** les 5 lignes (réel, mélangé, aléatoire, unigramme, bigramme) sont dans le **même tableau** que les modèles ; ordre attendu des descripteurs sensibles à l'ordre : réel ≈ modèle ≫ mélangé ; charge **identique** entre réel et mélangé (composition seule).

**T62 — Ablations (si temps) (Youssef)**
Faire varier **un facteur à la fois** (profondeur 1/2/3 pour le LSTM, 2/4/6 pour le Transformer ; dimension ; dropout 0,1/0,2/0,3) avec **3 graines**, PP val uniquement. Un petit tableau + une figure.

**À quoi ça sert en mots simples :** Justifier un choix d'architecture devant le professeur, avec des chiffres plutôt que « on a pris ça parce que c'est classique ».

**Validation :** un seul facteur change par ligne ; moyenne ± écart-type ; les écarts plus petits que l'écart-type sont déclarés **non significatifs**.

**T63 — Analyse de mémorisation (si temps) (Sofiene)**
Écart perte train vs val par epoch ; identité au plus proche voisin vs **température** ; 10 générations les plus proches du train, affichées et commentées.

**À quoi ça sert en mots simples :** Si la perplexité s'améliore mais que les sorties se rapprochent du train, le modèle apprend par cœur au lieu de généraliser.

**Validation :** une figure (identité max vs température) ; conclusion écrite, même négative.

### E. Comparaison, test et clôture

**T64 — Tableau comparatif et figures (Moi)**
Un tableau unique (lignes : références T61 + LSTM + Transformer ; colonnes : PP val ± σ, validité, unicité, KS longueur, JS composition, KS descripteurs, % ≥ 0,90 au train, paramètres, temps). Figure : **échelle de perplexités** (aléatoire / unigramme / bigramme / LSTM / Transformer), à la même échelle que `ecart_baselines.png`. Choix du **réglage d'échantillonnage final** par modèle, **sur val**, avec le critère écrit à l'avance en T48.

**À quoi ça sert en mots simples :** Une page qui permet à M. Abdallah de voir en 30 secondes ce qui marche, ce qui ne marche pas, et à quel prix.

**Validation :** toutes les cellules régénérables par un script ; seeds écrites ; aucune valeur recopiée à la main.

**T65 — Évaluation sur test, UNE SEULE FOIS (Moi)**
Avant : réglages et choix de modèle **signés** (commit + PR relue). Script `eval_test.py` : vérifie le hash de `split_test.csv` (= `test_frozen_sha256`), **refuse de s'exécuter si `results/week3/TEST_LOCK.json` existe**, calcule PP test (+ métriques de génération pour le réglage figé), écrit `TEST_LOCK.json` (date, commit, hash, résultats). Aucun réglage n'est modifié après.

**À quoi ça sert en mots simples :** C'est l'examen final : on le passe une fois. Si le résultat déçoit, on le rapporte tel quel et on explique, on ne « corrige » pas jusqu'à ce que ça passe.

**Validation :** `TEST_LOCK.json` présent, horodaté, avec le commit ; hash du test conforme ; seconde exécution refusée (test) ; écart PP val → PP test commenté (un écart important = signal à discuter).

**T66 — Inspection humaine et rapport (Tous)**
Lecture, par le binôme, de **20 séquences** par modèle (réglage final, hors toxines) : régularités, répétitions, motifs absurdes (ex. poly-Leu), motifs plausibles (alternance hydrophobe/cationique). Rapport `Docs/Semaine3_rapport.md` : objectif, méthode, tableau, figures, limites, **options exclues et pourquoi**.

**À quoi ça sert en mots simples :** Les métriques disent « ça ressemble statistiquement » ; les yeux détectent ce que les métriques ratent. Mais c'est qualitatif : on ne conclut rien sur l'activité biologique.

**Validation :** limites déclarées (aucune mesure d'activité ; pas de LLM/RAG/LoRA ; pas de conditionnement) ; **aucune phrase du type « ce peptide est antimicrobien »** ; chiffres du rapport = chiffres du tableau T64.

**T67 — Relance complète et tag `v0.3` (Tous)**
Comme T30/T46 : depuis `data/raw/`, régénération de `data/processed/`, puis ré-entraînement CPU d'**au moins une config** et régénération d'un jeu de générations ; mêmes comptages, mêmes familles, **mêmes pertes bit à bit**, mêmes séquences générées pour une même seed. Tag `v0.3`. **Moi** merge vers `main`.

**Validation :** les 4 étudiants obtiennent les mêmes valeurs ; manifeste identique octet pour octet (hash des configs, seeds, générations ; hashes de checkpoints **à égalité sur CPU**, sinon tolérance documentée pour GPU).

**T68 — Séminaire interne et présentation (Tous)** — ordre de passage : socle, configs et chargeur (Yassine), architectures et entraînement (Youssef), échantillonnage et évaluation (Sofiene), synthèse, test verrouillé et limites (Moi).

---

## 5. Formules à implémenter

| Quantité | Définition exacte |
|----------|-------------------|
| Entrée / cible (teacher forcing) | entrée = `BOS, x₁…x_L` ; cible = `x₁…x_L, EOS` ; positions `PAD` ignorées |
| Entropie croisée | `H = −(1/N) Σ log p(y_t \| x_<t)` sur les N positions **valides** (résidus + EOS), en nats |
| Perplexité | `PP = exp(H)` ; modèle uniforme sur 21 tokens ⇒ PP = 21 |
| Softmax avec température | `p_i = exp(z_i/τ) / Σ_j exp(z_j/τ)` ; `τ → 0` : glouton ; `τ = 1` : distribution du modèle |
| Top-k | ne garder que les k logits les plus grands (reste = −∞), renormaliser |
| Top-p | trier par probabilité décroissante, garder le plus petit préfixe de masse cumulée ≥ p, renormaliser |
| Masque PAD/BOS à l'échantillonnage | `logits[PAD] = logits[BOS] = −∞` à **chaque** pas |
| Attention causale | `softmax(QKᵀ/√d_k + M) V`, `M_ij = −∞` si `j > i` |
| Gradient clipping | si `‖g‖₂ > c`, `g ← g · c/‖g‖₂` (c = 1,0) |
| Validité | `% séquences passant filters.py` (alphabet canonique, 10–50, EOS présent) |
| Unicité | `|séquences distinctes| / N` |
| Identité Levenshtein | `1 − d_Lev(a,b) / max(|a|,|b|)` (mesure du clustering à 90 %, S2) |
| Distance au train | `id_max(g) = max_{t ∈ candidats 3-mers} identité(g, t)` ; recopie probable si ≥ 0,90 |
| KS (deux échantillons) | `D = max_x \|F₁(x) − F₂(x)\|` (fonctions de répartition empiriques) |
| Jensen–Shannon | `JS(P,Q) = ½ KL(P‖M) + ½ KL(Q‖M)`, `M = (P+Q)/2` |
| Moyenne ± écart-type (3 graines) | écart-type d'échantillon (`ddof = 1`) |
| Seeds | entraînement : 42 (référence), 43, 44 ; échantillonnage : seed documentée par run |

---

## 6. Tests obligatoires

| Test | Résultat attendu |
|------|------------------|
| Config | clé manquante/inconnue ⇒ erreur ; config copiée dans le dossier du run |
| `set_seed` | deux appels ⇒ mêmes tirages (`random`, `numpy`, `torch`) |
| Chargeur | entrée/cible décalées ; `PAD` ignoré ; même seed ⇒ batchs identiques bit à bit ; aucun accès à val/test depuis le train |
| Round-trip (régression S2) | `decode(encode(s)) == s` ; `len(VOCAB) == 23` |
| Perplexité | modèle uniforme 21 tokens ⇒ 21 ; PP baseline unigramme reproduite à 10⁻⁶ ; invariance au nombre de `PAD` |
| LSTM / Transformer — causalité | changer le token en t ne modifie aucun logit en < t |
| Sur-apprentissage d'un mini-batch | perte → ~0 en quelques centaines de pas (les deux modèles) |
| PAD sans effet | ajouter du padding ne change pas les logits des positions valides |
| Entraînement bit-à-bit | 2 runs CPU, même config/seed ⇒ pertes identiques |
| Reprise | run interrompu puis repris ⇒ même trajectoire |
| Test verrouillé | `train.py` refuse `split_test.csv` ; `eval_test.py` refuse une 2ᵉ exécution et un hash ≠ `test_frozen_sha256` |
| Échantillonneur | `PAD`/`BOS` jamais tirés (10 000 tirages) ; mêmes seeds ⇒ mêmes séquences ; top-k = 1 ⇒ glouton ; fréquences ≈ théoriques sur unigramme |
| Métriques de génération | val vs val ⇒ KS = 0, JS = 0 ; val réelle ⇒ validité 100 % |
| Distance au train | copie exacte ⇒ 1,0 ; val vs train < 0,90 |
| Références | charge réel = charge mélangé ; moment hydrophobe mélangé ≪ réel |
| Génération reproductible | mêmes seeds ⇒ CSV identique octet pour octet |
| Biosécurité | séquences `toxin_flag` absentes des exports et figures |
| Manifeste | régénération ⇒ manifeste identique ; hash ≠ ⇒ test échoue sans bump de version |

---

## 7. Pièges à éviter

- **PP du modèle ≠ PP des baselines** : mêmes tokens comptés (résidus + EOS, sans PAD ni BOS), même fonction (T51). C'est le piège n° 1 de la semaine.
- **PAD dans la perte** : sans `ignore_index`/masque, le modèle apprend à prédire du padding et la perplexité devient trop belle. Vérifier par le test d'invariance au padding.
- **Masque causal faux** : le modèle « triche » en regardant le futur, la perte train chute anormalement, mais la génération est absurde. Le test de causalité n'est pas optionnel.
- **Test set utilisé pour choisir** : même de façon indirecte (choisir la température qui « marche mieux sur test »). Tout réglage sur val ; test ouvert une fois (T65).
- **Mémorisation prise pour du succès** : PP val excellente mais générations proches du train. Toujours regarder T60.
- **Effondrement de modes** : toutes les séquences se ressemblent (riches en Leu/Lys). La **unicité** et l'**entropie des 3-mers** le détectent ; la validité seule ne le voit pas.
- **Validité « corrigée » en douce** : tronquer ou compléter les séquences pour qu'elles passent le filtre masque un défaut du modèle. Mode libre = référence ; mode contraint rapporté à part.
- **Reproductibilité GPU** : `nn.LSTM` cuDNN et certains noyaux sont **non déterministes**. La preuve bit-à-bit se fait **sur CPU** ; en GPU, tolérance documentée.
- **Budget inégal** : un Transformer entraîné plus longtemps ou avec plus de paramètres ne prouve rien. Écarts de config listés (T56).
- **Comparer à une seule graine** : un écart plus petit que l'écart-type de 3 graines n'est pas un résultat.
- **Réglage d'hyperparamètres sans trace** : toute valeur finale doit être dans un YAML versionné, avec le log du run qui la justifie.
- **Modifier V = 23 en S3** : la diffusion discrète masquée demandera un token `[MASK]` (V = 24) plus tard. Ne **pas** l'ajouter maintenant ; ce sera une version de tokenisation à part, avec test de non-régression.
- **Notebook contre module** : la logique vit dans `peplib` et ses tests ; le notebook appelle seulement.
- **Visualisation 3D interdite** (rappel) : matplotlib 2D uniquement.
- **UTF-8** : fichiers sans BOM, sans mojibake (vérifier chaque écriture, comme en S1–S2).
- **Seeds** : 42 en référence, 43 et 44 en confirmation ; seed d'échantillonnage écrite dans chaque CSV généré.

---

## 8. Biosécurité, provenance et licences

- Un modèle génératif entraîné sur des AMP peut produire des séquences ressemblant à des **toxines ou hémolysines** présentes dans les bases amont. Le drapeau toxines (T58) est appliqué à **chaque** génération ; les séquences signalées ne figurent ni dans les figures, ni dans les exemples, ni dans les exports.
- **Les séquences générées ne sont pas distribuées.** `data/generated/` n'est pas partagé hors du groupe sans décision écrite du coordinateur et de l'encadrant ; `PROVENANCE.md` consigne ce choix.
- **Aucune affirmation d'activité biologique.** Charge, pI, GRAVY et moment hydrophobe sont des **proxys physico-chimiques**, pas des mesures d'activité. Toute phrase « ce peptide est antimicrobien » est interdite sans validation en laboratoire qualifié ; c'est une limite déclarée du rapport.
- **Licences** : le corpus d'entraînement reste soumis aux décisions de T32 (APD3, AMPChallenge, 13 bases amont). Un modèle entraîné sur un corpus dont la redistribution est restreinte n'est **pas** automatiquement redistribuable : les checkpoints restent locaux (`.gitignore`), le manifeste ne contient que leurs hashes.
- `PROVENANCE.md` : ajouter un bloc « Modèles et générations » (configs, seeds, hashes, date, commit, statut de redistribution).

---

## 9. Critères de réussite de la semaine

**Barrières de performance (mesurées sur val ; test ouvert une seule fois) :**

| Critère | Seuil |
|---------|-------|
| PP val LSTM et Transformer | **< 19,49** (baseline unigramme `[à vérifier]`) ; objectif secondaire : ≤ bigramme `[valeur du manifeste S2]` |
| Validité (mode libre) | > 95 % `[objectif ; à discuter si hors d'atteinte]` |
| Unicité (N = 1 000) | > 90 % |
| Générations avec identité ≥ 0,90 au train | 0 % visé ; ≤ 1 % toléré, documenté |
| Résultats | moyenne ± écart-type sur 3 graines |

**Critères de livraison :**
- Les tâches T48–T68 sont cochées, ou les restantes (T62, T63) reportées avec une raison écrite.
- Chaque résultat est régénérable depuis une config versionnée + une seed.
- Test ouvert une seule fois (`TEST_LOCK.json`), écart val/test commenté.
- Tableau comparatif avec références (aléatoire, unigramme, bigramme, val réelle, mélangée).
- Rapport avec limites et options exclues ; manifeste et `PROVENANCE.md` à jour ; tag `v0.3` poussé.
- **Un seuil non atteint n'est pas un échec du projet** : c'est un résultat à rapporter et analyser honnêtement.

---

## 10. Ce que le superchef (moi) doit présenter à M. Abdallah

M. Abdallah doit voir des **visuels et des chiffres de garde-fou**, pas du code.

### Visuels à présenter
- **Échelle de perplexités** : aléatoire / unigramme / bigramme / LSTM / Transformer (mêmes échelles que `ecart_baselines.png`)
- **Courbes train/val** (une par modèle) montrant l'early stopping
- **Tableau comparatif** : validité, unicité, distance au train, descripteurs, avec lignes de référence
- **Distribution de l'identité au plus proche voisin du train** : générations vs val
- **Histogrammes de longueurs** et **descripteurs** : réel / mélangé / modèles
- *(si temps)* ablations et figure identité vs température

### Documents à présenter
- Protocole d'évaluation figé (T48) et `TEST_LOCK.json`
- Configs YAML des deux modèles
- Rapport `Semaine3_rapport.md` (limites et options exclues incluses)
- `PROVENANCE.md` et `corpus_manifest.json` mis à jour

### Preuve de reproductibilité à montrer
- Tags `v0.1`, `v0.2` et `v0.3` poussés
- Les 4 étudiants obtiennent les mêmes pertes (CPU) et les mêmes séquences générées pour une même seed

**On doit savoir expliquer chacun de ces éléments sans ouvrir le code**, surtout : (i) pourquoi LLM/RAG/LoRA sont exclus (contrainte *et* fuite), (ii) pourquoi la perplexité compte les mêmes tokens pour tous, (iii) pourquoi le test n'est ouvert qu'une fois, (iv) pourquoi une bonne perplexité ne suffit pas (mémorisation, effondrement de modes), (v) ce que **ne** démontrent **pas** les proxys physico-chimiques.

---

## 11. Commandes de démarrage rapide (interface cible)

Ouvrir un terminal **à la racine du projet** (`C:\Users\DELL\Projects\Projet_Peptide`) :

```bash
# 1. Environnement (déjà créé en S1–S2)
.venv\Scripts\Activate.ps1

# 2. Tests (2 interpréteurs)
pytest
python -m pytest

# 3. Régénérer le corpus traité (NE PAS re-télécharger)
python -m peplib.filters

# 4. Entraîner (interface cible ; adapter aux noms réels des modules)
python -m peplib.train --config configs/lstm_base.yaml --seed 42
python -m peplib.train --config configs/transformer_base.yaml --seed 42

# 5. Générer et évaluer (réglage d'échantillonnage choisi sur val)
python -m peplib.sampling --run lstm_seed42 --n 1000 --temperature 1.0 --seed 7
python -m peplib.eval_gen --run lstm_seed42

# 6. Test final — UNE SEULE FOIS, par le coordinateur
python -m peplib.eval_test --lstm lstm_seed42 --transformer transformer_seed42
```

Résultat attendu : tests verts ; pertes et séquences identiques d'une relance à l'autre (CPU).

---

## 12. Gabarits

**Journal d'expérience (`results/week3/runs.csv`, une ligne par run) :**

```csv
run_id,model,config_hash,seed_train,params,epochs_done,best_epoch,val_ppl,train_ppl,seconds,device,commit,checkpoint_sha256
lstm_seed42,lstm,<hash>,42,<n>,<n>,<n>,<x>,<x>,<s>,cpu,<sha>,<hash>
```

**Bloc `PROVENANCE.md` — modèles et générations :**

```markdown
# Modèles et générations — S3

- **Corpus d'entraînement** : corpus_v0.2 (split_train), hash : [SHA-256]
- **Modèles** : LSTM (config [hash]), Transformer (config [hash]), 3 graines chacun (42, 43, 44)
- **Poids pré-entraînés** : AUCUN (entraînement from scratch)
- **Appareil / bibliothèques** : [CPU ou GPU], PyTorch [version], Python [3.11 / 3.14]
- **Checkpoints** : locaux, non commités ; hashes dans corpus_manifest.json
- **Générations** : data/generated/ ; N = 1 000 par run ; non distribuées
- **Drapeau toxines** : [nb signalées / nb total], séquences exclues des figures
- **Test set** : ouvert une fois le [date], commit [sha], TEST_LOCK.json
- **Limites** : proxys physico-chimiques uniquement ; aucune validation d'activité biologique
```

---

## 13. Passerelle vers la suite (hors périmètre S3)

Rien de ce qui suit n'est codé cette semaine ; ces points **orientent les choix de S3**.

| Étape | Contenu | Ce que S3 prépare |
|-------|---------|-------------------|
| Suivante | **Diffusion discrète masquée** (absorbing state, type D3PM/MDLM) : Transformer **bidirectionnel**, embedding du pas de temps, token `[MASK]` | chargeur, boucle d'entraînement, évaluation et lignes de référence réutilisables ; V = 24 en version de tokenisation à part |
| Ensuite | **Conditionnement** : (1) token de condition en tête, (2) modulation adaLN/FiLM, (3) classifier-free guidance | conditions candidates : longueur, charge, classes de descripteurs (toujours disponibles) ; activité/organisme **seulement si fiables** (en-têtes AMPChallenge, métadonnées APD3 à auditer) |
| Plus tard | VAE | benchmark à budget égal des 4 familles (mois 8) |

**Principe : ne jamais ajouter deux nouveautés à la fois.** Si la diffusion conditionnée échoue, on doit savoir si le problème vient de la diffusion (étape non conditionnée d'abord) ou du conditionnement. Évaluations spécifiques à prévoir : **fidélité à la condition** (écart entre valeur demandée et valeur réelle), **diversité à condition fixée**, absence de recopie du train.
