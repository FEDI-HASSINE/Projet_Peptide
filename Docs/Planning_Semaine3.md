# Planning Semaine 3 — 4 jours (vendredi 30/10 – lundi 02/11/2026), réunion mardi 03/11

**Groupe :** Moi (Chef/Coordinateur), Yassine Achour (Données), Youssef Jarraya (Modèles), Sofiene Ouanes (Évaluation)
**Objectif :** LSTM et Transformer entraînés from scratch (aucun poids pré-entraîné), générations évaluées contre des références, test ouvert **une seule fois**, tag `v0.3`, prêts pour la réunion avec M. Abdallah le mardi.

> **Contraintes de disponibilité (Yassine) :** vendredi, plus rien à partir de **16:00** ; samedi, plus rien à partir de **12:00** ; dimanche et lundi : disponible `[à confirmer]`.
> Ses livrables critiques sont donc **terminés avant ses coupures**, et il ne lance aucun entraînement après.
> Détails techniques de chaque tâche : `Semaine3_plan.md` (T48–T68).

---

## Vendredi 30/10 — Fondation : protocole, chargeur, premier modèle

| Moment | Qui | Tâche | Livrable |
|--------|-----|-------|----------|
| Matin (avant 11:00) | **Moi** | T48 — Clôture S2 (tag `v0.2` vérifié), **hash du test gelé**, protocole d'évaluation écrit et commité | `test_frozen_sha256` dans `corpus_manifest.json`, protocole relu par un autre étudiant |
| Matin (Meet 15 min) | **Tous** | Lancement : on relit le protocole, les seuils et les dépendances | Tout le monde sait quoi livrer et pour quand |
| Matin | **Yassine** | T49 — Configs YAML versionnées + `set_seed` + validation des clés | `configs/lstm_base.yaml`, `configs/transformer_base.yaml`, tests de config |
| Matin | **Youssef** | T52 — LSTM `nn.Module` + tests sur tenseurs synthétiques (sans attendre le chargeur) | `peplib/models/lstm.py`, tests de forme et de causalité |
| Matin | **Sofiene** | T51 — Perplexité partagée (tokens = résidus + EOS), reproduction du 19,49 de la baseline | `peplib/metrics_lm.py`, test PP uniforme = 21 |
| Après-midi (**avant 14:00**) | **Yassine** | T50 — Chargeur : `make_loader` étendu (padding PAD, masque, bucketing optionnel) + tests bit-à-bit | `peplib/dataset.py`, PR ouverte **avant 14:00** |
| Après-midi (jusqu'à 16:00) | **Yassine** | Début du script `run_all_seeds` (seeds 42/43/44) et du squelette de T58 ; **note de passage** écrite | Note de passage dans la PR (ce qui reste, comment lancer) |
| Après-midi | **Moi** | Relecture et merge des PR de Yassine (T49, T50) | `main` à jour avec config + chargeur |
| Après-midi | **Youssef** | Branchement du LSTM au chargeur dès 14:00 ; début T53 (boucle, early stopping, checkpoints) | Boucle d'entraînement en cours |
| Après-midi | **Sofiene** | T57 — Échantillonneur (température, top-k, top-p, PAD/BOS à −∞), **testé sur les baselines unigramme/bigramme** | `peplib/sampling.py`, tests |
| Soir | **Youssef** | T53 terminée, **premier epoch chronométré** (CPU ou Colab T4 ?), lancement du **LSTM graine 42** | Run 42 lancé, décision CPU/Colab consignée |
| Soir | **Sofiene** | Début T59 — Métriques de génération testées sur val réelle | Validité, unicité, longueurs sur val |

**Ce soir :** test gelé, configs et chargeur mergés dans `main`, LSTM 42 lancé, décision CPU/Colab notée.

---

## Samedi 31/10 — Transformer et évaluation (Yassine s'arrête à 12:00)

| Moment | Qui | Tâche | Livrable |
|--------|-----|-------|----------|
| Matin (jusqu'à 12:00) | **Yassine** | T58 — Stockage des générations (`data/generated/`), re-filtrage `filters.py`, drapeau toxines, extension du manifeste | `data/generated/` + colonnes `valid`, `toxin_flag` |
| Matin (jusqu'à 12:00) | **Yassine** | Script de lancement des 3 graines + mini-mode d'emploi | PR avec mode d'emploi, **ouverte avant 12:00** |
| Matin (Meet 15 min, ~11:30) | **Tous** | **Passation de Yassine** : il explique ses scripts, les autres posent leurs questions avant sa coupure | Tout le monde sait lancer les graines et le stockage |
| Matin | **Youssef** | Lecture du LSTM 42 (PP val < 19,49 ?) ; début T55 — Transformer causal | Courbes train/val LSTM 42 |
| Matin | **Sofiene** | T59 terminé ; T61 — Lignes de référence (aléatoire, unigramme, bigramme, val réelle, val mélangée) | `results/week3/` : tableau de référence |
| Matin | **Moi** | Relecture T53, T57, T59 ; squelette du script de tableau comparatif (T64) | PR relues, squelette T64 |
| Après-midi | **Youssef** | T55 terminé + **test du masque causal** ; lancement du **Transformer graine 42** | `peplib/models/transformer.py`, run 42 lancé |
| Après-midi | **Sofiene** | T60 — Distance au train : d'abord val vs train (référence), puis première génération du LSTM 42 | `% ≥ 0,90` au train, référence val < 0,90 |
| Après-midi | **Moi** | Relecture de T58 ; contrôle des premières courbes | Validation ou demande de corrections |
| Soir | **Youssef** | Lancement des graines **43 et 44 du LSTM** (script de Yassine) | Runs 43/44 en cours |
| Soir | **Sofiene** | Premiers tableaux de métriques LSTM 42 | Tableau (validité, unicité, JS, KS, distance au train) |
| Soir | **Moi** | Squelette du rapport | `Docs/Semaine3_rapport.md` (plan) |

**Ce soir :** Transformer 42 lancé, LSTM 43/44 en cours, évaluation LSTM 42 disponible.

---

## Dimanche 01/11 — Générations, grille d'échantillonnage, gel des réglages

| Moment | Qui | Tâche | Livrable |
|--------|-----|-------|----------|
| Matin | **Youssef** | Suivi Transformer 42 ; lancement des graines **43 et 44 du Transformer** | Courbes Transformer 42, runs 43/44 lancés |
| Matin | **Sofiene** | Grille d'échantillonnage réduite **sur val** pour le LSTM (T ∈ {0,8 ; 1,0 ; 1,2} × top-p ∈ {1,0 ; 0,9}) | Tableau de la grille LSTM |
| Matin | **Yassine** | T58 sur les générations LSTM : re-filtrage, toxines, hashes | CSV de générations + hashes dans le manifeste |
| Matin | **Moi** | Vérification de la reproductibilité bit-à-bit sur 1 config (CPU) | 2 runs → mêmes pertes |
| Après-midi | **Sofiene** | Même grille pour le Transformer ; T59/T60 sur les deux modèles | Tableaux de métriques des 2 modèles |
| Après-midi | **Yassine** | T58 sur le Transformer 42 ; `PROVENANCE.md` « Modèles et générations » | `PROVENANCE.md` mis à jour |
| Après-midi | **Moi** | T64 — Tableau comparatif et figures (échelle de perplexités) | `results/week3/` : tableau + figures |
| Après-midi | **Youssef** | Ajustements **uniquement si divergence** (perte `NaN`, instabilité) ; T62 si temps | Runs stables |
| Soir (Meet 20 min) | **Tous** | **Gel des réglages** : choix d'échantillonnage et de modèles selon le critère écrit en T48 | Commit/PR « choix signés », relu par un autre étudiant |
| Soir | **Yassine / Youssef** | Relectures croisées (T60, T59, T53) | PR relues |

> **Règle :** après le gel de dimanche soir, **aucun réglage ne change**. Le test est ouvert lundi, sans retouche.

---

## Lundi 02/11 — Test unique, rapport, tag

| Moment | Qui | Tâche | Livrable |
|--------|-----|-------|----------|
| Matin | **Youssef** | Vérification des 6 runs (3 graines × 2 modèles), courbes finales | Tableau des runs `runs.csv` |
| Matin | **Sofiene** | Métriques finales avec réglages figés pour toutes les graines (T59/T60) | Moyennes ± σ des métriques |
| Matin | **Yassine** | Générations finales pour toutes les graines : re-filtrage, toxines, hashes ; manifeste final | `corpus_manifest.json` final |
| Matin | **Moi** | Agrégation moyenne ± σ ; dernières PR mergées | `main` à jour avant le test |
| Midi (~12:00) | **Moi** | **T65 — Test ouvert UNE SEULE FOIS** (`eval_test.py`) | `TEST_LOCK.json` (date, commit, hash, résultats) |
| Après-midi | **Sofiene** | T66 — Inspection humaine de 20 séquences par modèle (hors toxines) ; section « évaluation » du rapport | Notes qualitatives |
| Après-midi | **Youssef** | T66 — Section « modèles » du rapport ; relance bit-à-bit d'**une** config (T67) | Même perte qu'avant |
| Après-midi | **Yassine** | `PROVENANCE.md` final ; relance de contrôle depuis `data/raw/` | Même comptages, mêmes familles |
| Après-midi | **Moi** | T66 — Rapport final, commentaire de l'écart val/test ; T67 — **tag `v0.3`**, merge vers `main` | `Docs/Semaine3_rapport.md`, tag `v0.3` poussé |
| Soir | **Tous** | T68 — Répétition de la réunion (30 min) | Support de présentation prêt avant minuit |

---

## Mardi 03/11 — Réunion M. Abdallah uniquement

| Moment | Qui | Tâche | Livrable |
|--------|-----|-------|----------|
| Réunion | **Tous + M. Abdallah** | Présentation des visuels (pas du code) | Échelle de perplexités, courbes train/val, tableau comparatif avec références, `TEST_LOCK.json`, limites |
| Ordre de passage | **Yassine** → **Youssef** → **Sofiene** → **Moi** | Socle, configs et traçabilité → architectures et entraînement → échantillonnage et évaluation → synthèse, test verrouillé et limites | Chacun explique sa partie sans ouvrir le code |

**Aucune autre tâche n'est planifiée le mardi.**

---

## Règles du planning

- **Personne ne travaille sur `main` directement.** Chaque étudiant pousse sur sa branche `feat/nom-tâche` (ex. `feat/youssef-lstm`, `feat/sofiene-sampling`, `feat/yassine-config-loader`) et ouvre une PR.
- **Une PR ne merge que si :** un autre étudiant l'a approuvée + `pytest` vert (Python 3.11 **et** 3.14) + moi (chef) ai validé.
- **Le test n'est ouvert qu'une fois**, lundi vers 12:00, par le chef. Aucun réglage n'est choisi sur le test.
- **Tout ce qui est critique pour les autres est livré par Yassine avant ses coupures** (vendredi 14:00 pour le chargeur, samedi 12:00 pour T58 et le script des graines).
- **Une seule nouveauté à la fois** : chargeur testé avant modèle, modèle testé avant entraînement, entraînement validé avant génération.
- **Si une tâche bloque**, on le dit au meet suivant (vendredi matin, samedi 11:30, dimanche soir), pas mardi.

---

## Plan de repli (rappel)

| Situation | Décision |
|-----------|----------|
| CPU trop lent | Colab T4 dès samedi, checkpoints sur Drive ; les règles (test, protocole) ne changent pas |
| Graines 43/44 non terminées lundi midi | Test avec les graines terminées ; écrire « n = 1 » ou « n = 2 » dans le rapport |
| Transformer instable (`NaN`) | Baisser le taux d'apprentissage, warmup, réduire la taille ; documenter |
| PP val ne bat pas l'unigramme | À rapporter honnêtement avec l'analyse ; ne pas toucher au test |
| Livrable de Yassine en retard à 16:00 vendredi / 12:00 samedi | Sofiene ou Youssef reprennent T50 en priorité ; T58 peut glisser à dimanche/lundi |
| Reporté sans bloquer mardi | Ablations (T62), analyse de mémorisation (T63) : à noter comme suite, avec la raison |
