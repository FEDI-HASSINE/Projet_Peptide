# Planning Semaine 2 — 5 jours (vendredi 23 – mardi 27/10/2026)

**Groupe :** Moi (Chef/Coordinateur), Yassine Achour (Données), Sofiene Ouanes (Modèles), Youssef Jarraya* (Évaluation)
**Objectif :** Corpus enrichi **v0.2** (UniProt + APD3 + base AMP Challenge fusionnés sans doublons), **clusters à 90 %**, **splits étanches**, **tokenisation**, **baselines**, tout testé, tracé et mergé. Prêts pour la réunion M. Abdallah le mardi.

---

## Vendredi 23/10 — Clôture S1 + démarrage APD3 (travail individuel)

| Moment | Qui | Tâche | Livrable |
|--------|-----|-------|----------|
| Matin | **Tous** | T27 — Demi-page d'analyse (1 paragraphe par étudiant, 6 questions) | Demi-page rédigée et relue |
| Matin | **Moi** | T29 — Tag `v0.1` + relance complète depuis `data/raw/` (sans re-télécharger) | Tag poussé, mêmes comptages chez les 4, `pytest` vert |
| Après-midi | **Yassine** | **T32 — Licences : GO / NO-GO** pour les deux sources : APD3 (termes d'usage/redistribution) **et base AMP Challenge** (`antibacterial.fasta`, dépôt BSD-3 ; 13 bases amont à consigner) + téléchargement des 2 bruts + entrées `PROVENANCE.md` | Décisions écrites + fichiers bruts dans `data/raw/` avec SHA-256 (hash challenge : `cbbeac64…21e9`) |
| Après-midi | **Sofiene** | T35 (début) — Union-find : clusters à 90 % sur le corpus v0.1 | `peplib/cluster.py` en cours |
| Après-midi | **Youssef** | T33 — BLOSUM62 saisie depuis la littérature + tests sur scores connus | Table dans `constants.py` + `pytest` vert |
| Après-midi | **Moi** | Création des branches `feat/` de la semaine 2 + règles de PR rappelées | Branches prêtes |

**Ce soir :** Les deux licences sont tranchées. Si APD3 = NO-GO → filet : la base AMP Challenge suffit comme corpus AMP (elle contient DRAMP, DBAASP, APD…). Aucun merge sans les deux suites pytest vertes.

> **Règle :** aucune fusion de corpus tant que le clustering à 90 % n'est pas testé — c'est lui qui empêche de compter deux fois le même peptide.

---

## Samedi 24/10 — Meet de récap + connecteur APD3 + clusters

| Moment | Qui | Tâche | Livrable |
|--------|-----|-------|----------|
| Matin (Meet) | **Tous** | Récap de la clôture S1, répartition des tâches, questions | Tout le monde sait quoi faire |
| Matin | **Yassine** | T34 — **Connecteurs APD3 + AMP Challenge** → schéma commun (`sequence, source, id` + métadonnées). APD3 : **résidus modifiés** (amidation C-terminale, D-acides aminés, B/X/U/O) — décision écrite ; Challenge : parse des en-têtes (`len/charge/disulfide/activity/dbs`, 39 448 séquences, 0 résidu modifié) | `peplib/download_apd3.py` + `peplib/download_challenge.py` + décisions documentées + tests |
| Matin | **Sofiene** | T35 (suite) — Clusters terminés : tests (corpus connu, cas limites), distribution des tailles de familles | `clusters_90.csv` + figure + nb de familles |
| Après-midi | **Youssef** | T36 — Needleman–Wunsch global avec pénalités affines (Gotoh) | `peplib/align.py` + alignements de référence recalculés à la main |
| Après-midi | **Moi** | Reviews + merge des PR de la journée | `main` toujours vert |

**Ce soir :** Les clusters sont stabilisés et le connecteur APD3 fonctionne sur les données brutes.

---

## Dimanche 25/10 — Charge maximale (pour alléger lundi)

| Moment | Qui | Tâche | Livrable |
|--------|-----|-------|----------|
| Matin | **Yassine** | T37 — Doublons exacts inter-sources (colonne `sources`, atomes `{UniProt, APD3, AMPChallenge}`, jusqu'à 7 combinaisons) + **corpus fusionné v0.2** | `corpus_v0.2.csv`, invariant : aucune séquence perdue |
| Matin | **Yassine** | T38 — Rapport de recouvrement **trilatéral** UniProt ↔ APD3 ↔ AMPChallenge | Chiffres (3 paires + triple : exacts / quasi-doublons / uniques) + **figure** |
| Matin | **Sofiene** | T39 — Split **stratifié par source d'ingestion (3)** 90/5/5 **par familles** + **preuve d'étanchéité** (identité max train/test < 90 %, aucune famille à cheval) | 3 fichiers de split + rapport chiffré + test anti-fuite |
| Matin | **Youssef** | T36 (suite) — Smith–Waterman local | `peplib/align.py` complet + tests |
| Après-midi | **Youssef** | T40 — Tokenisation **V = 23** (20 AA + PAD/BOS/EOS) + round-trip encode/decode sur 1 000 séquences | Tokeniseur + test round-trip 100 % |
| Après-midi | **Moi** | T41 — `corpus_manifest.json` (total, **comptes par source d'ingestion (3)**, nb de familles, hashes CSV + 2 bruts AMP, version `v0.2`) + test anti-drift | Manifeste + test |
| Soir | **Moi** | Relance des analyses sur v0.2 (JS, descripteurs, redondance) + figures mises à jour | Chiffres v0.2 cohérents |

---

## Lundi 26/10 — Allégé (finalisation uniquement)

| Moment | Qui | Tâche | Livrable |
|--------|-----|-------|----------|
| Matin | **Yassine** | T37 (suite) — Vérifier l'invariant : `sources` traçable, chaque source extractible (`uniprot.csv`, `amp_apd3.csv`, `challenge.csv`, `corpus_v0.2.csv`) + `PROVENANCE.md` mis à jour (deux licences + 13 bases amont, hashes, version v0.2) | PROVENANCE v0.2 + pipeline régénère les 4 vues |
| Matin | **Sofiene** | T42 — Dédoublonnage quasi-doublons inter-sources via les clusters (1 représentant par famille) | Liste de représentants + mapping famille → sources |
| Après-midi | **Youssef** | T43 — **Échelle de baselines** : aléatoire → unigramme (21 tokens) → **bigramme** (21×21), perplexité en nats/lettre, seeds écrites | 3 baselines + tableau comparatif |
| Après-midi | **Sofiene** | T44 / T45 (si temps) — Validation croisée des descripteurs (APD3 publie charge/pI/masse ; **les en-têtes du challenge donnent la charge de 39 448 séquences**) + **témoins négatifs** (500 AMP mélangés au hasard) | 2 figures de contrôle qualité |
| Après-midi | **Tous** | T47 — Séminaire interne 30 min : présentation de la v0.2 (visuels, pas de code) | Présentation prête |
| Soir | **Moi** | T46 — Relance complète depuis `data/raw/`, tag `v0.2`, merge final vers `main` | Tag `v0.2` poussé, dépôt `main` à jour, `pytest` vert |

---

## Mardi 27/10 — Réunion M. Abdallah uniquement

| Moment | Qui | Tâche | Livrable |
|--------|-----|-------|----------|
| Réunion | **Tous + M. Abdallah** | Présentation des visuels de la v0.2 | Recouvrement trilatéral, distribution des familles, preuve d'étanchéité des splits, échelle de baselines, figures relancées, bloc PROVENANCE challenge |

**Aucune autre tâche n'est planifiée le mardi.**

---

## Livrables et preuves attendus (definition of done)

| Livrable | Preuve exigée |
|----------|---------------|
| Corpus v0.2 fusionné | Colonne `sources` (3 atomes : UniProt / APD3 / AMPChallenge), doublons exacts éliminés, aucune séquence perdue |
| Clusters à 90 % | `clusters_90.csv` + nombre de familles + distribution des tailles |
| Splits | 90/5/5 stratifiés par source d'ingestion (3) + **identité max train/test < 90 %** (chiffre rapporté) |
| Tokeniseur | Round-trip encode/decode 100 % sur 1 000 séquences |
| Baselines | Perplexité (nats/lettre) aléatoire / unigramme / bigramme, seeds écrites |
| Manifeste | `corpus_manifest.json` + test anti-drift |
| PROVENANCE | Licences APD3 + AMPChallenge consignées (BSD-3 + 13 bases amont), hashes, version du corpus |
| Tag | `v0.2` poussé sur `main` |

## Ce qu'on ne fait PAS cette semaine

- Aucun modèle (LSTM / Transformer) — c'est la semaine 3.
- Aucune base en plus des trois sources validées (BIOPEP, DADP, etc. en dehors) — les deux AMP validées (APD3 + challenge, fourni par l'encadrant) suffisent ; on n'ouvre pas autre chose.
- Aucun outil externe (BLAST, MMseqs2, cd-hit, HMMER) — alignement et clustering écrits par nous.
- Aucune affirmation biologique : les descripteurs sont des proxys, pas des mesures d'activité.

## Règles du planning

- **Personne ne travaille sur `main` directement.** Chaque ami pousse sur sa branche `feat/nom-tâche` et ouvre une PR.
- **Une PR ne merge que si :** un autre étudiant l'a approuvée + `pytest` vert (Python 3.11 **et** 3.14) + moi (chef) ai validé.
- **Les licences sont des GO / NO-GO** (APD3 et AMPChallenge) : sans décision écrite, pas de téléchargement.
- **Une seule nouveauté à la fois** : clusters testés avant fusion, fusion avant baselines.
- **UTF-8 propre, sans BOM, sans mojibake** ; seed 42 partout ; aucun chemin absolu dans le code.
- **Si une tâche bloque**, on en parle au meet du samedi matin, pas mardi.
