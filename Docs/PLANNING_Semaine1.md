# Planning Semaine 1 — 5 jours (vendredi 16 – mardi 20/10/2026)

**Groupe :** Moi (Chef/Coordinateur), Yassine Achour (Données), Sofiene Ouanes (Modèles), Youssef Jarraya* (Évaluation)  
**Objectif :** Corpus propre, `peplib` testé, tag `v0.1`, prêts pour la réunion M. Abdallah le mardi.

---

## Vendredi 16/10 — Fondation (travail individuel, pas de meet)

| Moment | Qui | Tâche | Livrable |
|--------|-----|-------|----------|
| Matin | **Moi** | T1 — Dépôt Git, arborescence, `.gitignore`, `requirements.txt`, `README.md` | Repo local créé, dossiers `data/`, `peplib/`, `tests/`, `notebooks/` prêts |
| Matin | **Youssef** | T2 — `pip install -r requirements.txt`, venv activé | Environnement fonctionnel |
| Matin | **Youssef** | T3 — `set_seed(42)` écrit et testé | Fonction seed qui redonne les mêmes tirages |
| Matin | **Sofiene** | T4 — `constants.py` saisi (masses, pKa, Kyte-Doolittle, Eisenberg) | `peplib/constants.py` avec références citées |
| Matin | **Youssef** | T5 — Premier test pytest (constantes) | `pytest` vert sur les constantes |
| Après-midi | **Tous** | T6 — Lecture de fond 2h (acides aminés, liaison peptidique, AMP) | Chacun sait expliquer pourquoi les AMP sont cationiques et amphipathiques (leli ma3andouch wa9t tache twali samedi sbe7 ) |

**Ce soir :** Repo poussé sur GitHub, 3 amis ajoutés collaborateurs, branche `main` protégée.

---

## Samedi 17/10 — Meet de récap + démarrage acquisition

| Moment | Qui | Tâche | Livrable |
|--------|-----|-------|----------|
| Matin (2h) | **Tous** | Recherche individuelle (complément T6) | eli 5dem vendredi(T6) c bon ken y7eb yzid pas probleme   |
| Matin (Meet) | **Tous** | Récap des 2h de recherche, clarification des tâches, questions | Tout le monde sait quoi faire |
| Après-midi | **Yassine** | T7 — Requête UniProt (fichier brut `data/raw/uniprot_*.tsv`) | Fichier brut téléchargé |
| Après-midi | **Yassine** | T8 — Base secondaire (APD3 / DBAASP / DRAMP) | Fichier brut dans `data/raw/` |
| Après-midi | **Yassine** | T9 — Source non-pharma (BIOPEP-UWM / BACTIBASE) | Fichier brut dans `data/raw/` |
| Après-midi | **Yassine** | T10 — `PROVENANCE.md` commencé | Une entrée par source avec SHA-256 |
| Après-midi | **Sofiene** | T18 — Descripteurs écrits à la main (charge, pI, GRAVY, moment hydrophobe, Boman, masse) | `peplib/descriptors.py` |
| Après-midi | **Moi** | Création branches `feat/` + règles de protection `main` | Branches prêtes pour les 3 amis |

> **Règle :** T7 (UniProt) commence **uniquement après le meet de samedi matin**, pas le vendredi soir.

---

## Dimanche 18/10 — Charge maximale (pour alléger lundi)

| Moment | Qui | Tâche | Livrable |
|--------|-----|-------|----------|
| Matin | **Yassine** | T11 — CSV par source (`sequence`, `source`, `id`) | `data/processed/` contient les CSV |
| Matin | **Yassine** | T12, T13, T14 — `dropna()` avant `astype(str)`, filtre regex 10–50, table des rejets | `conservés + rejetés = total brut` |
| Matin | **Sofiene** | T20 — PWM alignée N/C sans troncature | `peplib/pwm.py` |
| Matin | **Sofiene** | T21 — Figures (distribution, heatmap PWM, courbe info) | Graphiques `matplotlib` dans `results/` |
| Matin | **Youssef** | T15 — Test du filtrage (NaN, vide, X, U, longueurs 9 et 51) | `pytest` vert sur les rejets attendus |
| Après-midi | **Yassine** | T16 — Histogrammes + composition `f(a)` | `sum(f(a)) == 1` par source |
| Après-midi | **Yassine** | T17 — Matrice JS divergence entre sources | `JS(p,p)=0`, bornée par ln 2 |
| Après-midi | **Yassine** | T22 — Doublons exacts + taux de redondance | Taux par source et entre sources |
| Après-midi | **Sofiene** | T23 — Quasi-doublons (index 3-mers + Levenshtein) | `Levenshtein("kitten","sitting")==3` |
| Après-midi | **Yassine** | T25 — Taux de quasi-doublons à 90 % | Rapport écrit + justification clustering mois 2 |
| Après-midi | **Youssef** | T19 — Test descripteurs sur peptides à la main | `GRAVY("A")=1.8`, masse `"G"`≈75,0320, etc. |
| Après-midi | **Youssef** | T24 — Test redondance + entropie de composition | `H = -Σ f(a) log f(a)` |
| Soir | **Tous** | T26 — Notebook `NB01` (appels `peplib` uniquement) | Notebook qui tourne sans logique métier dans les cellules |

---

## Lundi 19/10 — Allégé (finalisation uniquement)

| Moment | Qui | Tâche | Livrable |
|--------|-----|-------|----------|
| Matin | **Tous** | T27 — Demi-page d'analyse (1 paragraphe par étudiant) | 6 questions d'analyse traitées |
| Après-midi | **Tous** | T28 — Séminaire interne 30 min (répétition pour M. Abdallah) | Présentation prête |
| Après-midi | **Tous** | T29 — Tag `v0.1`, relance complète depuis `data/raw/` | Mêmes comptages chez les 4, `pytest` vert |
| Après-midi | **Moi** | Merge final vers `main` uniquement | Dépôt `main` à jour, tag `v0.1` poussé |

---

## Mardi 20/10 — Réunion M. Abdallah uniquement

| Moment | Qui | Tâche | Livrable |
|--------|-----|-------|----------|
| Réunion | **Tous + M. Abdallah** | Présentation des visuels (pas du code) | Histogrammes, heatmap PWM, table des rejets, matrice JS, `PROVENANCE.md`, demi-page d'analyse |

**Aucune autre tâche n'est planifiée le mardi.**

---

## Règles du planning

- **Personne ne travaille sur `main` directement.** Chaque ami pousse sur sa branche `feat/nom-tâche` et ouvre une PR.
- **Une PR ne merge que si :** un autre étudiant l'a approuvée + `pytest` vert + moi (chef) ai validé.
- **Si une tâche bloque**, on en parle samedi matin au meet, pas mardi.