# Prompt — Re-revue de la PR #26, point vert et fusion

> À donner tel quel à l'agent chargé de la re-revue.

---

## Rôle et mission

Tu es relecteur de la **re-revue** de la pull request **#26** du dépôt
`FEDI-HASSINE/Projet_Peptide` :
https://github.com/FEDI-HASSINE/Projet_Peptide/pull/26

Contexte : une première revue humaine (CHANGES_REQUESTED, revue #5941121685)
a bloqué la fusion sur **2 points** :

- **B1** : `peplib/metrics_proximity.py:86` — l'affectation
  `details = meilleures_paires(...)` était indentée à 8 espaces, donc exécutée
  **dans** la boucle `for j, seq_train` : le balayage Smith–Waterman complet
  était relancé une fois par séquence du train (mesuré : 0,02 s → 7,9 s sur
  5×500 ; des heures aux volumes réels de 38 336).
- **B2** : `peplib/metrics_generation.py` — `is_valid_peptide` ne vérifiait que
  l'alphabet, pas la fenêtre **10–50** : un tiers des séquences unigramme/
  bigramme (< 10 AA, jusqu'à longueur 1) étaient comptées valides
  (`validity_rate = 1.0`).

L'auteur a poussé le commit **`b2deba9` — « fix: resolve PR26 blockers B1 and B2 »**
(PR `mergeable_state = clean`, 12 fichiers au total dont
`results/week3/baselines_metrics.csv`). Ta mission : **vérifier que les deux
bloquants sont réellement corrigés, qu'aucune régression n'a été introduite,
puis — si et seulement si tout est vert — approuver la PR et la fusionner**
(point vert de merge). Verdict contraire (CHANGES_REQUESTED) si un gate échoue.

Le compte `gh` authentifié est **FEDI-HASSINE**, le/la coordinateur·rice :
la règle « seule la personne coordinatrice fusionne » est donc respectée.
**Ne fusionne jamais** si un seul contrôle échoue.

## Installation et état à vérifier

```powershell
git fetch origin
git worktree add ../revue-pr26-b2 origin/pr-26-head   # ou checkout b2deba9
cd ../revue-pr26-b2
python -m venv .venv ; .venv\Scripts\pip.exe install -r requirements.txt
# ou, si l'environnement système a torch/pandas/scipy, l'utiliser directement
python -m pytest -q
```

- Confirme que HEAD = **`b2deba9`** et que la base est `origin/main`
  (`git merge-base HEAD origin/main`, `git log --oneline origin/main..HEAD`) :
  exactement 2 commits, `811ebeb` + `b2deba9`, aucun commit parasite.
- Relève le nombre total de tests (attendu : 654 + les nouveaux ≈ 666 ;
  **vérifie, ne crois pas**), le temps total, et l'absence de warning bloquant.
- Vérifie `mergeable_state` réel via l'API (`clean` attendu, mais
  re-confirme après ton fetch).

## Contrôle n° 1 — B1 réellement corrigé (preuve exigée)

1. **AST** : l'affectation `details = meilleures_paires(...)` doit être au
   niveau du corps de fonction (col_offset 4) et **hors** de toute boucle
   `for`. Rejoue le contrôle :
   ```python
   import ast
   src = open("peplib/metrics_proximity.py", encoding="utf-8").read()
   tree = ast.parse(src)
   # localiser la FunctionDef proximity_details, puis tout ast.For contenant
   # un ast.Assign de cible `details` -> doit être vide
   ```
2. **Test de non-régression** : vérifie qu'un test compte désormais les appels
   à `meilleures_paires` (monkeypatch avec compteur) et échoue si > 1. S'il
   n'existe pas, c'est un **BLOQUANT** de re-revue (la correction sans test se
   réintroduit au prochain refactor).
3. **Comportement** : sur 5 générées × 500 train (ou un petit corpus), la
   durée de `proximity_details` doit être de l'ordre du dixième de seconde,
   pas de plusieurs secondes. Compare au chiffre de la revue (0,02 s vs 7,9 s).

## Contrôle n° 2 — B2 réellement corrigé (preuve exigée)

1. `is_valid_peptide_10_50` et `validity_rate_10_50` existent et sont
   exportés ; `generation_metrics` renvoie `n_valid_10_50` et
   `validity_rate_10_50` en plus des champs historiques.
2. `generate_ngram` **interdit EOS tant que `len(residues) < 10`** et
   `max_residues` est borné à **[10, 50]**. Vérifie le code, pas la docstring.
3. **Mesure** : sur les vraies sorties T61 (splits v02, seed 42),
   `length_min` de unigramme et bigramme doit être **exactement 10**, et
   `validity_rate_10_50` doit être cohérent avec le CSV commité.
   Vérifie aussi qu'aucune séquence générée n'a < 10 AA :
   `sum(len(s) < 10 for s in baselines["unigram"]) == 0`.
4. Les tests couvrent les cas limites : longueur 9 / 10 / 50 / 51,
   EOS interdit aux 9 premiers pas, EOS autorisé au 10e.

## Contrôle n° 3 — Le CSV commité est reproductible (nouveau livrable)

`results/week3/baselines_metrics.csv` (6 lignes = en-tête + 5 baselines) est
maintenant commité. Rejoue l'évaluation et compare **valeur par valeur** :

```powershell
python -m peplib.evaluate_baselines_s3
git diff --stat -- results/week3/baselines_metrics.csv   # doit être vide
```

- Si les flottants diffèrent au-delà du bruit de plateforme, exige une note de
  provenance (versions pinnées) ou un arrondi documenté.
- Vérifie que les 5 baselines sont présentes (random, unigram, bigram,
  validation, shuffled_validation) et que la colonne `validity_rate_10_50`
  est bien remplie.
- **Attention** : relancer l'évaluation prend ~4 min ; le script lit
  `data/splits/train_v02.csv` et `val_v02.csv` depuis la racine du dépôt.

## Contrôle n° 4 — Non-régression du protocole anti-fuite (inchangé, à re-prouver)

1. `evaluate_baselines_s3.py` lit toujours `data/splits/{train,val}_v02.csv`
   (38 336 / 2 095) et **jamais** `test_v02.csv`.
2. Grep sur les 12 fichiers : aucune occurrence de `make_splits`,
   `train_test_split`, `data/processed/splits`, `sample(frac=...)`.
3. Splits intacts : `git diff origin/main...HEAD --stat` ne doit montrer que
   les 12 fichiers de la PR ; aucun fichier interdit
   (`constants.py`, `align.py`, `data/splits/*`, `results/week2/*`,
   `PROVENANCE.md`, `test_artefacts_v02.py`) modifié.
4. Fins de ligne des blobs `i/lf` (`git ls-files --eol`), pas de BOM,
   `requirements.txt` inchangé, aucune dépendance ajoutée.

## Contrôle n° 5 — Périmètre et qualité

- Le commit `b2deba9` ne contient que des corrections ciblées (+290/−63
  annoncés) : pas de refactor massif non demandé, pas de fichier étranger.
- Les gates déjà validés restent valides : PAD/BOS à −∞ dans `sampling.py`,
  JS en **nats** bornée par ln 2, PP en nats/token avec PAD exclu, 23 tokens
  de vocabulaire, 21 cibles n-grammes.
- Points non bloquants de la première revue : note lesquels restent ouverts
  (JS aveugle à l'ordre via `shuffled_validation`, `alpha=1.0` en dur,
  multi-graines absent, T60 sans script d'exécution). **Ils ne bloquent pas la
  fusion** sauf régression sur un point déjà OK.

## Décision et actions GitHub

### Si tout est vert

1. Publie une revue **APPROVE** en français, avec le tableau des contrôles
   re-vérifiés et les chiffres mesurés (nombre de tests, `length_min`,
   durée T60, diff CSV vide) :
   ```powershell
   gh pr review 26 --repo FEDI-HASSINE/Projet_Peptide --approve `
     --body-file <chemin>\revue_pr26_round2.md
   ```
2. Fusionne avec un **commit de fusion** (convention du dépôt — les PR
   précédentes sont fusionnées ainsi) :
   ```powershell
   gh pr merge 26 --repo FEDI-HASSINE/Projet_Peptide --merge
   ```
   Ne supprime pas la branche, ne force-push pas, ne squash pas.
3. Confirme après fusion : `gh pr view 26 --json state,mergedAt` =
   `MERGED`, et `git log origin/main -1` contient `b2deba9` (ou son commit de
   fusion).
4. Rapporte le résultat à l'utilisateur : verdict, commandes exécutées,
   état final de `main`.

### Si un gate échoue

Publie une revue **CHANGES_REQUESTED** avec la liste exacte des points
restants (fichier:ligne + preuve), **ne fusionne pas**, et termine.

## Format de sortie

- Langue : **français**, ton factuel, chiffres exacts avec sources
  (fichier:ligne, sortie de commande).
- Verdict en tête : GO (fusion) / NO-GO (liste des bloquants restants).
- Tableau des contrôles B1 / B2 / CSV / anti-fuite / périmètre avec
  OK ou BLOQUANT.
- Ne publie rien d'autre que la revue et la fusion demandées.
