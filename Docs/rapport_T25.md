# T25 — Taux de quasi-doublons à 90 % d'identité

**Auteur :** Youssef Jarraya (Évaluation) — semaine 1, 18/10/2026

## Résultats

Corpus filtré (longueur 10–50 aa) : `data/processed/corpus_filtered_all.csv`.

### Redondance exacte (T22)

| Périmètre | Séquences | Uniques | Redondance |
|---|---|---|---|
| uniprot | 11 821 | 8 911 | 24,6 % |
| amp | 1 084 | 1 026 | 5,4 % |
| bactibase | 177 | 161 | 9,0 % |
| toutes sources | 13 082 | 9 021 | 31,0 % |

### Quasi-doublons à 90 % (T25)

| Mesure | Valeur |
|---|---|
| Séquences uniques | 9 021 |
| Paires avec identité ≥ 90 % | 5 675 |
| Séquences impliquées dans au moins une paire | 2 988 |
| **Taux de quasi-doublons** | **33,1 %** |

Identité = `1 − Levenshtein(a, b) / max(|a|, |b|)`. Une paire est un quasi-doublon si l'identité est ≥ 0,90.

## Vérification indépendante

Le résultat de `peplib.dedup.compute_near_duplicates_90` a été recalculé séparément :
- index de 3-mers, avec un filtre sûr : deux séquences à ≥ 90 % d'identité partagent au moins `D − 3·d_max` 3-mers distincts (`D` = nombre de 3-mers distincts de la plus longue, `d_max = ⌊0,1 × longueur max⌋`) ;
- distance de Levenshtein exacte sur les paires retenues.

Résultat identique : 5 675 paires, 2 988 séquences, taux 0,3312.

## Analyse

Un tiers du corpus est constitué de quasi-doublons. Avec un découpage aléatoire train / validation / test, des séquences quasi identiques se retrouveraient des deux côtés : le modèle serait évalué sur ce qu'il a déjà vu (fuite de données) et les performances seraient surestimées.

## Justification : clustering reporté au mois 2

En semaine 1, on mesure seulement l'ampleur du problème. Le regroupement complet (composantes connexes, union-find à 90 %) est nécessaire au mois 2 pour deux raisons :
1. Les quasi-doublons forment des chaînes (A proche de B, B proche de C, mais A pas proche de C). Compter des paires ne suffit pas, il faut des clusters.
2. Le découpage train / validation / test (90/5/5) doit se faire par clusters disjoints, avec alignement (Needleman-Wunsch / Smith-Waterman, BLOSUM62), pour garantir l'absence de fuite.
