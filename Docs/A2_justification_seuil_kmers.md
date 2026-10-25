# A2 — Justification du seuil de partage des 3-mers

**Responsable :** Sofiene  
**Date :** 25 octobre 2026  
**Objet :** justification du préfiltrage des paires de peptides à 90 % d'identité.

## Principe

Le préfiltrage repose sur un index inversé associant chaque 3-mer aux séquences qui le contiennent. Une paire est conservée si elle partage au moins `m` 3-mers distincts.

Pour une séquence de longueur \(L\), le nombre d'occurrences de k-mers est \(L-k+1\). Une opération d'édition affecte au plus \(k\) occurrences. Si la distance de Levenshtein entre deux séquences est au plus \(d\), au moins \(\max(0,L-k+1-kd)\) occurrences restent préservées, en prenant \(L\) comme longueur de la séquence de référence.

Pour des peptides de longueurs 10 à 50 AA, au seuil de 90 % d'identité Levenshtein normalisée par la plus grande longueur, cette borne est strictement positive avec \(k=3\). Cela garantit l'existence d'au moins un 3-mer commun.

**Choix : `min_shared=1`.**

L'index actuel compte les 3-mers distincts, et non leurs occurrences. Il serait donc incorrect de déduire directement un seuil supérieur à 1 du nombre d'occurrences préservées.

## Validation

Le préfiltre a été comparé à une recherche exhaustive sur un échantillon reproductible de 200 séquences uniques (seed 42). Toutes les paires dont l'identité Levenshtein est supérieure ou égale à 90 % ont été conservées par le préfiltre.

Le test `test_kmer_rappel_echantillon` est réussi.

Un cas de répétitions (`AAAAAAAAAA` / `AAAAAAAACA`) a également montré que `min_shared=3` manque une paire admissible à 90 %. Le seuil est désormais paramétrable dans `compute_near_duplicates_90()`.

## Limite

Cette démonstration concerne l'identité globale fondée sur Levenshtein. Elle ne garantit pas automatiquement le rappel pour l'identité de couverture issue de Smith–Waterman demandée pour le corpus v0.2. Une justification supplémentaire est nécessaire avant d'appliquer ce préfiltre au clustering de couverture.

**Conclusion :** `min_shared=1` est un seuil conservateur justifié pour le critère Levenshtein dans le domaine de longueurs étudié. Le rappel empirique observé sur l'échantillon est de 100 %.
