# A6 — Borne exacte du préfiltre 3-mers pour le critère de couverture

**Analyse :** Sofiene Ouanes (via l'agent de code), lundi 26/10/2026.
**Prérequis :** mise en garde §5 de `AJUSTEMENTS_Semaine2_reconciliation_Abdallah.md`
et demi-page A2 (PR #18). Cette page re-démontre la borne pour le critère de
**couverture** — elle **ne réutilise pas** la formule `max_distance((1−τ)·max(len))`
de la PR #16, valable uniquement pour l'identité Levenshtein globale.

## 1. Notations

Soit a et b deux peptides de longueurs `|a|`, `|b|`, `m = min(|a|, |b|)`, et
supposons qu'il existe un alignement local (SW, o = 11, e = 1, BLOSUM62) :

- `W` = nombre de colonnes de l'alignement,
- `i` = colonnes identiques (résidus réels égaux, pas de gap),
- `F_a` = nombre de résidus réels de a dans l'alignement (idem `F_b`).

Le critère de liaison A5 exige :

    identité   = i / W ≥ τ            avec τ = 0,90
    couverture = max(F_a, F_b) / m ≥ τ

## 2. Lemme — 3-mers intacts

Un **3-mer intact** est 3 colonnes consécutives sans gap dont les résidus
coïncident. Un tel 3-mer correspond à une même sous-chaîne de longueur 3
présente dans a et dans b.

Supposons le critère satisfait pour cet alignement. Alors :

1. colonnes non identiques (mismatch, insertion ou délétion) :
   `W − i ≤ (1 − τ)·W` ;
2. toute fenêtre de k colonnes consécutives qui ne contient aucune colonne non
   identique délivre un k-mer intact ;
3. le nombre de fenêtres k-totales vaut `W − k + 1 ≥ W − k`, et chaque colonne
   non identique appartient à au plus k fenêtres, donc
   `intacts ≥ (W − k) − k·(1 − τ)·W = W·(1 − k(1 − τ)) − k` ;
4. couverture ≥ τ ⇒ `W ≥ max(F_a, F_b) ≥ τ·m` (l'alignement couvre au moins
   `τ·m` colonnes non-gap de la séquence la plus courte).

En combinant :

    intacts ≥ τ·m·(1 − k(1 − τ)) − k.

Avec τ = 0,90 et k = 3 :

    intacts ≥ 0,9·m·(1 − 0,3) − 3 = 0,63·m − 3.
    (le −k final approche les fenêtres manquantes en bord de fragment)

Pour m ≥ 10 (fenêtre du corpus), `intacts ≥ 0,63·10 − 3 = 3,3` : **tout
alignement satisfaisant le critère impose des occurrences intactes de 3-mers
en nombre supérieur ou égal à ⌊0,63·m − 3⌋** — l'arrondi à l'entier inférieur
est conservatif : la borne reste sous la valeur réelle et ne rogne jamais
le rappel.

## 3. Borne supérieure utilisable pour couper (coupure sonore)

Pour un 3-mer s, le nombre d'occurrences alignées par paires de s ne peut
dépasser `min(c_a(s), c_b(s))` où `c_x(s)` est le nombre d'occurrences de s
dans x. La somme

    U(a, b) = Σ_s min(c_a(s), c_b(s))

est donc une **borne supérieure** du nombre total d'occurrences alignées de
3-mers intacts (elle peut surestimer : des occurrences détachées dans
l'alignement sont comptées ; elle ne peut pas sous-estimer).

**Rejet justifié** : si `U(a, b) < K(m)` avec `K(m) = max(0, ⌊0,63·m − 3⌋)`,
alors aucun alignement ne satisfait (identité ≥ 0,90 ∧ couverture ≥ 0,90) :
la paire est rejetée **sans** Smith-Waterman, et ce rejet est sonore
(aucun faux négatif possible : U est une borne supérieure).

Toute paire avec `U(a, b) ≥ K(m)` est conservée comme candidate (le sur-ensemble
peut contenir des faux positifs, jamais de faux négatifs : rappel garanti 1,0).

## 4. Rapprochement avec le paramétrage de la PR #18 (A2)

La demi-page A2 (PR #18) démontre, pour l'identité Levenshtein **globale**,
que m k-mers partagés garantissent le rappel au seuil 0,90 (une opération
d'édition casse au plus k k-mers). Ici la quantité à garantir est différente :
un **fragment** couvrant `τ·m` résidus, à l'identité locale `τ`. La même
méthode (comptage des fenêtres cassées) est appliquée, mais avec deux
corrections qui changent la formule :

- la longueur de référence est la **couverture** (`W ≥ τ·m`) et non la longueur
  totale : c'est le terme `τ·(1 − k(1−τ))` (proche de τ²) que la formule
  `max_distance` de la PR #16 n'avait pas ;
- le comptage des occurrences est **multiset** (`c_x(s)`) et non par types :
  un k-mer répété peut être intact plusieurs fois, et couper sur les types
  seuls permet des faux négatifs (un fragment répétitif partage un seul type
  mais 20 occurrences intactes).

L'arbitrage « types vs occurrences » annoncé dans DECISIONS.md est tranché ici
au profit des **occurrences (multiset)**, faute de quoi la borne ne trouve pas
d'équivalent — les candidats de la PR #18 sont conservés en mesure
complémentaire.

## 5. Portée et limites

- La borne est démontrée pour le critère (identité, couverture) ≥ (τ ; τ) :
  elle est ici appliquée en τ = 0,90 et en τ = 0,70 (remplacer τ dans §2).
- Elle est **démontrée** : aucune méthode « epsilon » ou coefficient approché
  n'est utilisée ; seul l'arrondi ⌊·⌋ conservatif est appliqué.
- Le critère de liaison sur couverture est le seul qui détecte le motif §5
  (fragment exact, sous-chaîne d'un peptide d'entraînement, identité
  globale 0,4) — cf. test de veille dans `tests/test_identite.py`.
