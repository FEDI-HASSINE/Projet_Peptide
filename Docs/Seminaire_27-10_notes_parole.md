# Séminaire mardi 27/10 — ce qu'on montre au Prof. Abdallah (Semaine 2)

> Ordre : moi (vue générale) → Yassine (données) → Sofiene (familles/étanchéité) → Youssef (évaluation).
> Figures prêtes dans `wtsof\results\week2\` : `cluster_sizes_v02.png`,
> `courbe_fuite.png`, `ecart_baselines.png`, `temps_align.csv` (tableau).
> Chiffres exacts = sorties de commandes commitées (PR #21 + PR #20 fusionnée).

---

## 0. Moi — vue générale (2 min)

**En une phrase :** la semaine 2 a transformé les téléchargements bruts en un
**corpus verrouillé → des familles → des splits étanches → des baselines
chiffrées**, le tout sans outil externe (BLAST interdit) et avec preuve à chaque
étape.

*Image simple à dire :*
> « On a d'abord **lavé et rangé la cuisine** avant de cuisiner : un corpus
> propre (v0.2), rangé en **familles**, séparé en **trois casiers** qui ne
> communiquent pas (train/val/test), et un **premier goûteur très simple**
> (n-gram) pour fixer la marque à battre. Aucune étape n'est “crue” : chacune
> a une empreinte SHA-256 et un test qui pleure si on la touche. »

Trois résultats à poser dès le début (les "3 piliers") :
1. **Corpus v0.2** : 42 575 séquences uniques multi-sources, verrouillé par
   empreintes SHA-256 (tableau + manifeste).
2. **Familles v0.2** : 20 673 familles pour une seule règle —
   *figure `cluster_sizes_v02.png`*.
3. **Étanchéité mesurée, pas déclarée** : F(0,90) = 0 —
   *figure `courbe_fuite.png`*.

Le fil rouge : « **chaque critique qu'on peut nous faire, on l'a anticipée
par un test verrouillé** » (fuite, familles à cheval, grande famille dans le
mauvais split, réutilisation de code non prouvé).

---

## 1. Yassine — données et biosécurité (3–4 min)

- Corpus v0.2 multi-sources : **AMPChallenge 33 974 / UniProt 9 021 /
  APD3 2 796** ; vues par source + vues croisées (T39).
- Toxines (clause W4, train) : **unknown 29 517 / not_flagged 7 017 /
  flagged 1 785** ; règle : `unknown` = traité comme `flagged` pour tout
  usage génératif (mode d'échec sûr) ; **0 séquence** du corpus v0.2 absente
  du signalement. Mention honnête : limite G4 — 110 UniProt `unknown` sans
  annotation toxine, pas de rattrapage possible sans outil.

*Image simple à dire :*
> « C'est comme un entrepôt de bouteilles : celles marquées “poison”
> (flagged) sont mises de côté, et celles **sans étiquette** (unknown),
> on les traite **comme du poison par précaution** — on ne les donne jamais
> à un modèle qui invente de nouvelles molécules. »

- Artefacts immuables : `data/raw/` non re-touché, SHA-256 consignés dans
  `PROVENANCE.md`.

## 2. Sofiene — familles, étanchéité (5 min, les deux figures)

*(1) D'abord le coût, pour justifier la méthode — `temps_align.csv` :*
- SW sur **906 294 025** paires tous-contre-tous = **63 h (m=10) à 1 400 h
  (m=50)** ⇒ impossible ⇒ préfiltre mathématiquement borné.
- Vérification exhaustive d'abord (A9) : NW/SW testés contre l'énumération
  de **1 000 paires** ≈ tous les chemins ≤ longueur 4 — divergence = 0 ;
  `align.py` non modifié d'un octet.

*Image simple pour le coût :*
> « Aligner deux peptides, c'est les caler lettre par lettre comme deux
> textes. 906 millions de paires × un millième de seconde ≈ **deux mois
> de calcul**. Donc on ne regarde la loupe que sur les candidats sérieux. »

*(2) Le préfiltre + le clustering — `cluster_sizes_v02.png` :*
- Borne démontrée (occurrences intactes ≥ ⌊0,63·m − 3⌋) appliquée aux
  3-mers partagés : 116 206 731 couples qui partagent au moins un 3-mer →
  **1 061 432 candidats** seulement → SW complet → **20 673 familles**,
  plus grande **4 541 (10,67 %)**, 17 592 singletons.
- Critère : identité ≥ 0,90 **ET** couverture ≥ 0,90 sur le meilleur
  alignement local (l'identité seule rate les fuites-fragment).

*Image simple pour la borne (le point le plus “math” de la séance) :*
> « Pour que deux séquences ressemblent à 90 %, elles doivent garder les
> mêmes morceaux de 3 lettres : on a PROUVÉ que si la ressemblance est
> possible, la paire passe forcément dans nos candidats. On ne saute
> personne — c'est ça, une coupure **saine**. »

*Image simple pour les familles :*
> « Une famille = le même peptide retrouvé plusieurs fois, recopié ou
> modifié d'une lettre. Si on laissait le même peptide dans train ET dans
> test, ce serait une question posée à l'examen déjà vue en révision. »

*(3) Splits + fuite — `courbe_fuite.png` :*
- 90/5/5 **par familles entières** : train 38 336 / val 2 095 / test 2 144
  (±212,9 respecté), familles > 1 % forcées en train, stratification par
  source mesurée (val : 1 743/459/141 vs attendu 1 699/451/140), zéro
  famille à cheval.
- **F(0,90) = 0,000000** : aucune des 2 095 séquences de val n'a un lien
  ≥ (90 % ; 90 %) contre train ; le **contrôle aléatoire fuit** : 55,8 %
  (0,90) et 70 % (0,70).
- F(0,70) exact du split par familles achevé ce soir (6 710 520
  alignements ; borne actuelle 0,333).

*Image simple pour la figure :*
> « Deux barres : notre split étanche = **0 sosie** ; un mixage aléatoire
> = **56 % des séquences de val ont un sosie dans le train**. Sans la
> courbe, on dit “on a fait attention” ; avec la courbe, on le PROUVE. »

## 3. Youssef — évaluation, tokenizer, baselines (4 min, figure + tableau)

- A9-mesure : invariants sur **500 paires réelles** du corpus — **0
  violation** (`align_invariants_a9.csv`) ; lui énumère l'optimum, moi je
  mesure les invariants : deux contrôles indépendants.
- Tokenizer **V = 23** (20 AA + PAD/BOS/EOS), round-trip 100 %.

*Image simple pour le tokenizer :*
> « Chaque peptide devient une suite de nombres fixes : un “début de
> phrase” (BOS), la liste des lettres, le “point final” (EOS), et du
> matelas (PAD) pour rendre toutes les phrases de même longueur. »

- n-gram ordres 0–3, **21 cibles (20 AA + EOS)** — ancre exacte : uniforme
  = **ln 21 ≈ 3,0445** nats/token ; α ∈ {0,01 ; 0,1 ; 1} choisi sur val.

*Image simple pour un n-gram :*
> « C'est la saisie prédictive des SMS : l'ordre 0 devine chaque lettre au
> hasard parmi 21 ; l'ordre 1 regarde les fréquences des lettres ; l'ordre 3
> lit les DEUX lettres précédentes avant de deviner. Plus l'ordre monte,
> plus le modèle est “intelligent” — mais plus il a besoin d'exemples. »

- *Figure `ecart_baselines.png` + `baseline_ngram_A7.csv` (vrais splits)* :

  | split | ordre | CE nats/token | PP |
  |---|---|---|---|
  | familles | 0 | 3,0445 | 21,00 |
  | familles | 1 | 2,9698 | 19,49 |
  | **familles** | **3** | **2,9717** | **19,53** |
  | aléatoire | 3 | 2,7990 | 16,43 |

  **Le point à vendre :** sur le split par familles le trigramme
  **régresse** (2,97 > 2,96 du bigramme) = preuve indirecte que le split
  empêche la fuite : comme l'examen (val) vient de familles différentes,
  mémoriser deux lettres de contexte **n'aide pas** . Sur le split
  aléatoire, le trigramme “progresse” seulement parce qu'il a déjà vu la
  famille dans train. GATE mois-2 fixée : battre **PP < 19,49**.

*Image simple pour la régression du trigramme :*
> « Un trigramme qui mémorise “après W-V vient D” devient bon SI train et
> val parlent du même style. Sur notre split, val vient d'une famille que
> le trigramme n'a jamais vue : sa mémoire le dessert — courbe verte
> (familles) qui remonte, courbe orange (aléatoire) qui descend. »

- Datasets PyTorch (A8) : padding à droite + labels −100 sur PAD ;
  reproductibilité **bit à bit** (même seed = mêmes batchs, testé).

## 4. Moi — clôture (2 min) : ce qu'on demande au prof

1. **Arbitrer la définition d'identité** : tuple (identité, couverture) sur
   le meilleur SW (nôtre, imposée ; ancre fragment 20/50 → 0,4 / 1,0) vs
   float min/max-len sans gaps (PR #20) — les deux coexistent sans crash,
   un seul nom doit rester (décision mardi).
2. **Valider la gate mois-2** : PP < 19,49 (unigramme 21 cibles, split familles).
3. Informer que T44/T45 sont reportés **jeudi** (raison écrite : budget SW
   consommé par le clustering réel), et que la variante couverture 0,70 sur
   le corpus entier est reportée pour la même raison (le 0,70 est mesuré
   exactement par la courbe de fuite).

---

## Annexe — Questions probables du prof et réponses courtes

- **« Pourquoi pas BLAST/MMseqs2 ? »** → Interdits par les règles du projet ;
  on a tout codé (Gotoh par-dessus BLOSUM62) et **prouvé** par énumération
  exhaustive. Le prix, on le paie en heures de calcul — assumé (63→1 400 h
  sans préfiltre ; ~1 000 alignements/s avec).
- **« Une famille de 4 541, c'est un bug ? »** → Non : les bases AMP contiennent
  beaucoup de variantes proches d'un même peptide star. La règle “famille
  > 1 % ⇒ train” garantit qu'aucun split ne peut en contenir un morceau.
- **« Pourquoi ne pas exclure les toxines ? »** → Elles sont marquées
  (flagged) et **exclues du génératif** ; `unknown` est traité comme
  flagged (échec sûr) ; la décision d'exclusion finale reste ouverte (§6),
  c'est un choix à discourir, pas à enterrer.
- **« F(0,70) = 0,333, si le split fuit ? »** → Non : c'est la borne basse
  mesurée sur un stock réduit de candidat ; la mesure exacte (tous les
  candidats 0,70) tombe ce soir. La preuve d'étanchéité stricte est F(0,90)
  = 0 — et 0,70 est un seuil “cousin”, il est NORMAL d'y voir des cousins.
- **« Vous avez utilisé scikit-learn ? »** → Non ; tout écrit à la main
  (NumPy/SciPy/pandas + un peu de PyTorch pour les tenseurs) ; cross-check
  possible en mois 3 si demandé.

---
*Chiffres collés depuis : PR #21 (`feat/sofiene-a5-a6-v02`), PR #20 (fusionnée),
`splits_manifest.json`, `PROVENANCE.md`.*
