# Ajustements Semaine 2 — réconciliation avec le plan de M. Abdallah

**Date de rédaction :** dimanche 25/10/2026
**Réunion M. Abdallah :** mardi 27/10/2026
**Documents de référence :**
- Plan de l'encadrant : `Semaine2_peptides_alignement_splits.docx` (tâches T1–T31, 3 étudiants)
- Plan équipe : `Docs/Semaine2_peptides_peplib_4etudiants.md` (tâches T30–T47, 4 étudiants)
- Planning : `Docs/PLANNING_Semaine2.md` (vendredi 23 – mardi 27)

**Objet :** décrire les ajustements à apporter à notre plan pour satisfaire les
demandes de l'encadrant sans perdre le travail déjà réalisé, les assigner,
fixer les échéances et l'ordre de travail.

---

## 1. Constat

Les deux plans décrivent la même semaine (alignement à la main, clusters,
splits étanches, tokenisation, baselines, tag `v0.2`, aucun réseau de neurones).
Le plan de l'encadrant ajoute trois éléments que le nôtre ne couvre pas encore :

1. **DECISIONS.md** — conventions figées et signées *avant* tout calcul (T2).
2. **La fuite comme objet de mesure** — split aléatoire de contrôle (T16) et
   courbe de fuite F(0,70) / F(0,90) (T18), qui répondent à sa « question de la
   semaine » : un split aléatoire surestime-t-il la qualité, et de combien selon
   l'ordre n du modèle à comptages ?
3. **Le critère d'identité de couverture** (alignement local SW normalisé par la
   plus courte séquence) comme critère de redondance, l'identité globale (NW)
   rapportée en complément.

À l'inverse, notre plan est plus large sur les données : APD3 + AMP Challenge
2027, licences GO/NO-GO, recouvrement trilatéral, dédoublonnage inter-sources,
`corpus_manifest.json`. Rien chez l'encadrant ne contredit ces ajouts : le
corpus v0.2 à trois sources est un sur-ensemble de son « union des sources ».
Il faut le présenter comme tel mardi, pas comme une déviation.

**Principe à poser mardi :** « Ton plan est notre plan, l'ordre du plan annuel
est respecté (mois 2 : NW/SW + clusters + splits étanches + baselines) ; nous
l'avons élargi au corpus à 3 sources car l'union des sources est ton critère ;
voici la matrice de correspondance T1–T31 ↔ nos tâches/PR et DECISIONS.md à
confirmer. »

---

## 2. Ce qui est déjà couvert (à assumer tel quel mardi)

| Plan encadrant | Plan équipe | État au 25/10 |
|---|---|---|
| T3 BLOSUM62 saisie | T33 (Youssef) | en cours (dimanche) |
| T4–T6 NW linéaire, Gotoh, SW | T36 (Youssef) | en cours (NW, puis SW dimanche/lundi) |
| T12–T13 union-find, clusters 0,90 | T35 (Sofiene) | fusionné sur PR #16 (Levenshtein, corpus v0.1) |
| T15 split 90/5/5 par clusters, seed 42 | T39 (Sofiene) | prévu lundi (après recalcul des clusters sur v0.2) |
| T20–T21, T23 tokenisation V = 23, round-trip 1 000 | T40 (Youssef) | prévu lundi |
| T19 manifeste des splits | T41 (Moi) | corpus_manifest.json prévu (T19 l'exprime déjà, courte réécriture) |
| T25–T26, T27 baselines, nats/token | T43 (Youssef) | partiel (aléatoire/unigramme/bigramme — à étendre, cf. A7) |
| T31 tag `v0.2`, relance | T46 (Moi) | prévu lundi soir |
| Drapeaux de toxines | toxin_flags.json (semaine 1) | à relancer sur le corpus fusionné (Yassine) |
| Biosécurité, séminaire 30 min | T47 + §8 plan équipe | inchangé |

Les plus-values équipe à marquer « hors tâches T1–T31 » dans la matrice :
**T42** (dédoublonnage par représentants de famille), **T44** (validation des
descripteurs contre APD3 et les en-têtes du challenge), **T45** (témoins
négatifs), recouvrement trilatéral, GO/NO-GO licences. Aucune n'est contredite
par le plan de l'encadrant.

---

## 3. Ajustements à exécuter

Chaque ajustement porte un identifiant **A1–A10**. Un ajustement est
« terminé » quand sa preuve de validation est visible (PR relue, test vert,
figure ou fichier généré) — pas quand le code est écrit.

| ID | Ajustement | Assigné | Échéance | Validation exigée |
|----|------------|---------|----------|-------------------|
| **A1** | **DECISIONS.md** : reprendre la section 3 du plan de l'encadrant (conventions : BLOSUM62, gap `c(k) = o + e·(k−1)` avec o = 11, e = 1, identité de couverture 0,90, clustering par composantes connexes, 3-mers k = 3, split 90/5/5 par clusters ± 0,5 pt avec clusters > 1 % en train, seed 42, V = 23, max_len 52, évaluation en nats/token sur 21 cibles, test lu une seule fois) **et** ajouter nos conventions (stratification par source d'ingestion, corpus v0.2 à 3 sources, filtres semaine 1), avec la mention « proposition — à confirmer mardi ». Mise à jour datée et justifiée si l'encadrant modifie mardi. | Moi (coordinateur) | **dimanche 25/10 soir** | PR ouverte et relue ; l'encadrant peut la lire avant la réunion |
| **A2** | **Demi-page de dérivation du seuil m** (T10 encadrant) : dériver par écrit la borne de 3-mers partagés qui garantit le rappel à 0,90 (l'indice est déjà dans son document : une opération d'édition casse au plus k k-mers), ensuite seulement la valeur concrète. La fonction de préfiltre devient paramétrable par m. Nos mesures empiriques (PR #16) servent de validation. | Sofiene (rédaction + paramétrage dans `peplib/dedup.py`), relue par Youssef | **lundi 26/10 matin** | demi-page + test du rappel 100 % sur le sous-échantillon |
| **A3** | **Alignement des noms de modules** sur la liste de l'encadrant : `peplib/matrices.py`, `align.py`, `cluster.py`, `splits.py`, `tokenizer.py`, `datasets.py`, `ngram.py`, avec `data/splits/`, `splits_manifest.json`, `results/week2/`, `notebooks/NB02_peptides.ipynb`. À la création, utiliser le nom de l'encadrant ; si un module existant doit changer de nom, garder un alias d'import jusqu'à la fin de semaine pour ne pas casser les branches en cours. | Moi (cohérence globale) ; chaque auteur applique à la création de son fichier | à la création (dimanche–lundi) | les chemins attendus existent dans l'arborescence |
| **A4** | **Rapport d'une page + squelette NB02** : conventions, tableau des baselines, courbe de fuite, réponse mesurée à la question de la semaine, limites. Les chiffres viennent de `results/week2/` et du manifeste — rien à recalculer. NB02 appelle `peplib` uniquement (pas de `ids[:max_len]`, aller-retour sur 1 000 séquences). | Moi (rapport), Moi + Tous (NB02) | **lundi soir**, finalisé mardi matin après les figures | rapport d'une page rendu ; NB02 relancé sans erreur |
| **A5** | **Identité de couverture + recalcul des clusters sur v0.2** : Youssef écrit `identite_couverture` (quelques lignes par-dessus SW) et `identite_globale` ; Sofiene recalcule les clusters du corpus v0.2 avec le critère couverture à 0,90. Le clustering Levenshtein v0.1 (PR #16) reste la première version, conservée et rapportée dans DECISIONS.md comme mesure complémentaire — c'est le T14(b) de l'encadrant. **Voir la mise en garde technique du §5 avant tout recalcul.** | Youssef (fonctions d'identité), Sofiene (recalcul) | SW + identités : **dimanche soir** ; recalcul clusters : **lundi matin** | 2 x `identite_*` testées ; tableau clusters v0.2 nb/max/singletons/fraction du plus grand (T13) |
| **A6** | **`peplib/splits.py` + split aléatoire de contrôle + courbe de fuite** : split par clusters 90/5/5 (déjà prévu T39) **plus** le split aléatoire de contrôle (T16) de mêmes tailles ; courbe de fuite (T18) : pour chaque séquence de val, meilleure identité de couverture avec train parmi les candidats du préfiltre ; fractions F(0,90) et F(0,70) pour les deux splits. Figure dans `results/week2/`. | Sofiene | **lundi 26/10** (split par clusters l'après-midi, contrôle + courbe le soir) | F(0,90) = 0 pour le split par clusters ; valeur rapportée pour le split aléatoire ; mêmes tailles de jeux |
| **A7** | **Baselines étendues à l'ordre 3 avec α choisi sur val** (T26) : uniforme + ordres 0, 1, 2, 3, entropie croisée en nats/token sur 21 cibles (20 AA + EOS, BOS et PAD jamais cibles), lissage additif avec α ∈ {0,01 ; 0,1 ; 1} choisi sur val ; tableau final évalué sur le split par clusters **et** le split aléatoire (T27), figure de l'écart en fonction de n. | Youssef | **lundi soir ; repli : ordres 0–2 mardi, ordre 3 jeudi matin, annoncé tel quel à l'encadrant** | probabilités de chaque contexte sommant à 1 ; uniforme = ln 21 ≈ 3,0445 ; test lu une seule fois |
| **A8** | **`peplib/datasets.py`** (T22) : Dataset PyTorch, collate_fn avec padding à droite + masque + longueurs, DataLoader avec générateur seedé, **aucune boucle d'entraînement** (conforme aux outils autorisés de l'encadrant et au plan annuel). | Youssef | **lundi si le temps ; sinon différé jeudi avec raison écrite (critère prévu dans le plan de l'encadrant)** | batchs de forme (B, T_max), PAD seulement à droite ; même seed, mêmes batchs |
| **A9** | **Vérification exhaustive + budget de temps** (T8, T9) : énumérer tous les alignements de séquences de longueur ≤ 4, comparer NW/SW à l'optimum énuméré sur ≥ 500 paires (seed fixée) ; temps moyen d'alignement pour longueurs 10, 30, 50 + extrapolation du coût tous-contre-tous sur le corpus ; budgets écrits dans DECISIONS.md. **Vérification faite par quelqu'un d'autre que l'auteur des alignements.** | Sofiene (vérification exhaustive des alignements de Youssef), Youssef (tableau temps), Moi (budgets consignés) | **lundi** | pytest vert ; tableau temps + budgets dans DECISIONS.md |
| **A10** | **Matrice de correspondance T1–T31 ↔ tâches équipe ↔ PR** : une page qui montre à l'encadrant que chaque tâche de son plan a un propriétaire et un état chez nous, plus la ligne « hors T1–T31 » pour T42/T44/T45 et le corpus v0.2. | Moi | **dimanche soir** (avec A1) | la page existe, chiffres à jour au lundi soir |

---

## 4. Ordre de travail et échéances (vue calendaire)

### Dimanche 25/10 (ce soir)

| Ordre | Qui | Tâche | Échéance |
|-------|-----|-------|----------|
| 1 | Moi | A1 DECISIONS.md + A10 matrice de correspondance | ce soir |
| 2 | Moi | Clôture S1 : merge PR #15 puis PR #16 (revue approuvante + dérogation SonarCloud acteé), commit des documents Semaine 2, tag `v0.1` | ce soir ou demain matin — **avant** tout merge v0.2 |
| 3 | Youssef | Fin de SW (T6) puis `identite_couverture` / `identite_globale` (T7) — prérequis de A5 | ce soir |
| 4 | Yassine | Fin de T37/T38 (corpus v0.2, recouvrement trilatéral) — déjà au planning | ce soir |
| 5 | Sofiene | A2 : début de la demi-page du seuil m + paramétrage du préfiltre | ce soir (fin demain matin) |

### Lundi 26/10

| Ordre | Qui | Tâche | Dépend de | Échéance |
|-------|-----|-------|-----------|----------|
| 6 | Sofiene | A5 : recalcul des clusters v0.2 au critère couverture 0,90 (+ T14 : à 0,70, à identité globale, (o,e) ∈ {(11,1),(10,1),(8,2)}) | A5 (SW + identités) | matin |
| 7 | Sofiene | A2 : demi-page seuil m terminée + rappel mesuré | préfiguré ce soir | matin |
| 8 | Youssef | A9 : vérification exhaustive des alignements de Youssef (Sofiene vérifie, Youssef mesure : rôles inversés, 1 fois) | T5–T7 | matin |
| 9 | Sofiene | A6 : `splits.py` — split par clusters 90/5/5 stratifié par source, puis split aléatoire de contrôle | 6 | après-midi |
| 10 | Youssef | `tokenizer.py` (T20–T21) puis `ngram.py` ordres 0–3, α sur val (A7) | corpus v0.2 figé | après-midi → soir |
| 11 | Yassine | PROVENANCE.md v0.2 (2 licences + 13 bases amont), vues par source, relance toxin_flags sur v0.2 | merge v0.2 | après-midi |
| 12 | Sofiene | A6 : courbe de fuite (val × train, deux splits) | 9 | soir |
| 13 | Moi | A4 : rapport une page + squelette NB02 ; merge du jour figures PR ; tag `v0.2` après relance complète | 6–12 | soir |
| 14 | Youssef | A8 datasets.py si le temps | 10 | soir (optionnel) |

### Mardi 27/10 (réunion M. Abdallah)

| Ordre | Qui | Tâche | Échéance |
|-------|-----|-------|----------|
| 15 | Moi | Générer les 3 figures de garde-fou : distribution des tailles de familles, courbe de fuite (clusters vs aléatoire), échelle des baselines — depuis les fichiers du lundi, **sans recalcul** | mardi matin |
| 16 | Moi | Finale : test **lu une seule fois par configuration**, tableau final (T27) — les configurations doivent être **figées lundi soir**, pas mardi (consigne explicite de l'encadrant) | mardi matin |
| 17 | Tous | Réunion : pitch en 3 phrases (§1), matrice A10, DECISIONS.md à confirmer, visuels | réunion |

### Jeudi 29/10 (repli si différé)

| Qui | Tâche |
|-----|-------|
| Youssef | ordre 3 des baselines et/ou `datasets.py` si différé, avec raison écrite dans DECISIONS.md |

**Invariants pendant toute l'exécution :** seed 42 partout ; personne ne
travaille sur `main` ; une PR mergée exige un relais approuvant + `pytest` vert
(Python 3.11 et 3.14) ; UTF-8 sans BOM ; le tag `v0.1` doit exister **avant**
le tag `v0.2` ; aucun merge de corpus sans clusters testés (sans conflit).

---

## 5. Mise en garde technique — identité de couverture vs Levenshtein

- Un fragment exact de 20 résidus d'un peptide de 50 a une identité globale
  Levenshtein de 0,4 (sous le seuil 0,90, la paire passe le test) mais une
  **identité de couverture de 1,0** : fuite flagrante. C'est exactement le motif
  écrit par l'encadrant (« un fragment exact d'un peptide d'entraînement est
  une fuite que l'identité globale ne détecte pas »). Les deux mesures
  cohabitent et la différence est documentée dans DECISIONS.md — on ne jette
  pas le travail de la PR #16, on le garde comme mesure complémentaire.
- **Ne pas réutiliser la formule `max_distance((1−τ)·max(|a|,|b|))` de la PR #16
  telle quelle pour le critère de couverture.** La borne exacte au seuil τ est
  différente pour la couverture SW ; elle doit être re-démontrée avec le seuil
  m des k-mers partagés (c'est précisément l'objet de la demi-page A2). Toute
  borne approximée non justifiée contredirait la règle du projet : pas
  d'epsilon numérique sans démonstration.
- Le recalcul des clusters porte sur **v0.2** (les clusters v0.1 ne couvrent que
  9 021 séquences sur les ~10⁴–10⁵ du corpus fusionné ; règle T39 déjà écrite
  dans la docstring du module et dans PROVENANCE.md, qui rappelle que les
  clusters versionnés v0.1 doivent être recalculés sur tout nouveau corpus).

---

## 6. Arbitrages à trancher en équipe avant mardi (décision d'équipe, pas du coordinateur seul)

1. **Division du travail 3 → 4 :** l'encadrant répartit T1–T31 sur 3 rôles
   (E1 données, E2 modèles, E3 évaluation). Proposer mardi : E1 ≈ Yassine,
   E2 ≈ Youssef, E3 ≈ Sofiene (il vérifie le clustering et la fuite — tâches E1
   vérifiées par E3, conformément à sa règle « l'auteur ne valide pas sa
   brique »), coordinateur = revue de tous les PR + merges. À valider par les
   quatre, puis présenter à l'encadrant.
2. **Noms de modules :** suivre les noms de l'encadrant dès la création
   (recommandation du coordinateur : oui — friction gratuite évitée), ou garder
   les noms actuels en documentant la correspondance dans le README.
3. **Ordre 3 et `datasets.py` avant ou après mardi :** si après, l'annoncer tel
   quel à l'encadrant (« ordres 0–2 mardi, suite jeudi ») — son propre critère
   accepte le report avec raison écrite.
4. **Langue des commits et commentaires de PR** (français vs anglais) :
   toujours en suspens ; à trancher maintenant, toute prose écrite cette
   semaine en dépendra.

---

## 7. Rappels pour la réunion

- L'objectif de l'encadrant, cité dans son message : « apprendre les techniques
  informatiques et AI adapté au design des peptides thérapeutiques » — la
  courbe de fuite et le tableau des baselines répondent directement à sa
  question empirique ; c'est la démonstration la plus convaincante à montrer.
- Montrer des visuels et documents, jamais du code (demande explicite) :
  figures `results/week2/`, DECISIONS.md, matrice de correspondance, rapport
  d'une page, PROVENANCE.md, manifeste.
- Savoir expliquer sans le code : pourquoi le cluster protège l'examen, pourquoi
  la stratification par source colle des peptides AMP dans le test, ce que
  signifie nats/token.
- La question JIT/NumPy en anti-diagonales (T9) est explicitement une décision
  de l'encadrant : la poser mardi au lieu de trancher seuls.
