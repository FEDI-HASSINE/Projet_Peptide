# RAPPORT CRITIQUE — PROJET PEPTIDES (SEMAINE 3 + PROJECTION)

> **Objet revu** : `main` @ `df7f493`, tag `v0.3-semaine3` (commit de code `362c7e5`), semaine 3 du plan 12 mois.
> **Date de revue** : 23 novembre 2026.
> **Méthode** : lecture intégrale du code `peplib/`, des configs, des splits et des artefacts ; exécution de la suite de tests ; **expériences de contrôle refaites par mes soins** (baselines n-gram, baseline positionnelle, comptage de paramètres, latence de génération, vérification de l'initialisation, audit des exemples publiés, ré-entraînement LSTM seed 42). Les chiffres cités sans mention « mesuré » proviennent du dépôt ; ceux marqués **mesuré** ont été recalculés ici, sur les mêmes splits, avec le code du dépôt.
> **Périmètre** : ce qui est dans le dépôt. Ce qui n'y est pas (checkpoints, `runs.csv`, script des figures) est traité comme non prouvé, pas comme acquis.

---

## 0. RÉSUMÉ EXÉCUTIF (5 lignes pour un CTO)

1. L'ingénierie est propre pour un projet étudiant (689 tests verts, seed, config, splits étanches à 90 %) — **mais 0 test ne couvre la boucle d'entraînement**, le seul composant qui produit les chiffres.
2. **Aucun artefact de la semaine 3 n'est reproductible depuis le dépôt** : pas de `runs.csv`, pas de checkpoint, aucune entrée `PROVENANCE.md`, et les deux figures publiées (courbe de perte, exemples générés) ne sont produites par aucun script du dépôt.
3. **La métrique phare est plausible, mais sa figure est fausse** : mon ré-entraînement (seed 42, config publiée) atteint 15,38 en 6 époques, cohérent avec 15,21 ; en revanche la courbe de validation publiée finit à `val_loss ≈ 2.06` (PP ≈ 7,9), ce qui n'est pas ce run. Le chiffre tient, la preuve graphique non.
4. **2 des 4 peptides « générés » publiés sont des séquences du train** (`GLFDIVKKVVGALG` à l'identique, `KWKLFKKIEKVGQN` = préfixe exact d'une séquence du train). Aucun contrôle de nouveauté n'a été exécuté, alors que `peplib.metrics_proximity` existe et n'est jamais appelé hors tests.
5. **La génération des modèles entraînés n'existe pas** : `generate_sequences` n'est appelé que dans les tests, sur des modèles à initialisation aléatoire. « 100 % valides » est une tautologie du rejet d'échantillons, pas un résultat.

**Verdict** : `GO conditionnel` pour le mois 4, avec **8 correctifs P0 hors modélisation (≈ 4-5 jours-homme)** à passer avant d'ouvrir la comparaison D3PM. Détail en §5.

---

## 1. CRITIQUE TECHNIQUE SEMAINE 3

### 1.0 Ce que j'ai vérifié moi-même (base factuelle)

| Vérification | Résultat mesuré |
|---|---|
| `pytest -q` | **689 passed en 72 s** — le chiffre annoncé est exact |
| Tests touchant `train_model` / `_run_epoch_*` | **0** (`grep` sur `tests/`) |
| `generate_sequences` sur un modèle entraîné, hors tests | **0 occurrence** dans le dépôt |
| Paramètres LSTM (config publiée) | **930 455** |
| Paramètres Transformer (config publiée) | **802 967** (+15,9 % pour le LSTM) |
| Écart-type de l'embedding après `_init_weights` | **0,115** (attendu 0,02) ; norme L1 de la ligne PAD = **12,79** (attendu 0) |
| Latence de génération CPU, batch=1 | LSTM **307 ms/séquence** ; Transformer **1 564 ms/séquence** |
| Perplexité val des baselines recalculées | unigramme add-1 **19,49** · bigramme **19,27** · trigramme **19,53** · interpolation Jelinek-Mercer ordre 2 **19,19** · unigramme positionnel **18,65** · « position+AA précédent » **18,74** |
| Ré-entraînement LSTM seed 42 (config publiée, run interrompu après l'époque 6) | val_pp époques 1→6 : **16,96 / 16,31 / 15,54 / 15,39 / 15,55 / 15,38** → le 15,21 publié est cohérent ; la courbe de validation de la figure (2,69 → 2,06, soit PP 7,9) ne l'est pas |
| Composition du corpus (train) | AMPChallenge **73,6 %** · UniProt 17,5 % · APD3 1,4 % · recouvrements 7,6 % · BACTIBASE/BIOPEP **0 %** |
| Longueur train / val / test | 22,81 / 21,24 / 20,91 AA (médiane 20 partout) |
| Toxines flaguées | train **1 493 (3,9 %)** · val 161 (**7,7 %**) · test 132 (6,2 %) ; **79 % du corpus a un statut « unknown »** |
| Chevauchement de clusters train/val et train/test | **0** (confirmé) — mais `F_70 = 0,607` sur le split familles (mesure de la semaine 2) |

---

### 1.1 Architecture

| # | Point | Verdict | Constat et justification |
|---|---|---|---|
| 1 | Masque causal du Transformer | ✅ | Deux niveaux de preuve concordent : le test comportemental `tests/test_transformer.py:72` (modifier le token futur ne change pas les logits passés) et ma sonde directe (batch paddé, 4 couches : 0 NaN dans les logits, 0 paramètre avec gradient NaN). Le piège classique du `TransformerEncoder` PyTorch — lignes entièrement masquées → `softmax(-inf) = NaN` qui contamine les gradients via `LayerNorm.backward` — **n'est pas déclenché ici** (vérifié : 0 NaN dans les logits, 0 paramètre avec gradient NaN sur un batch où 15 positions sur 40 sont du padding). Attention : ce comportement dépend du *backend* d'attention choisi par PyTorch ; la garantie doit rester le test comportemental, pas une propriété supposée du framework. |
| 2 | Weight tying (Transformer) | ❌ | `transformer.py:103` lie `head.weight = token_emb.weight`, puis `_init_weights()` (`:119-126`) parcourt `self.modules()`, trouve `nn.Linear` (la tête) et applique `xavier_uniform_` **sur le tenseur partagé** : l'initialisation `N(0, 0,02)` et la remise à zéro de PAD sont écrasées. Mesuré après construction : `token_emb.std = 0,115` au lieu de 0,02, `||PAD||₁ = 12,79` au lieu de 0. Conséquence : (a) l'init « small embeddings » recommandée pour le tying (Press & Wolf 2017) est perdue, (b) la ligne PAD n'est plus nulle et **dérive** car le gradient de la tête s'écrit dans le même tenseur alors que `padding_idx` ne bloque que le gradient d'embedding. Le docstring affirme l'inverse de ce que fait le code. |
| 3 | Weight tying (LSTM) | ❌ | Le tableau de livrables annonce « weight tying » pour le LSTM. C'est **faux** : `configs/lstm_base.yaml:11` pose `tie_weights: false`, `train.py` ne lit jamais cette clé, et les dimensions sont incompatibles (embed 23×128 vs tête 256×23) sans projection. Le rapport `RAPPORT_SEMAINE3.md` reprend cette affirmation. Erreur documentaire, pas neutre : elle sert d'argument de qualité. |
| 4 | Capacité et comparabilité des deux références | ⚠️ | 930 455 paramètres (LSTM) vs 802 967 (Transformer), dropout 0,2 vs 0,1, 2×256 vs 4×128/FFN 512. **Le LSTM gagne avec 16 % de paramètres en plus et deux fois plus de dropout.** Aucune ablation, aucun contrôle de budget. Pour un projet dont la thèse est « benchmark 4 familles à budget égal », la comparaison des deux références est déjà non contrôlée — et c'est celle qui sera citée. |
| 5 | Architecture Transformer | ⚠️ | Encodage positionnel appris, `norm_first=False` (post-LN). Le commentaire `transformer.py:92` (« norm_first désactivé : incompatible avec enable_nested_tensor ») est **factuellement exact** (PyTorch désactive le *fast path* nested tensor si `norm_first=True`, vérifié dans le source de `TransformerEncoder.__init__`) — mais **sans objet ici** : ce fast path exige `mask is None` et un padding mask left-aligned, et il est désactivé en mode entraînement (`TransformerEncoder.forward`, conditions `mask is not None` et `first_layer.training`). Or le modèle passe **toujours** un masque causal dense (52×52). Pre-LN était donc disponible à coût nul, et c'est la variante stable sans warmup long. Le choix post-LN n'est ni documenté comme un choix ni testé. |
| 6 | Qualité de la référence LSTM | ⚠️ | Dropout appliqué trois fois (après embedding, dans le LSTM, après le LSTM) à p=0,2 sur 38 k séquences. Plausible, mais **c'est un hyperparamètre non exploré** : le « LSTM > Transformer » peut être un artefact de régularisation, pas d'architecture. Sans sweep (même grossier, 3 valeurs), la conclusion n'est pas défendable en review. |
| 7 | Coût de décodage | ❌ | Pas de KV-cache, pas de génération par batch : `generate.py:96-101` refait un forward complet sur tout le préfixe à chaque token, une séquence à la fois. Mesuré : **1,56 s/séquence** pour le Transformer sur CPU (LSTM 0,31 s). Pour 10⁴ séquences : 4,3 h contre 51 min. La grille d'échantillonnage publiée (3 températures × 2 top-k × 2 top-p = 12 configs) × 10⁴ séquences ≈ **100 h CPU par famille**. C'est un bloqueur du mois 8, pas un détail. |
| 8 | Absence de *seam* pour le conditionnement | ⚠️ | `generate_sequences(model, n, temperature, top_k, top_p, max_len, seed, device)` n'a ni préfixe, ni attribut, ni log-probabilité de sortie. Le mois 9 (génération conditionnée + Pareto) devra réécrire l'API ; l'anticiper coûte 20 lignes maintenant. |

**Note architecture : ⚠️ (6,5/10).** Le cœur (causalité, masquage) est correct et testé ; deux bugs réels (init liée, légende LSTM) et un choix de décodage qui ne passe pas l'échelle.

---

### 1.2 Entraînement

| # | Point | Verdict | Constat et justification |
|---|---|---|---|
| 1 | Padding | ❌ | `train.py:102-103` instancie `PeptideDataset(..., pad_to=52)` : **toutes** les séquences sont paddées à 52 tokens, alors que la longueur moyenne est de 24,8 tokens (22,81 résidus + BOS/EOS). **2,10× de calcul jeté**, et 4× sur le coût d'attention du Transformer. Le chargeur dynamique existe (`datasets.py:263`, `collate_peptides` + `make_loader_s3`) et **n'est pas utilisé par le trainer**. Le rapport présente le padding dynamique comme livré : il l'est, mais pas branché. |
| 2 | Validation de config | ❌ | `peplib/config.py:39` implémente une validation stricte (clés inconnues refusées, `vocab_size=23`, `max_len=50`, `d_model % n_heads`) — et `train.py` **ne l'appelle jamais** : il fait `yaml.safe_load` (`train.py:39-41`) et lit les clés avec `_cfg(..., default=...)`. Une faute de frappe (`d_hiden`, `wight_decay`) est silencieusement remplacée par la valeur par défaut. Les 689 tests valident une fonction que le chemin d'entraînement n'utilise pas : **fausse assurance**. |
| 3 | Clés de config mortes | ⚠️ | `tie_weights`, `train.seed`, `train.optimizer`, `early_stopping.metric`, `sampling_grid`, `warmup_steps` (LSTM) : lues par personne. Mettre `optimizer: sgd` ne change rien. Une config qui ment est pire qu'une config absente. |
| 4 | AdamW / weight decay | ⚠️ | `train.py:284` applique `weight_decay=0,01` à **tous** les paramètres, y compris LayerNorm et biais (init LN à 1,0 → contraction), et les embeddings. Standard : deux groupes, `wd=0` sur biais/LN/embeddings. Impact faible mais gratuit à corriger. |
| 5 | Politique de learning rate | ⚠️ | LSTM : 1e-3 constant, aucun scheduler. Transformer : warmup linéaire 200 pas (`0,67` époque) puis **constant**, sans décroissance. Combiné à `patience=5`, l'arrêt précoce se déclenche pendant un plateau à LR maximale : on jette des époques qui auraient convergé avec un decay cosinus. C'est le correctif le plus rentable en qualité par ligne de code. |
| 6 | Early stopping | ✅ | Patience 5 sur la loss de validation moyennée par token (~48 k tokens pour la val) : le bruit d'estimation est faible (σ sur la PP finale = 0,06 sur 3 seeds), donc le critère est stable. Correct. |
| 7 | Sélection sur la même val que la métrique publiée | ⚠️ | Le « meilleur modèle » est choisi sur `val_loss`, puis la performance est **rapportée sur cette même val**. Biais de sélection optimiste (min sur ~11 époques), de l'ordre de 0,5-1 % en PP, non quantifié. Il faut soit publier aussi la dernière époque, soit réserver une val de sélection. Le test reste verrouillé : c'est bien. |
| 8 | `--resume` | ❌ | `_load_ckpt` recharge `best.pt` (**pas** le dernier), ne restaure pas l'état du scheduler, remet `no_improve=0` et repart de `best_epoch+1`. Reprendre après un arrêt précoce rejoue donc les époques sans mémoire de la patience et peut boucler. `torch.load` est appelé sans `weights_only=True` (exécution de code arbitraire sur un checkpoint non fiable). Le notebook `NB03` documente `--out-dir`, option qui **n'existe pas** dans `train.py`. |
| 9 | Traçabilité de l'expérience | ❌ | `runs.csv` contient `model, seed, epoch, train_loss, val_loss, val_pp, elapsed_s` — pas de hash de config, pas de SHA git, pas d'horodatage, pas de device, pas de version torch. Et **le fichier n'est pas dans le dépôt** (`results/week3/` ne contient que `baselines_metrics.csv` et `visuals/`). `PROVENANCE.md` n'a **aucune section semaine 3** (la dernière est « Variantes du clustering v0.2 »). Le projet impose un enregistrement de provenance par artefact : il n'est pas fait pour le seul artefact qui compte. |
| 10 | Couverture de test du trainer | ❌ | **Aucun test n'appelle `train_model`** ni `_run_epoch_train`. Pas de test de reprise, pas de test « 1 époque sur données jouets → la loss baisse », pas de test de non-régression sur la PP d'un batch fixe. « 689 tests » couvre tout sauf l'entraînement. |
| 11 | **Cohérence chiffres / courbe** | ❌ | La figure `2_courbes_perte.png` affiche `val_loss = 2,06` à l'époque 11 → **PP = e^2,06 = 7,85**. Le tableau publie **15,21 ± 0,06** → val_loss = 2,72. J'ai rejoué le run (config publiée, seed 42, CPU) : époques 1→6 : **16,96 / 16,31 / 15,54 / 15,39 / 15,55 / 15,38**. Donc **le tableau est crédible** (15,21 = meilleure époque sur ~11, ce que ma trajectoire atteint), et c'est **la série de validation de la figure qui est fausse** : elle n'est pas issue de ce run. Détail aggravant : sur la figure, `train_loss` (2,78 en époque 1) colle à mon run (2,7884), alors que `val_loss` (2,69) ne colle pas (2,8309). Une figure de résultats dont une des deux courbes n'a pas été produite par le code est un défaut de publication, même si le chiffre final est juste. |

**Note entraînement : ❌ (4/10).** La boucle est lisible et correcte algorithmiquement ; elle est non testée, non validée par son propre schéma et non traçable. Le 15,21 est plausible — je retrouve 15,38 à l'époque 6, et la config mène bien vers 15,2 sur ~11 époques — mais aucun artefact du dépôt ne permet de le vérifier, et la figure associée est fausse.

---

### 1.3 Génération

| # | Point | Verdict | Constat et justification |
|---|---|---|---|
| 1 | PAD/BOS à −∞ | ✅ | `sampling.py:65-66` avant température, avec garde-fou NaN/+inf et vérification qu'il reste un token autorisé. Testé (`test_generation.py` / `test_sampling.py`). Correct. |
| 2 | Température / top-k / top-p | ✅ | Ordre standard (T → k → p), renormalisation après masquage, `top_p` implémenté comme « plus petit ensemble de cumul ≥ p » (`cumsum − p_i ≥ p` → rejet). Équivalent à HF. Validation des types et des bornes. Bon. |
| 3 | Graine de génération sur GPU | ❌ | `generate.py:110` : `generator=gen if device.type == "cpu" else None`. Sur **T4 — la cible annoncée du projet** — la graine est ignorée et `torch.multinomial` tire dans le RNG global : génération **non reproductible** et non isolée entre appels. Le rapport affirme « `torch.Generator` seedé » sans nuance. Correctif : `torch.Generator(device=logits.device)`, passer `generator` dans tous les cas, ajouter un test GPU. |
| 4 | « 100 % valides » | ❌ | Tautologie. `generate_sequences` rejette en interne tout échantillon invalide et boucle jusqu'à `n` valides (`max(10n, 5000)` tentatives) : par construction, 100 % des séquences **renvoyées** sont valides. La métrique mesurable est le **taux de validité brut = valides / tentatives**, et le taux de rejet par cause ; ni l'un ni l'autre n'est rapporté. C'est exactement le piège que GuacaMol/MOSES imposent d'éviter (rapporter validité *avant* filtrage). |
| 5 | Asymétrie avec les baselines | ⚠️ | Les n-grammes interdisent EOS avant 10 AA (`baselines_generation.py:154`) ; le réseau, non — il est filtré **après** coup. Les colonnes `validity_rate`, `length_mean`, `js_divergence`, `ks_length` ne portent donc pas sur le même protocole. La comparaison publiée n'est pas valide, même si l'écart est favorable au réseau. |
| 6 | Aucune génération des modèles entraînés | ❌ | `generate_sequences` n'est appelé que dans `tests/test_generate.py`, sur des modèles **aléatoires** (embed 16, hidden 32). Aucun script du dépôt n'ouvre un checkpoint et n'écrit un échantillon. `results/week3/visuals/3_exemples_generes.png` n'est produit par **aucun** code versionné. |
| 7 | Exemples publiés | ❌ | Audit des 4 séquences de la figure contre `train_v02.csv` : `GLFDIVKKVVGALG` = **séquence du train, à l'identique** ; `KWKLFKKIEKVGQN` = **préfixe exact (14/38 AA)** d'une séquence du train ; les 2 autres sont absentes du corpus. Donc 50 % des exemples présentés comme « générés » sont des copies (ou des troncatures de copies) du train, sans annotation. Deux hypothèses, toutes deux graves : (a) le modèle mémorise à température 1,0 — et c'est alors le résultat le plus important de la semaine, non mesuré ; (b) les exemples ont été piochés dans le corpus. Dans les deux cas, la publication est trompeuse. `peplib/metrics_proximity.py` existe, est testé, et **n'a jamais été exécuté sur une sortie de modèle**. |
| 8 | Longueur / EOS | ⚠️ | Pas de contrainte de longueur minimale à la génération, pas d'EOS forcé à 50 : la longueur est celle qu'apprend la tête EOS. La distribution de longueur des générés n'est comparée à celle du train que pour les baselines (`baselines_metrics.csv`), pas pour les modèles. |
| 9 | Robustesse de l'API | ⚠️ | Si `n` séquences valides ne sont pas atteintes, la fonction retourne **moins** de `n` sans avertissement ; `model.eval()` n'est pas restauré ; aucun log du taux de rejet. Toute métrique aval hérite silencieusement d'un échantillon conditionné sur le succès. |

**Note génération : ❌ (3,5/10).** L'échantillonneur est correct ; le protocole d'évaluation de la génération n'existe pas, et les preuves publiées ne résistent pas à un audit de 3 minutes.

---

### 1.4 Reproductibilité

| # | Point | Verdict | Constat et justification |
|---|---|---|---|
| 1 | `set_seed` | ✅ | `seed.py` couvre `random`, `numpy`, `torch`, `cuda.manual_seed_all` et **`torch.use_deterministic_algorithms(True)`**. Au-dessus de la moyenne des projets étudiants. |
| 2 | Déterminisme CUDA | ⚠️ | `use_deterministic_algorithms(True)` ne suffit pas : la doc PyTorch exige `CUBLAS_WORKSPACE_CONFIG=:4096:8` (sinon erreur ou non-déterminisme sur les GEMM) et `torch.backends.cudnn.benchmark = False`. Aucun des deux n'est posé. Les runs de `NB03` sont sur T4 : la revendication « bit-à-bit » n'est démontrée **que sur CPU**. |
| 3 | Portée de « bit-à-bit » | ⚠️ | Mon ré-entraînement (même config, seed 42, même machine) donne `train_loss = 2,7884` en époque 1 quand la figure publiée affiche 2,78 : identique à l'arrondi près, **pas bit-à-bit** (threading MKL, ordre de réduction). La revendication doit être reformulée : « même machine, même build, même nombre de threads ». |
| 4 | DataLoader | ✅ | `torch.Generator` seedé sur le shuffle, `num_workers=0` : ordre des batchs reproductible. Testé. |
| 5 | Environnement | ⚠️ | `requirements.txt` épinglé (torch 2.14.1, numpy 2.5.4, pandas 3.0.6) mais **sans hash**, sans Docker/conda-lock, sans `torch.__version__` dans les résultats, sans CI. « Python 3.13 + 3.11 » n'est attesté par aucun artefact exécutable dans le dépôt. |
| 6 | Artefacts | ❌ | Pas de checkpoint (`.gitignore` les exclut — choix défendable), **mais alors il faut committer les hashes, la commande exacte, le device et le temps par époque**. Rien de tout cela n'existe. Ni `runs.csv`, ni entrée PROVENANCE, ni model card. Le contrat du projet (« chaque artefact traçable ») est rompu précisément pour la semaine 3. |

**Note reproductibilité : ⚠️ (5,5/10).** Les mécanismes sont là ; la preuve ne l'est pas.

---

### 1.5 Évaluation

| # | Point | Verdict | Constat et justification |
|---|---|---|---|
| 1 | Métrique unique | ❌ | La perplexité de validation est **la seule** métrique publiée pour les modèles. Elle mesure l'ajustement à la distribution du corpus, pas la qualité de génération. Aucune des métriques minimales d'un benchmark génératif n'est calculée sur les modèles : validité **brute**, unicité à n fixé, nouveauté vs train, Self-BLEU / Distinct-n, correspondance des distributions de propriétés (charge, GRAVY, moment hydrophobe, pI), validité conditionnelle. GuacaMol (Brown et al., 2019) et MOSES (Polykovskiy et al., 2020) ont standardisé ces métriques il y a six ans ; le projet n'en implémente aucune pour ses modèles. |
| 2 | Baseline « barre à battre » | ❌ | 19,49 = unigramme **add-1 sur 21 cibles**. Mesures refaites avec le code du dépôt, mêmes splits : bigramme 19,27, trigramme 19,53, interpolation Jelinek-Mercer 19,19, **unigramme positionnel (P(aa | position), 20 lignes, sans apprentissage) 18,65**. La « barre » publiée est donc **0,84 PP trop haute**. Le gain réel est de 18,4 % par rapport à la meilleure baseline triviale, pas les 22 % suggérés. À noter : `results/week2/baseline_ngram_A7.csv` contient déjà l'ordre 1-3 sur le split familles (19,49 / 19,27 / 19,53) et sur le split aléatoire de contrôle (18,65 / 17,95 / 16,43) ; le tableau de la semaine 3 ne retient que **l'unigramme du split familles**. La comparaison des réseaux aux n-grammes est donc à portée de copier-coller, et elle n'a pas été faite. |
| 3 | Interprétation du gain | ⚠️ | Le fait que toutes les baselines locales plafonnent vers 18,6-19,5 alors que le LSTM atteint 15,2 suggère que le gain vient de l'**inférence de composition globale** (le modèle identifie la classe AMP/cationique dès les premiers résidus et prédit le reste), pas d'une « grammaire » locale. C'est une hypothèse testable en une heure (PP en fonction de la longueur du préfixe ; entraînement sur séquences mélangées comme contrôle négatif). Le rapport affirme « ils ont vraiment appris la langue des peptides » sans le tester. |
| 4 | Seeds et significativité | ⚠️ | 3 seeds, moyenne ± σ : LSTM 15,21 ± 0,06 vs Transformer 15,41 ± 0,16. Écart 0,20 pour une σ poolée ≈ 0,12 : **non significatif** avec n=3 et sans test apparié. Le rapport ne classe pas, mais la figure invite à conclure « LSTM meilleur ». Il faut 5 seeds et un intervalle bootstrap, ou au minimum un test de Welch avec la réserve explicite. |
| 5 | Décomposition | ❌ | PP globale seulement. Pas de PP par tranche de longueur, par source (AMPChallenge vs UniProt vs APD3), ni par statut toxine. Avec 73,6 % d'AMPChallenge courts, la PP globale est dominée par un seul dataset ; un modèle excellent sur les AMPs courts et mauvais sur UniProt peut afficher 15,2. |
| 6 | Métriques de génération des baselines | ⚠️ | `baselines_metrics.csv` compare random / unigramme / bigramme / val / val-mélangée — **sans les modèles**. De plus, la baseline « random » reçoit **la longueur exacte de chaque séquence de référence** (`baselines_generation.py:60`) : son `KS_length = 0` est vrai par construction et son `JS` n'est pas comparable. L'unicité à n = 2095 dans un espace de 20²⁰ séquences vaut 1,0 par construction : la colonne est inutile. |
| 7 | Test set | ✅ | `test_v02.csv` n'est lu par aucun script de la semaine 3 (`evaluate_baselines_s3.py` ne lit que train/val). Le verrou est respecté. C'est un vrai point fort, rare. |
| 8 | Absence de mesure de fuite sur le généré | ❌ | La semaine 2 a produit une courbe de fuite, `F_70 = 0,607` sur le split familles — autrement dit **60,7 % des séquences de val ont un homologue ≥ 70 % d'identité dans le train**. « `F(0,90) = 0` » ne dit rien de la mémorisation, qui se produit à des seuils bien plus bas (les analyses de ProGen/ProtGPT2 utilisent 50-70 %). Appliquer `proximity_summary` aux générés est indispensable et déjà outillé. |

**Note évaluation : ❌ (3/10).** C'est le point faible majeur du projet : la partie « données » est mature, la partie « mesure de ce que le modèle produit » n'existe pas.

---

### 1.6 Tests

| # | Point | Verdict | Constat et justification |
|---|---|---|---|
| 1 | Volume et discipline | ✅ | 689 tests, 72 s, verts sur ma machine. Tests de causalité comportementale, de masque de padding, de seed de DataLoader, de sampling, de config, d'artefacts. Pour un projet à 4 étudiants en semaine 3, c'est nettement au-dessus de la moyenne. |
| 2 | Le trou : l'entraînement | ❌ | Zéro test sur `train_model`, `_run_epoch_train`, `_run_epoch_val`, `_save_ckpt`, `_load_ckpt`. Le composant qui produit les chiffres n'est pas testé. Un test d'intégration de 30 s (200 pas sur 2 000 séquences, assertion `val_pp < 20` + hash de la loss finale) aurait attrapé le bug de config et le bug de padding. |
| 3 | Tests aveugles aux bugs trouvés | ❌ | Trois bugs réels passent la suite : (a) init liée écrasée — aucun test n'asserte `token_emb.weight[PAD].abs().sum() == 0` après construction ; (b) générateur GPU ignoré — les tests ne tournent que sur CPU et n'assertent que « la même seed donne le même résultat » sur CPU ; (c) `train.py` n'utilise pas `load_config` — le test de config passe, le trainer l'ignore. |
| 4 | Tests faibles | ⚠️ | `test_causal_mask_generated_correctly` asserte `mask[i,j] == 0.0 or not mask[i,j].isinf()` : toujours vrai pour toute valeur finie. `test_generate` asserte la validité des sorties — donc **codifie la tautologie** du §1.3.4. `test_count_parameters` borne entre 1e3 et 1e7 : inutile. |
| 5 | Ce qui manque | ⚠️ | Pas de test de reprise (`--resume`), pas de test de génération conditionnée (n'existe pas), pas de test de latence/mémoire, pas de test de non-régression inter-versions (golden file), pas de test de biosécurité sur les sorties générées, pas de CI (aucun workflow GitHub Actions : rien n'empêche de merger du rouge). |

**Note tests : ⚠️ (6/10).** Bonne couverture périphérique, angle mort exactement là où sont les risques.

---

### 1.7 Qualité de code

| # | Point | Verdict | Constat et justification |
|---|---|---|---|
| 1 | Structure, typage, docstrings | ✅ | Modules courts, responsabilité unique, `from __future__ import annotations`, type hints quasi systématiques, docstrings numpy-style. `peplib/` est lisible par un nouveau venu en 30 minutes. |
| 2 | Duplication | ⚠️ | **Quatre** constructeurs de loader (`load_splits`, `make_loader`, `make_loader_s3`, `_make_loaders` dans `train.py`), un seul utilisé. `metrics_generation.is_valid_peptide` réimplémente `peplib.filters.classify_sequence` avec des règles légèrement différentes (alphabet seul vs alphabet + longueur) : deux sources de vérité pour la notion de « valide ». |
| 3 | Code mort / aspirations | ⚠️ | `baseline_a7.py` s'annonce « stub de livraison… sera rempli » et contient une implémentation complète ; `config.py` inutilisé par le trainer ; `tie_weights` fantôme ; `--out-dir` documenté et inexistant. Un relecteur ne peut pas distinguer ce qui est opérationnel de ce qui est décoratif. |
| 4 | Encodage / prose | ⚠️ | Moebake cp1252 visible dans `PROVENANCE.md`, `Docs/clustering_justification_t25.md`, `figures.py` ; les `print()` sont volontairement sans accents, les docstrings non — hétérogène. Pour des livrables académiques, l'encodage doit être UTF-8 partout, accents compris. |
| 5 | Logging / observabilité | ⚠️ | `print` uniquement, pas de `logging`, pas de barre de progression, pas de sortie structurée (JSONL). `runs.csv` est appendé sans verrou ni en-tête versionné. |
| 6 | Hygiène | ⚠️ | `torch.load` sans `weights_only` ; **`AGENTS.md` — le fichier d'instructions et de contraintes du projet — est listé dans `.gitignore` (`git check-ignore` le confirme) : il n'est pas versionné, donc absent d'un clone frais** ; `.pytest_cache` présent dans l'arbre. |
| 7 | Config-driven | ⚠️ | L'intention est bonne (YAML, schéma), l'exécution est partielle : la moitié des clés ne sont pas lues, et le validateur n'est pas sur le chemin critique. |

**Note code : ⚠️ (6,5/10).**

---

### 1.8 Synthèse de notation + Top 3 des correctifs critiques

| Dimension | Note | Poids | Motif principal |
|---|---|---|---|
| Architecture | 6,5/10 ⚠️ | 15 % | Causalité correcte ; init liée cassée ; LSTM non lié annoncé lié ; décodage non scalable |
| Entraînement | 4,0/10 ❌ | 20 % | Padding fixe 2,1×, config non validée, resume cassé, zéro test, chiffres non reproductibles |
| Génération | 3,5/10 ❌ | 15 % | Seed GPU ignorée, validité tautologique, aucun échantillon de modèle entraîné, 2/4 exemples copiés du train |
| Reproductibilité | 5,5/10 ⚠️ | 15 % | Mécanismes présents, preuves absentes (ni runs.csv, ni provenance, ni hash) |
| Évaluation | 3,0/10 ❌ | 20 % | PP seule, baseline « barre » surévaluée de 0,84 PP, aucune métrique générative |
| Tests | 6,0/10 ⚠️ | 5 % | 689 verts mais 0 sur le trainer |
| Qualité de code | 6,5/10 ⚠️ | 10 % | Lisible, mais dupliqué et partiellement décoratif |
| **Global pondéré** | **4,5/10** | | |

**Top 3 des correctifs critiques (avant tout travail sur D3PM)**

1. **Reconstruire la chaîne de preuve** (½ journée). Committer `runs.csv` complet (avec SHA git, hash de config, device, version torch, horodatage), les hashes SHA-256 des checkpoints, une section PROVENANCE semaine 3 (commande exacte, seed, matériel, date), et **corriger la figure de perte** (la série de validation publiée n'est pas celle du run : 2,69 en époque 1 contre 2,8309 mesuré). Tant que ce point n'est pas réglé, tout ce qui suit est invérifiable.
2. **Faire exister la génération des modèles entraînés** (1 journée). Écrire `scripts/generate_week3.py` : charge un checkpoint, génère *N* séquences **sans filtrage**, logge `n_attempts`, `n_valid_raw`, `n_rejected` par cause, puis applique `proximity_summary` (seuils 0,5/0,7/0,9) et les métriques de diversité. C'est aussi ce qui dira si `GLFDIVKKVVGALG` est de la mémorisation ou une copie manuelle.
3. **Corriger et tester les deux bugs silencieux** (½ journée) : `_init_weights` qui écrase l'embedding lié (asserter `PAD` nul et `std ≈ 0,02`), et `generator` ignoré hors CPU (asserter la reproductibilité sur le device courant). Plus `torch.backends.cudnn.benchmark = False` et `CUBLAS_WORKSPACE_CONFIG`.

---

## 2. PROJECTION MONDE RÉEL

### 2.1 Tableau par livrable

| Livrable | Utilité réelle pour un labo | Gap vers production | Intégration biologiste | Conformité / audit | Échelle | Maintenance |
|---|---|---|---|---|---|---|
| **LSTM causal** | Faible en soi : un LM de 0,9 M paramètres sur 38 k séquences ne battra jamais ProGen2 (des milliards de paramètres, des centaines de millions de séquences). Utile comme **contrôle négatif** : établit ce qu'une architecture séquentielle simple extrait du corpus. Valeur réelle = la référence interne du benchmark, pas un générateur de candidats. | Aucun mécanisme de packaging (pas de `save/load` documenté, pas d'export), pas de versionnage de modèle, pas de test sur données hors distribution, pas de mesure de calibration. | Un biologiste ne peut pas l'utiliser : il faudrait un checkpoint, une CLI qui prend `--n 100 --seed 42` et écrit un FASTA avec charge/GRAVY/pI par séquence et un champ `generated_by`. Rien de tout cela n'existe aujourd'hui. | Model card absent, data card absent, pas de provenance par séquence générée (impossible de dire quel checkpoint/seed/config a produit une séquence — bloquant pour un audit). Données publiques : pas de problème RGPD. | CPU seulement démontré ; 0,31 s/séquence en batch 1 (mesuré) → 10⁴ séquences en 51 min, 12 configs en 10 h. Pas de batch, pas de KV-cache, pas de quantification. | Réentraînement manuel, aucun registre de modèles, aucune détection de dérive de distribution. |
| **Transformer causal** | Idem, plus coûteux : 1,56 s/séquence en batch 1 (mesuré), soit **5× le LSTM** pour une PP légèrement pire. Aucun avantage démontré à ce stade, mais la famille est nécessaire au benchmark (et au conditionnement du mois 9). | Pas de KV-cache / `past_key_values` : l'inférence interactive ou par lot est hors de portée en l'état. Pas d'export ONNX/TorchScript. | Même absence d'interface. | Idem. | 4,3 h pour 10⁴ séquences CPU ; sur T4, ~10-20 min avec un décodeur batché + cache — c'est le correctif qui débloque l'échelle. | Idem. |
| **Boucle d'entraînement** | C'est le livrable le plus « réutilisable » : AdamW, early stop, checkpoints, CSV. Un labo peut effectivement l'utiliser comme *baseline harness*. | Manque : param groups, scheduler cosinus, reprise correcte, validation de config sur le chemin critique, logging structuré (JSONL), tracking d'expériences (MLflow/W&B), déterminisme CUDA complet, tests d'intégration. | Indirect (outil interne ML, pas biologiste). | Pas d'audit trail (ni git SHA, ni hash de config, ni version de données dans les logs) → **non conforme aux exigences de reproductibilité d'un dossier réglementaire** (FDA *AI/ML in Drug Development*, 2023 ; EMA *Reflection paper on AI*, 2024). | Mono-nœud, mono-device, pas de DDP/AMP/grad accumulation. Suffisant pour 10⁴ séquences, insuffisant dès 10⁶. | Aucun mécanisme de réentraînement automatisé, de comparaison de versions, de rollback. |
| **Génération (`generate` + `sampling`)** | Cœur de valeur : c'est ce qui produit des candidats. Correct sur le plan algorithmique. | **Seed ignorée sur GPU** (bug), pas de batch, pas de cache, pas de sortie structurée (FASTA+CSV+logprob), pas de taux de rejet, pas de garde-fou de longueur, pas de filtre de biosécurité à la sortie, pas de déduplication. | Doit devenir : `generate --checkpoint X --n 10000 --temp 1.0 --out candidates.fasta --meta candidates.csv` avec colonnes `sequence, charge, gravy, pI, boman, len, best_identity_to_train, model, seed, config_hash`. Aujourd'hui : une fonction Python qui retourne une `list[str]`. | Aucune métadonnée par séquence, donc **traçabilité nulle** d'un candidat jusqu'au modèle. | Voir latence ci-dessus ; bloquant à 10⁴ et au-delà. | — |
| **Configs YAML + schéma** | Bonne pratique ; facilite la reproductibilité et la revue. | Schéma non branché sur le trainer ; 6 clés mortes ; pas de hash de config dans les résultats. Un schéma qui n'est pas enforce est de la documentation, pas de la validation. | Indirect. | Un hash de config est exactement ce qu'un auditeur demande ; il manque. | — | — |
| **DataLoader reproductible** | Utile ; le padding dynamique existe. | Pas utilisé par le trainer (2,1× de calcul jeté) ; pas de `pin_memory`, `persistent_workers`, bucketing par longueur. | — | — | Gain immédiat 2× en temps d'époque. | — |
| **Multi-seeds + agrégation** | Bon réflexe méthodologique (3 seeds, moyenne ± σ). | 3 seeds est le minimum ; pas d'IC bootstrap, pas de test apparié, pas de correction pour comparaisons multiples (4 familles × plusieurs métriques au mois 8). | — | — | — | Le script réécrit `runs_summary.csv` (écrasement) sans versionner ; pas d'historique. |
| **Tests (689)** | Excellent signal de maturité pour un jury ou un recruteur. | Ne couvrent pas l'entraînement, le GPU, la reprise, ni la génération réelle. Pas de CI. | — | Un auditeur regarde les tests : il verra 689 verts et conclura à tort que le trainer est couvert. | — | Pas de tests de non-régression sur les artefacts. |
| **Métriques + baseline** | La PP et la baseline sont correctes *en soi* ; la comparaison est publiable après correction de la barre. | Métrique unique ; baseline surévaluée de 0,84 PP ; aucune métrique générative ; pas de décomposition par source/longueur. | Une fiche « ce que ces chiffres ne disent pas » serait nécessaire pour un lecteur biologiste : PP ≠ activité, PP ≠ qualité de génération. | — | — | La PP n'est pas monitorable en production (pas de données de référence à mesure que le corpus évolue). |
| **Tag Git `v0.3`** | Bonne discipline de versionnage. | Le tag pointe un état dont les artefacts de résultats sont absents ; le tag n'est donc pas auto-suffisant. | — | Un tag sans provenance ni hash de checkpoint n'a pas de valeur d'audit. | — | — |

### 2.2 Synthèse « prêt pour la production ? »

| Livrable | Prêt ? | Blockers (par ordre de gravité) |
|---|---|---|
| LSTM / Transformer | **Non** | Pas de checkpoint publié ni de model card ; pas de génération évaluée ; seed GPU cassée ; mémorisation non mesurée ; screening de biosécurité absent à la sortie. |
| Boucle d'entraînement | **Partiel** | Config non validée, resume cassé, pas de tracking d'expérience, déterminisme CUDA incomplet, zéro test. Utilisable en interne, pas en régime auditable. |
| Génération | **Non** | Tautologie de validité, pas de sortie structurée, pas de batch/cache (coût ×5 à ×10), seed GPU, pas de filtre toxine. |
| Configs | **Partiel** | Schéma non appliqué par le trainer ; clés mortes ; pas de hash dans les résultats. |
| DataLoader | **Partiel** | Padding dynamique non utilisé → 2× de coût ; pas de bucketing. |
| Multi-seeds | **Partiel** | 3 seeds, pas d'IC, pas de test apparié, pas de correction multiple. |
| Tests | **Partiel** | Aucun test du trainer, de la reprise, du GPU, de la génération réelle ; pas de CI. |
| Évaluation | **Non** | Une seule métrique, aucune mesure de diversité/noveté/fuite/propriétés sur les modèles, baseline biaisée, exemples non vérifiés. |
| **Ensemble « semaine 3 »** | **Non (pilote de recherche uniquement)** | La valeur industrielle actuelle est **nulle en tant que générateur** et **réelle en tant que harnais d'évaluation** — à condition de le terminer. Un CTO biotech n'investirait pas dans le générateur ; il pourrait financer le *benchmark* si la revue par les pairs est l'objectif. |

**Estimation de valeur, pour être concret.** Un criblage de 100 peptides synthétisés coûte ~5 000-15 000 € (synthèse + tests). Un modèle qui augmente le taux de succès de 2 % à 4 % double la valeur d'un écran, soit 5-15 k€ par campagne. Un LSTM de 0,9 M paramètres entraîné sur 38 k séquences publiques **n'atteint pas ce niveau de preuve** : rien dans le dépôt ne démontre qu'il produit des séquences non mémorisées, et 50 % des exemples publiés sont des copies du train. Le produit commercialisable ici n'est pas le générateur, c'est **la procédure d'évaluation reproductible et sans fuite** — un actif de type « GuacaMol/MOSES pour peptides », qui a une valeur académique et open-source réelle, et un intérêt industriel comme standard commun.

---

## 3. ANALYSE CRITIQUE GLOBALE DU PROJET (vision stratégique)

### 3.1 SWOT

| **Forces** | **Faiblesses** |
|---|---|
| • Splits par familles avec `F(0,90) = 0` **vérifié** (0 cluster partagé train/val et train/test) — c'est plus rigoureux que la majorité des papiers de génération de peptides, qui splittent au hasard.<br>• Mesure explicite de l'effet de fuite : trigramme 19,53 (familles) vs 16,43 (aléatoire) — **3,1 PP d'optimisme évités**, chiffre publiable.<br>• Test set réellement verrouillé (aucun script ne le lit).<br>• 689 tests, tag Git, configs YAML, seed partout : discipline rare à ce stade.<br>• Contrainte « from scratch » assumée : le benchmark sera comparable, pas dépendant de poids propriétaires.<br>• Un vrai plan 12 mois avec jalons par famille — la seule façon d'éviter le « tout à la fin ». | • **73,6 % du corpus vient d'un seul dataset (AMPChallenge)** — le déséquilibre pharma/non-pharma annoncé n'existe pas (BACTIBASE/BIOPEP ≈ 0 %).<br>• **Biosécurité : 1 493 séquences flaguées toxines restent dans le train, 79 % du corpus n'est pas criblé, aucun filtre à la génération.**<br>• Évaluation réduite à la PP ; aucune preuve de génération réelle ; 2/4 exemples publiés copiés du train.<br>• Chaîne de provenance rompue (pas de runs.csv, pas de checkpoint, pas d'entrée PROVENANCE).<br>• Budget « égal » jamais défini : les deux références diffèrent déjà de 16 % en paramètres.<br>• `F(0,90) = 0` mais `F(0,70) = 0,607` : la fuite à seuil bas, celle qui gouverne la mémorisation, n'est pas traitée. |
| **Opportunités** | **Menaces** |
| • Le benchmark « 4 familles à budget égal, splits sans fuite, proxies calibrés, échecs rapportés » **n'existe pas** pour les peptides ; c'est une contribution publiable (NeurIPS D&B, *Bioinformatics*) si les métriques de génération sont au niveau.<br>• Les proxies physico-chimiques + PWM par position déjà codés (semaine 1) permettent une évaluation *distribution matching* sans 3D — différenciant vs les benchmarks qui exigent AlphaFold.<br>• L'inférence de composition globale (hypothèse §1.5.3) est un résultat scientifique original et peu coûteux à démontrer.<br>• Ouverture possible : réplication de résultats publiés (ProGen2, PepGLAD) sur splits sans fuite — beaucoup de résultats SOTA sont probablement gonflés par la fuite. | • **Un reviewer demandera la comparaison à ProGen/ESM-2** : sans poids pré-entraînés autorisés, la réponse doit être « comparaison contrôlée, pas SOTA », assumée dès l'introduction, sinon rejet.<br>• **Dual-use** : un générateur d'AMPs cationiques entraîné sur un corpus contenant des conotoxines, sans filtre, est un objet sensible ; la publication doit inclure un screening et une politique de non-diffusion.<br>• Risque de « negative result » : si les 4 familles atteignent toutes 15 ± 0,3 en PP et que la seule métrique est la PP, le mois 8 conclura « pas de différence » — un an pour un résultat nul.<br>• Risque de calendrier : l'évaluation n'est prévue qu'au mois 8 ; la construire après coup introduit un biais d'ajustement rétroactif.<br>• Licence des données : APD6 à redistribution restreinte (déjà gitignoré) — l'usage pour entraîner doit être documenté, sinon problème de conformité à la publication. |

### 3.2 Positionnement scientifique (honnête)

- **Ce n'est pas du SOTA et ne peut pas l'être.** Un LM de 0,9 M paramètres sur 38 k séquences se situe au niveau 2018-2019. ProGen (Madani et al., 2023) utilise 1,2 B paramètres et ~280 M de séquences ; ESM-2 va jusqu'à 15 B paramètres ; RFdiffusion/ProteinMPNN exigent la structure 3D (hors périmètre). La comparaison directe est donc impossible et doit être refusée explicitement.
- **Ce qui est défendable** : un **protocole** de comparaison contrôlée (même corpus, mêmes splits, même budget, mêmes seeds, mêmes métriques) sur quatre familles, avec une mesure quantitative du biais de fuite. C'est plus proche de GuacaMol/MOSES que de ProGen : la contribution est le *benchmark*, pas le modèle.
- **Attention à la revue de littérature du brief** : `Riffusion` est un modèle de diffusion audio (spectrogrammes) — le citer comme SOTA de génération de peptides est une erreur de catégorie qui ferait rejeter un papier. `DiffDock` est du docking (hors périmètre, pas de 3D). `ProtDiff` cible le backbone 3D. Les références pertinentes pour un benchmark de séquences 1D sont : ProGen/ProGen2, ProtGPT2, ESM-2 (comme plafond non atteignable), PepGLAD/PepFlow/DiffPepBuilder (diffusion pour peptides), AMP-Diffusion (2024) pour la génération d'AMPs. Toute citation doit être vérifiée avant le manuscrit (le plan annuel indique lui-même que ses références sont « de mémoire »).

### 3.3 Les cinq questions qui décideront du mois 8

1. **Quelle unité de budget ?** Paramètres, FLOPs, tokens vus, ou heures GPU ? Les quatre familles n'ont pas le même coût par pas. Sans définition écrite **avant** le mois 4, la comparaison sera contestée. Recommandation : FLOPs d'entraînement mesurés (formule analytique + `torch.utils.flop_counter`) à ±5 %, plus le temps mural en annexe.
2. **Quelles métriques, gelées quand ?** Il faut figer la liste et le code d'évaluation **avant** de regarder les résultats des familles 3 et 4, sinon on choisit la métrique qui favorise la dernière famille implémentée. C'est le même raisonnement que le test set verrouillé, appliqué aux métriques.
3. **Quel contrôle de nouveauté ?** Seuils 0,5/0,7/0,9 vs train, avec la mesure déjà disponible (`F_70 = 0,607` sur la val) comme référence : un générateur qui produit 60 % de séquences ≥ 70 % identiques à un du train ne génère rien, il paraphrase.
4. **Quelle preuve que le corpus n'est pas un dataset unique ?** 73,6 % AMPChallenge : il faut soit rééquilibrer, soit rapporter toutes les métriques **par source**, sinon le benchmark mesure la reproduction d'une campagne de design de 2018.
5. **Quel sort pour la biosécurité ?** Exclure les 1 493 flaguées du corpus génératif, cribler les 79 % « unknown » (au moins par motif), et **filtrer les sorties** avant toute diffusion.

### 3.4 Score global : **5,0 / 10**

Justification par composante (pondération « qualité scientifique ») :

| Composante | Note | Commentaire |
|---|---|---|
| Protocole de données / anti-fuite | 8/10 | Le meilleur du projet ; au niveau d'un papier de benchmark. |
| Ingénierie logicielle | 7/10 | Lisible, testé, mais le trainer est le trou. |
| Rigueur expérimentale (training) | 4/10 | Hyperparamètres non explorés, budget non défini, chiffres non reproductibles. |
| Évaluation générative | 2/10 | Quasi inexistante ; c'est la raison principale de la note globale. |
| Positionnement / originalité | 6/10 | Le créneau existe ; l'exécution ne le matérialise pas encore. |
| Traçabilité / conformité | 3/10 | Provenance absente là où elle compte. |
| **Global** | **5,0/10** | Un projet bien parti, avec une fondation « données » solide et un trou béant au milieu : **personne ne peut dire aujourd'hui ce que les modèles génèrent**. |

À titre de comparaison : un projet de benchmark publiable à NeurIPS D&B demanderait ≥ 8/10 sur évaluation, provenance et budget ; un projet industriel pilote ≥ 7/10 partout. On est à mi-chemin sur l'ingénierie, très en retard sur la mesure — ce qui est normal en semaine 3, mais **doit être corrigé avant le mois 4**, pas au mois 8.

---

## 4. PLAN D'AMÉLIORATION PRIORISÉ (mois 4-12)

Convention : Impact / Effort / Risque sur 1-5 (5 = max). Priorité **P0** = bloquant avant d'ouvrir D3PM (mois 4) ; **P1** = nécessaire avant le benchmark (mois 8) ; **P2** = industrialisation / manuscrit (mois 10-12). L'ordre des familles du plan annuel est respecté : aucune ligne ne propose d'ouvrir une famille avant son mois.

| # | Amélioration | Catégorie | Impact | Effort | Risque | Prio | Mois | Dépendances |
|---|---|---|---|---|---|---|---|---|
| 1 | **Évaluer la génération des modèles entraînés** : `scripts/generate_week3.py` (validité brute, taux de rejet par cause, proximité 0,5/0,7/0,9, unicité, Distinct-n, JS/KS par source) | Évaluation | 5 | 2 | 1 | **P0** | 4 | checkpoints (à produire) |
| 2 | **Réconcilier les chiffres** : `runs.csv` complet committé, hashes SHA-256 des checkpoints, section PROVENANCE semaine 3, correction de la figure de perte | Traçabilité | 5 | 1 | 1 | **P0** | 4 | — |
| 3 | **Corriger la seed de génération GPU** (`torch.Generator(device=…)`, passer dans tous les cas) + test paramétré par device | Génération | 4 | 1 | 1 | **P0** | 4 | — |
| 4 | **Corriger l'init weight-tying** (ne pas réinitialiser la tête ; asserter PAD = 0, std ≈ 0,02) + test | Architecture | 3 | 1 | 1 | **P0** | 4 | — |
| 5 | **Brancher `load_config` dans `train.py`**, supprimer/implémenter les clés mortes (`tie_weights`, `optimizer`, `early_stopping.metric`, `train.seed`) | Entraînement | 3 | 1 | 1 | **P0** | 4 | — |
| 6 | **Padding dynamique dans le trainer** (brancher `make_loader_s3` + `pad_mask`) → 2× plus rapide | Entraînement | 3 | 2 | 1 | **P0** | 4 | — |
| 7 | **Baselines honnêtes dans le tableau canonique** : unigramme positionnel (18,65), interpolation n-gram (19,19), ordre 3 sur split familles (19,53) ; retirer la baseline « random à longueur oracle » | Évaluation | 4 | 1 | 1 | **P0** | 4 | — |
| 8 | **Biosécurité à la génération** : exclure les 1 493 flaguées du corpus génératif, cribler les « unknown » par motif, filtre de sortie + journal | Éthique / conformité | 4 | 2 | 2 | **P0** | 4 | politique écrite |
| 9 | **Définir le budget égal** (FLOPs d'entraînement mesurés, ±5 %) et publier le tableau des 4 familles avec le même budget | Benchmark | 5 | 2 | 1 | **P1** | 4-5 | #6 |
| 10 | **Décodage batché + KV-cache** (×50-100 sur la génération) | Génération | 4 | 3 | 1 | **P1** | 4-5 | #3 |
| 11 | **D3PM (diffusion discrète masquée)** + schedule de bruit + loss CE, avec la même évaluation que #1 | Modélisation | 5 | 5 | 3 | **P0** | 4-5 | #1, #9 |
| 12 | **Scheduler cosinus + param groups** (wd=0 sur biais/LN) + 5 seeds + IC bootstrap + test apparié | Entraînement | 3 | 1 | 1 | **P1** | 5 | — |
| 13 | **GDSBM / Sinkhorn-IPFP** (mois 5) avec couplage documenté et coût mesuré | Modélisation | 4 | 5 | 3 | **P1** | 5-6 | #11 |
| 14 | **Model card + data card + DVC/MLflow + Docker + CI GitHub Actions** | MLOps | 4 | 3 | 1 | **P1** | 5-6 | #2 |
| 15 | **Schrödinger Bridge discret** (mois 6) sur le même harnais d'éval | Modélisation | 4 | 5 | 3 | **P1** | 6-7 | #13 |
| 16 | **Flow matching + VAE from scratch** (mois 7), même budget | Modélisation | 4 | 4 | 2 | **P1** | 7 | #15 |
| 17 | **Conditionnement** (préfixe / attributs de propriétés) + CFG, et *seam* d'API dès maintenant (`prefix`, `attributes`, `return_logprobs`) | Génération | 5 | 4 | 2 | **P1** | 7-9 | #10 |
| 18 | **Tests OOD** : longueurs 51-80, familles nouvelles (clusters réservés), propriétés extrêmes (charge > +10, GRAVY < −1) | Généralisation | 4 | 2 | 2 | **P1** | 8 | #17 |
| 19 | **Audit de fuite générative à seuils multiples** (50/60/70/80/90 %) + comparaison au taux de base `F_70 = 0,607` | Évaluation | 4 | 2 | 1 | **P1** | 8 | #1 |
| 20 | **Tests de significativité** : 4 familles × k métriques, correction Holm-Bonferroni, rapport des échecs | Statistiques | 3 | 2 | 1 | **P1** | 8 | #12 |
| 21 | **Test set : gouvernance** (hash public, accès journalisé, une seule évaluation finale, rapport de tout ce qui a été essayé) | Benchmark | 4 | 1 | 1 | **P1** | 8 | — |
| 22 | **Quantification / export** (ONNX, INT8) + serveur d'inférence + CLI FASTA/CSV avec métadonnées par séquence | Industrialisation | 2 | 4 | 2 | P2 | 10-11 | #10 |
| 23 | **Monitoring dérive** : suivi de la distribution des propriétés générées vs corpus, alerte sur dérive ; réentraînement documenté | Maintenance | 3 | 4 | 2 | P2 | 11-12 | #14 |
| 24 | **Politique dual-use publiée** (screening, non-diffusion des toxines, divulgation responsable) + section éthique du manuscrit | Éthique | 4 | 2 | 2 | **P1** | 4 puis continu | #8 |

**Chemin critique des 2 semaines qui viennent (mois 4, avant D3PM)** : #2 → #1 → #3/#4/#5/#6 (parallélisables) → #7 → #9 → #11. Les lignes #1 à #8 représentent **≈ 4-5 jours-homme** au total et débloquent tout le reste.

---

## 5. VERDICT FINAL

**Go / No-Go pour le mois 4 : GO conditionnel.**
Démarrer D3PM le mois 4 est légitime (le plan annuel l'exige, et le harnais d'entraînement existe), **à condition** que les 8 lignes P0 hors modélisation (#1 à #8 du §4) soient fermées **avant** la première comparaison D3PM vs LSTM. Aucune ne demande de recherche : ce sont des correctifs de plomberie et des mesures. Si elles ne sont pas faites, D3PM sera entraîné dans le même angle mort — et le mois 8 comparera quatre familles avec une métrique qui ne mesure pas ce qui compte, sur des artefacts invérifiables.

**Risque majeur (un seul, celui qui peut tuer le projet) :**
> **Le projet peut produire un an de travail non publiable parce que la chaîne de mesure et de preuve est en retard sur la chaîne de modélisation.** Concrètement : (a) la courbe de validation publiée n'est pas celle du run qu'elle illustre ; (b) 50 % des exemples générés publiés sont des copies du train, sans contrôle de nouveauté ; (c) 60,7 % de la validation a un homologue ≥ 70 % d'identité dans le train alors que le projet revendique « zéro fuite » ; (d) 73,6 % du corpus vient d'un seul dataset. Chacun de ces quatre points, découvert par un reviewer plutôt que par l'équipe, invalide la conclusion principale. Le risque n'est pas technique, il est **épistémique** : l'équipe croit mesurer une chose et en mesure une autre.

**Recommandation #1 (à faire lundi matin, avant toute ligne de D3PM) :**
> **Écrire `scripts/generate_week3.py` et publier un tableau de génération pour les modèles réellement entraînés** — validité brute (avant filtrage), taux de rejet par cause, identité maximale vs train à 0,5/0,7/0,9, unicité à n fixé, Distinct-2/3, JS des propriétés (charge, GRAVY, moment hydrophobe, pI) et KS des longueurs, **par source**. Ce script unique répond simultanément aux points (a), (b) et (c) : il rend les chiffres vérifiables, il dit si les exemples publiés sont de la mémorisation, et il produit la référence à laquelle les trois familles suivantes seront comparées. C'est aussi le seul livrable de la semaine 3 qui manque pour que le mot « génératif » ait un sens.

**Recommandation #2 (gratuite, à faire dans la même journée)** : reformuler la revendication de reproductibilité en « même machine, même build, même nombre de threads », et ajouter `torch.backends.cudnn.benchmark = False` + `CUBLAS_WORKSPACE_CONFIG=:4096:8`. Une revendication exacte vaut mieux qu'une revendication absolue et fausse.

**Recommandation #3 (pour le professeur)** : la partie « données » de ce projet est au niveau d'un papier de benchmark ; la partie « évaluation de la génération » est au niveau d'un TP. Le déséquilibre est dangereux pour la note finale : un jury lira d'abord le tableau de résultats de génération, et il n'existe pas. Rééquilibrer les deux prochaines semaines vers la mesure, pas vers un troisième modèle.

---

## Annexe A — Preuves exécutées pour ce rapport

```text
# 1. Suite de tests
python -m pytest -q                     -> 689 passed in 72.36s

# 2. Couverture du trainer
grep -r "train_model|_run_epoch_train" tests/   -> aucun résultat

# 3. Bug d'initialisation liée (transformer.py)
CausalTransformer(...)  -> token_emb.weight.std() = 0.1150   (attendu 0.02)
                        -> token_emb.weight[PAD].abs().sum() = 12.79  (attendu 0)
                        -> head.weight is token_emb.weight = True
                        -> forward batch paddé : 0 NaN, 0 gradient NaN (masque OK)

# 4. Paramètres et latence (CPU, batch = 1)
LSTM        930,455 params   20 séquences en 6.15 s   -> 307 ms/séquence
Transformer 802,967 params   20 séquences en 31.27 s  -> 1564 ms/séquence

# 5. Baselines recalculées sur data/splits/{train,val}_v02.csv
ordre 0 (uniforme 21)                  PP = 21.00
unigramme add-1 (publié)               PP = 19.49
bigramme add-1                         PP = 19.27
trigramme add-1                        PP = 19.53
interpolation Jelinek-Mercer ordre 2   PP = 19.19
unigramme positionnel (sans modèle)    PP = 18.65
(position, AA précédent) add-1         PP = 18.74

# 6. Audit du corpus
train 38,336 / val 2,095 / test 2,144 = 42,575
AMPChallenge 73.6 % | UniProt 17.5 % | APD3 1.4 % | BACTIBASE/BIOPEP 0 %
toxines flaguées : train 1,493 (3.9 %) | val 161 (7.7 %) | test 132 (6.2 %) | unknown 79 %
clusters partagés train/val = 0 ; train/test = 0 ; F_70 (split familles) = 0.607

# 7. Audit des exemples publiés (visuals/3_exemples_generes.png)
KLLKALKAALLKA      -> absent du corpus
GLFDIVKKVVGALG     -> présent à l'identique dans train_v02.csv
KWKLFKKIEKVGQN     -> préfixe exact (14/38 AA) d'une séquence de train_v02.csv
IIGKLKSLFKKAAK     -> absent du corpus

# 8. Ré-entraînement LSTM seed 42 (config publiée, CPU, interrompu après l'époque 6)
époque 1 : train_loss = 2.7884 | val_loss = 2.8309 | val_pp = 16.96
époque 2 : train_loss = 2.6623 | val_loss = 2.7919 | val_pp = 16.31
époque 3 : train_loss = 2.5396 | val_loss = 2.7433 | val_pp = 15.54
époque 4 : train_loss = 2.4234 | val_loss = 2.7338 | val_pp = 15.39
époque 5 : train_loss = 2.3309 | val_loss = 2.7441 | val_pp = 15.55
époque 6 : train_loss = 2.2562 | val_loss = 2.7332 | val_pp = 15.38
figure publiée, époque 1 : train_loss = 2.78 | val_loss = 2.69
-> la courbe d'entraînement publiée est reproduite (2,78 vs 2,7884) ;
   la courbe de validation publiée ne l'est pas (2,69 vs 2,8309), et sa valeur
   finale (2,06 -> PP 7,9) est incompatible avec le 15,21 du tableau.
   Le tableau est plausible ; la figure de validation ne provient pas de ce run.

# 9. Justification technique du choix norm_first=False
inspect(nn.TransformerEncoder.__init__) -> le fast path nested tensor est bien
   désactivé si norm_first=True (commentaire du code exact) ; mais
inspect(nn.TransformerEncoder.forward) -> ce fast path exige mask is None et
   est désactivé en mode entraînement : il est donc inutilisé ici (masque causal
   dense 52x52 passé à chaque appel). Le compromis invoqué est sans objet.
```

## Annexe B — Références utilisées

- **Métriques de benchmark génératif** : Brown et al., *GuacaMol* (2019) ; Polykovskiy et al., *MOSES* (2020) — validité/unicité/novelty/Self-BLEU/Distinct-n, et la règle « validité brute avant filtrage ».
- **Weight tying et initialisation** : Press & Wolf (2017) ; Inan et al. (2017) — embeddings liés de faible amplitude.
- **Reproductibilité PyTorch** : documentation officielle « Reproducibility » (`torch.use_deterministic_algorithms`, `CUBLAS_WORKSPACE_CONFIG`, `cudnn.benchmark`).
- **Contexte biologique / SOTA** : Madani et al., *ProGen* (2023) ; Nijkamp et al., *ProGen2* (2023) ; Lin et al., *ESM-2* (2023) ; Ferruz et al., *ProtGPT2* (2022) ; Watson et al., *RFdiffusion* (2023) ; Dauparas et al., *ProteinMPNN* (2022) ; PepGLAD / PepFlow / DiffPepBuilder (2024) ; AMP-Diffusion (2024).
- **Gouvernance** : Mitchell et al., *Model Cards* (2019) ; Gebru et al., *Datasheets for Datasets* (2018) ; FDA, *AI/ML in Drug Development* (2023) ; EMA, *Reflection Paper on the Use of AI* (2024).
- **Statistique** : Efron & Tibshirani, bootstrap percentile ; Holm (1979) pour la correction des comparaisons multiples.

> **Limites de ce rapport** : je n'ai pas pu ré-entraîner les 6 runs complets ni générer depuis leurs checkpoints (non publiés). Les conclusions sur la PP finale reposent sur la cohérence interne des artefacts publiés, sur mon ré-entraînement partiel et sur les mesures de baselines. Si les checkpoints existent hors dépôt, les points §1.2.11 et §1.3.6 doivent être re-vérifiés avec eux — c'est précisément le correctif P0 n°2.
