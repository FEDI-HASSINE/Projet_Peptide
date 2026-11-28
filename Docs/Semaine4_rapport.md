# Semaine 4 · Données, mémorisation et provenance

Contribution de Yassine Achour au 28 novembre 2026. La partie mesurée couvre la représentation des données et le masquage. Le NELBO, la comparaison AR/diffusion et les mesures expérimentales de mémorisation sont non disponibles.

## Données et contrôles mesurés

Les SHA-256 et les effectifs de train (38 336) et validation (2 095) correspondent à `splits_manifest.json`. Le contrôle ouvre uniquement ces deux CSV. Les 38 336 lignes du train ont chacune passé l’aller-retour dans les deux vocabulaires. Le tokenizer S3 place PAD/BOS/EOS aux indices 0/1/2. Dans S4, les résidus occupent 0 à 19, PAD/BOS/EOS occupent 20/21/22 et MASK vaut 23. Importer les constantes du module correspondant à la famille.

`DiffusionDataset("train")` et `DiffusionDataset("val")` fournissent des lignes int64 de 51 positions, avec EOS après le dernier résidu, puis PAD. Le chargeur existant `peplib.datasets.make_loader` les assemble. `mask_batch(x0, t, generator)` renvoie le lot, un masque booléen et un entier commun au lot. Si le temps est omis, il est tiré uniformément de 1 à 50. EOS est éligible. PAD et BOS sont protégés. L’ajout de padding conserve les tirages sur les positions valides à graine égale.

| Étape | Probabilité attendue | Masquées sur 100 000 | Fréquence mesurée |
|---:|---:|---:|---:|
| 1 | 0,02 | 1 960 | 0,01960 |
| 10 | 0,20 | 20 056 | 0,20056 |
| 25 | 0,50 | 50 081 | 0,50081 |
| 50 | 1,00 | 100 000 | 1,00000 |

Le [tableau complet](../results/week4/masking_control.csv) couvre 50 étapes et cinq millions de positions synthétiques. À la graine 42, les 50 effectifs entrent dans leur intervalle prédictif binomial central, avec une couverture simultanée d’au moins 95 % par correction de Bonferroni. Sur 10 000 lots supplémentaires, aucun PAD ni BOS masqué. Ces mesures caractérisent le chargeur. La qualité du modèle reste à mesurer.

## Échantillonnage et audit

**T79.** `count_block_passes(model, K)` entoure le traitement d’un lot par T79 et compte les entrées dans `model(...)`. Consigner le K demandé, la taille du lot et les appels observés. Le test sur un module identité vérifie que trois appels restent comptés comme trois, même si K vaut 50. La méthode par blocs est une approximation. L’interface T78, les mesures de qualité et la comparaison statistique à K = T sont non disponibles.

**T83.** `memorization_report(records, train, flags)` accepte les colonnes `sequence`, `toxin_flag` et le registre v02 `sequence`, `toxicity_status`. Seuls les enregistrements explicitement négatifs et annotés dans le registre entrent dans la distribution. Un signalement portant sur un doublon exclut toutes ses occurrences. Les agrégats comprennent effectifs admissibles/exclus, histogramme, proportion à 0,90 et IC de Wilson à 95 %. Si aucun enregistrement n’est admissible, la proportion et son intervalle sont non disponibles. Utiliser exclusivement ces agrégats pour les CSV et figures. L’interface ne fournit aucun exemple individuel.

La mesure T83 applique l’identité de Levenshtein aux fragments sans gaps du meilleur alignement local de Smith–Waterman, puis prend le maximum sur le train. Un court fragment identique peut donc obtenir 1. La métrique S3 `best_min` combine identité et couverture et utilise un préfiltrage. L’équipe doit valider la définition retenue pour T83 et comparer validation/train avec les mêmes exclusions et le même dénominateur. Le taux obtenu sur le seul sous-ensemble annoté ne peut pas être étendu à toutes les générations. Un résultat négatif au filtre de toxines ne démontre pas l’innocuité. Les descripteurs physico-chimiques restent des proxys. Toute mesure d’activité exige une validation en laboratoire qualifié.

## Passage à l’équipe

- **Fedi · T69/T70.** Les contributions distantes antérieures de Yassine sont incluses dans main. Le tag disponible est `v0.3-semaine3` sur `7ee2040`. Le nom exact `v0.3` et le protocole S4 signé sont absents des références récupérées. Dans `results/week3/`, seules la table des baselines et quatre figures sont suivies. Réconcilier les trois CSV de clôture, la provenance et le seuil encore affiché à 19,49 dans le rapport S3.
- **Youssef · T73/T78/T81.** Confirmer les identifiants S4, le temps scalaire par lot et la concordance avec T73. Fournir l’interface T78 et les métadonnées des runs autorisés pour 42/43/44. Le dépôt inspecté ne fournit ni config S4 ni script d’entraînement diffusion. `verify_runs` vérifie champs, graines, unicité des epochs par graine et hash YAML canonique. `verify_checkpoints` compare les noms, tailles et SHA-256 du manifeste aux fichiers conservés dans le répertoire privé du run.
- **Sofiene et Yassine · T79/T83.** Comparer les statistiques à K = T sur les artefacts autorisés et publier séparément qualité et coût des blocs. Calculer la référence validation/train et fournir les agrégats de mémorisation pour le même périmètre d’annotations. Résultats, NELBO et équité du budget sont non disponibles.
- **Tous · T87/T88.** Relecture croisée, comparaison de deux exécutions CPU de référence et vérification sous Python 3.11 et 3.14 restent à faire. Le coordinateur conserve le gel, T86, le merge et le tag.

## Reproduire les contrôles de données

Depuis la racine, avec le Python du venv, `python -m scripts.check_week4_artifacts` contrôle les deux splits autorisés. `python -m scripts.check_masking_s4` imprime le CSV de contrôle. Le SHA-256 publié porte sur les octets UTF-8 avec fins de ligne LF.

Le contrôle ciblé couvre les trois fichiers de tests S4 ainsi que `tests/test_tokenizer.py` et `tests/test_biosafety.py`. Résultat final sous Python 3.13.14 et PyTorch 2.14.1+cpu : 67 réussis, zéro échec. La suite complète contient des tests qui référencent les splits verrouillés. Le fichier local non suivi est resté hors de la livraison.
