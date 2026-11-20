# NOTE DE PASSAGE -- Youssef -> Chef / Gestionnaire
# Semaine 3 -- Version v2 (corrections PR #25)
# Date : 20/11/2026

## Corrections apportees (PR #25 -> v2)

| Bloquant | Correction |
|----------|-----------|
| Splits naifs -> fuite | Supprime make_splits.py ; configs -> data/splits/*_v02.csv |
| runs.csv invalide | Reinitialise ; relancer les 6 runs sur v02 |
| Generation absente | peplib/generate.py + tests/test_generate.py |
| Configs conflictuelles | Schema imbrique (PR #24) adopte dans les deux configs |

## Etat actuel des runs (sur splits v02)

| Modele | Seed | Statut | Splits utilises |
|--------|------|--------|-----------------|
| LSTM | 42 | A LANCER | data/splits/train_v02.csv |
| LSTM | 43 | A LANCER | data/splits/train_v02.csv |
| LSTM | 44 | A LANCER | data/splits/train_v02.csv |
| Transformer | 42 | A LANCER | data/splits/train_v02.csv |
| Transformer | 43 | A LANCER | data/splits/train_v02.csv |
| Transformer | 44 | A LANCER | data/splits/train_v02.csv |

NOTE : Les runs precedents (val_pp=7.07 LSTM-42) etaient sur un split fuyard
       et ne sont pas valides. A recalculer.

**Baseline UNIGRAMME : val_pp = 19.49** (correction : c'est l'unigramme, pas le bigramme)

## Commandes pour relancer tous les runs

```bash
# Installer pyyaml si necessaire
.venv\Scripts\pip.exe install pyyaml

# Lancer les 6 runs (LSTM 42/43/44 puis Transformer 42/43/44)
.venv\Scripts\python scripts\run_all_seeds.py

# Ou separe :
.venv\Scripts\python scripts\run_all_seeds.py --only lstm
.venv\Scripts\python scripts\run_all_seeds.py --only transformer
```

## Critere de succes

- val_pp < 19.49 (baseline unigramme) pour au moins LSTM seed 42 sur val_v02
- results/week3/runs_summary.csv : tableau moyenne +/- ecart-type
- 17 tests pytest verts (test_lstm, test_transformer, test_generate)

## Fichiers livres dans cette PR

| Fichier | Role |
|---------|------|
| peplib/models/lstm.py | T52 -- LSTM causal |
| peplib/models/transformer.py | T55 -- Transformer causal |
| peplib/train.py | T53 -- boucle (AdamW, splits v02, generateur seede) |
| peplib/generate.py | NOUVEAU -- generation BOS->EOS, sampling.py integre |
| configs/lstm_base.yaml | Config imbrique, splits v02 |
| configs/transformer_base.yaml | Config imbrique, splits v02 |
| scripts/run_all_seeds.py | Lance 6 runs, agregation moy+/-std |
| tests/test_lstm.py | 8 tests LSTM |
| tests/test_transformer.py | 9 tests Transformer |
| tests/test_generate.py | NOUVEAU -- 5 tests generation |
