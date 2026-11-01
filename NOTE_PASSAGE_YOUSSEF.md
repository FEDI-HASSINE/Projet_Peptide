# NOTE DE PASSAGE — Youssef → Chef / Gestionnaire
# Semaine 3 — Entraînements restants (01/11/2026, 13h15)

## État actuel des runs

| Modèle | Seed | Statut | Best val_pp | Epochs |
|--------|------|--------|-------------|--------|
| LSTM | 42 | **TERMINE** | **7.07** | 30/30 |
| LSTM | 43 | INTERROMPU epoch 15 | 7.48 (en cours) | 15/30 |
| LSTM | 44 | A LANCER | — | 0 |
| Transformer | 42 | A LANCER | — | 0 |
| Transformer | 43 | A LANCER | — | 0 |
| Transformer | 44 | A LANCER | — | 0 |

**Baseline bigramme : val_pp = 19.49**
**LSTM 42 : val_pp = 7.07 → -64% ✅**

---

## Pré-requis

```bash
# Installer pyyaml si pas déjà fait
.venv\Scripts\pip.exe install pyyaml

# Vérifier que les splits existent (déjà créés et commités)
ls data/processed/splits/
```

---

## Commandes à lancer dans l'ordre

### 1. LSTM seed 43 (reprend depuis epoch 15 si best.pt présent, sinon repart de 0)
```bash
.venv\Scripts\python.exe -m peplib.train --config configs/lstm_base.yaml --seed 43 --resume
```
Durée estimée : ~15 epochs × 2 min = ~30 min (si résumé depuis epoch 15)

### 2. LSTM seed 44
```bash
.venv\Scripts\python.exe -m peplib.train --config configs/lstm_base.yaml --seed 44
```
Durée estimée : ~30 epochs × 2 min = ~60 min

### 3. Transformer seed 42
```bash
.venv\Scripts\python.exe -m peplib.train --config configs/transformer_base.yaml --seed 42
```
Durée estimée : ~30 epochs × 2-3 min = ~60-90 min

### 4. Transformer seed 43
```bash
.venv\Scripts\python.exe -m peplib.train --config configs/transformer_base.yaml --seed 43
```

### 5. Transformer seed 44
```bash
.venv\Scripts\python.exe -m peplib.train --config configs/transformer_base.yaml --seed 44
```

### OU — Tout lancer d'un coup (séquentiel)
```bash
# LSTM 43 (resume) + 44 + Transformer 42/43/44 d'un seul coup
.venv\Scripts\python.exe scripts\run_all_seeds.py --only lstm --resume  # fait 43(resume)+44
.venv\Scripts\python.exe scripts\run_all_seeds.py --only transformer     # fait 42+43+44
```

---

## Résultats en temps réel

Les runs s'enregistrent automatiquement dans `results/week3/runs.csv` (append).

```bash
Get-Content results\week3\runs.csv
```

---

## Critère de succès

- val_pp < 19.49 pour tous les runs (baseline bigramme) ✅ déjà atteint sur LSTM 42
- Checkpoints `best.pt` dans `checkpoints/{lstm,transformer}/seed_{42,43,44}/`
- `results/week3/runs.csv` avec les 6 runs complets (ou autant que possible avant lundi 12h)

---

## Où sont les fichiers

| Fichier | Description |
|---------|-------------|
| `peplib/models/lstm.py` | T52 — LSTM causal |
| `peplib/models/transformer.py` | T55 — Transformer causal |
| `peplib/train.py` | T53 — Boucle d'entrainement |
| `configs/lstm_base.yaml` | Config LSTM |
| `configs/transformer_base.yaml` | Config Transformer |
| `results/week3/runs.csv` | Log de tous les epochs |
| `checkpoints/lstm/seed_42/best.pt` | Meilleur LSTM 42 (local, non pushé) |
| `checkpoints/lstm/seed_43/best.pt` | Meilleur LSTM 43 epoch 15 (local, non pushé) |
| `tests/test_lstm.py` | 8 tests LSTM (tous ✅) |
| `tests/test_transformer.py` | 9 tests Transformer (tous ✅) |
