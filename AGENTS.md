# AGENTS.md — Projet Peptide

Research/training course (ISIMS, supervisor A. Ben Abdallah): 12 months to build a leak-free pharma + non-pharma peptide benchmark and compare four generative-model families — every line coded, tested and documented by the students.

## State of the repo — read this first
- **Specs only, no `peplib/` code yet.** Git is initialized on branch `main` (Week-1 T1 done); directory tree, `.gitignore`, `requirements.txt` and `README.md` exist. Python files and tests remain to be written; tag `v0.1` is not yet set.
- All specs are in **French**, now under `Docs/`. `Docs/Semaine1_peptides_peplib_4etudiants.md` is the most complete Week-1 text. `Docs/plus.txt` is scratch notes, **not** a spec.
- Empty directories are kept in git with `.gitkeep` placeholders — git does not track empty folders.

## Toolchain — hard constraints
- **Allowed:** Python ≥ 3.10, NumPy, SciPy, pandas, matplotlib, `requests` (UniProt REST), pytest, Git, PyTorch (`nn.Module`, autograd, optimizers only, `torch>=2.0`), LaTeX/Beamer (months 10–12). scikit-learn **only** to cross-check the team's own metrics (month 3).
- **Forbidden:** GNNs (`torch_geometric`), pretrained protein models (ESM, ProtT5), structure predictors/designers (AlphaFold, ESMFold, Boltz, RFdiffusion), PyMOL and any 3D viewer, BLAST/MMseqs2/cd-hit/HMMER, turnkey generative libraries. No pretrained weights ever.
- **Limitation to declare in every writeup:** with no structure predictors, structural constraints stay proxies (helicity, amphipathicity, HP lattice) and scorer outputs are not activity measurements. No biological claim without qualified-lab validation.

## Planned layout (Week-1 spec, `peplib/` still empty)
```
README.md  requirements.txt  PROVENANCE.md
data/raw/      # immutable downloads, SHA-256 logged
data/processed/  # CSVs: sequence, source, id
peplib/constants.py   # 20 AAs: monoisotopic mass, side-chain pKa, Kyte-Doolittle, Eisenberg
peplib/filters.py     # alphabet + length filter, returns kept + rejection dict by cause
peplib/stats.py       # composition, length histogram, KL/JS, entropy
peplib/pwm.py         # PWM from N/C termini, per-position entropy/info
peplib/descriptors.py # net charge, pI, GRAVY, hydrophobic moment, Boman, mass
peplib/dedup.py       # exact dups, inverted k-mer index, Levenshtein, identity
tests/  notebooks/NB01_peptides.ipynb  models/  results/
```
`tests/` gets one file per module. **No business logic in notebook cells** — `NB01` imports and calls `peplib` only.

## Commands
```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1     # Windows / PowerShell
pip install -r requirements.txt   # pinned with `pip freeze`
pytest                          # every function must have a unit test
```
Seed **42** everywhere: `random`, `numpy`, `torch.manual_seed` (+ `manual_seed_all` if CUDA). CPU only expected in month 1 — no training; note `torch.cuda.is_available()` in the README.

## Constants the tests assert — transcribe exactly
- Alphabet regex `^[ACDEFGHIKLMNPQRSTVWY]+$`; length window **10–50 inclusive**.
- Vocab = 20 AA + PAD + BOS + EOS = **23**. Month-2 unigram baseline uses **21** tokens (20 AA + EOS).
- Rejection causes apply in this order, each sequence counted once in its **first** cause: missing → empty → ambiguous (`X B Z J`) → rare (`U O`) → other characters → < 10 → > 50. Invariant: kept + rejected = raw total.
- Citation required per constants table: side-chain pKa, **N-terminus 9.69 / C-terminus 2.34** (justify, e.g. EMBOSS/ExPASy), monoisotopic residual mass (A 71.03711, G 57.02146), Kyte–Doolittle (A 1.8, G −0.4), Eisenberg. One typo in a scale invalidates every downstream result — cross-review and test against known values.
- Hand-check anchors: `GRAVY("A")=1.8`, `GRAVY("GG")=-0.4`, `GRAVY("AG")=0.7`, `mass("G")≈75.0320`, `mass("AA")≈160.0848` (Σ residual masses + H₂O 18.01056), `Levenshtein("kitten","sitting")==3`.
- **Log bases are not uniform**: per-position entropy/information and PWM use **log₂** (`log₂ 20 ≈ 4.3219`); Jensen–Shannon is in **nats**, bounded by ln 2. Don't unify them.

## Implementation gotchas (the spec's own trap list)
- **`dropna()` before any `astype(str)`** — the reverse turns `NaN` into `"NAN"`, which passes the regex. Order: `dropna()` → `strip()` → `upper()` → filters.
- PWM: index from the N-terminus (position 0) and from the C-terminus (position −1). **Never truncate to the shortest sequence** — truncation is not alignment.
- Mask exact zeros (`0 log 0 = 0`) or write down a justified pseudo-count; never let `-inf` reach a figure.
- Per-position information needs the small-sample correction `e(n) = (s−1)/(2 ln 2 · n)`, `s = 20`; without it `I(i)` is biased upward at small `n` (Schneider & Stephens 1990).
- Net charge at pH 7 via Henderson–Hasselbalch **including the N and C termini**; pI by bisection (charge strictly decreasing in pH); hydrophobic moment with δ = 100° for an α-helix.
- Near-duplicates: inverted index over shared 3-mers → candidate pairs, then exact Levenshtein; identity `= 1 − d/max(|a|,|b|)`. Identical sequences *are* candidates; sequences sharing no 3-mer are *not*.
- No unjustified numerical epsilon; document and validate every approximation.

## Git & review workflow
- Branch off `main` per person: `feat/<name>-<role>` (e.g. `feat/yassine-donnees`).
- Every PR reviewed by **someone other than the author**; `pytest` green before merging.
- **Only the coordinator merges into `main`** (Week-1 "superchef"; the annual plan defines 3 rotating roles E1 data / E2 models / E3 evaluation, rotating again at the month-6 and month-9 milestones).
- `.gitignore` covers `__pycache__/`, `*.pyc`, `.venv/`; raw data zips either gitignored or LFS if > 50 MB. Tag per milestone (`v0.1` end of Week 1).

## Data provenance & biosafety
- `data/raw/` is immutable; record each file's SHA-256 and write a `PROVENANCE.md` entry per source (name, URL/query, date, version, license, raw line count, hash, manual-download flag).
- Everything in `data/processed/` must be regenerable from `data/raw/` by code, columns `sequence, source, id`.
- Sources: UniProt REST (`reviewed:true` first, cursor pagination via the `Link` header), one AMP database (APD3 recommended; fall back to DBAASP or DRAMP if the license blocks it), one non-pharma source (BIOPEP-UWM, or BACTIBASE for bacteriocins). Check usage/redistribution terms *before* downloading.
- UniProt holds toxins (conotoxins): flag or exclude them from the generative corpus and **never distribute toxin-like generated sequences**.

## Where the project is
Month 1 of 12 (corpus + `peplib` v0.1, no training). **The month order is fixed, the calendar is not** — never open a new model family before the previous milestone. Then: 2 NW/SW + BLOSUM62, union-find clusters at 90% identity, leak-free 90/5/5 splits by clusters, LSTM/Transformer baselines · 3 Chou–Fasman/GOR, HP lattice, calibrated scorers · 4 masked discrete diffusion · 5 D3PM + discrete Sinkhorn/IPFP · 6 GDSBM reproduction + discrete Schrödinger bridge · 7 continuous flow matching + from-scratch VAE · 8 Dirichlet/discrete flow matching, 4-family equal-budget benchmark · 9 conditional generation + Pareto · 10 robustness, uncertainty, ablations · 11 validation plan + manuscript · 12 reproducibility, defense, preprint.

Expect only ~10⁴–10⁵ sequences: keep the Transformer small, regularize, early-stop, report intervals over seeds. Gates to respect later: causal-mask test (future tokens must not change past logits), PAD/BOS logits set to −∞ when sampling, the 21-token unigram baseline must be beaten, every experiment relaunchable from a versioned config, test set locked before final comparisons, equal compute budget, and **report failures, not just the best result**.

## Source documents
- `Docs/Plan_annuel_peptides_12_mois.docx` — 12-month plan, governance, risks, references (citations there are from memory; verify before quoting).
- `Docs/Semaine1_peptides_peplib.docx` / `Docs/Semaine1_peptides_peplib_4etudiants.md` — Week-1 tasks T1–T29, exact formulas (§5), mandatory tests (§6), traps (§7), `PROVENANCE.md` template, quick-start commands. Read these instead of re-deriving formulas.
- `Docs/Les 10 meilleures pistes non pharmaceutiques pour les peptides.odt` — the non-pharma tracks (biocosmetics P1, biopreservation P2, …) that define the non-pharma scorer brief.

## Open
Comment and commit language (French vs English) is undecided — ask before writing prose.
