# CV-level comparison: our engine vs OUTRIDER

Same 5-fold StratifiedKFold(by Batch_ID, seed=42) split on the 676 HC samples used in both. Reported on the **same 12,305-gene subset** (OUTRIDER's FPKM>=1-filtered set, which is a strict subset of our engine's 19,858) unless noted -- medians are the primary comparison since a handful of extreme-outlier genes dominate the means on both sides.

## 1. Z-score moments (held-out HC, pooled across 5 folds)

| metric | our engine (raw, pre-SHASH) | our engine (SHASH-corrected) | OUTRIDER (raw) |
|---|---|---|---|
| mean(Z) | 0.0017 | -- | -0.00001 |
| std(Z) | 1.039 | -- | 0.999 |
| skew(Z) | **-0.230** | 0.095 (full-universe number; see caveat) | **-2.03** |
| kurtosis(Z) | **0.560** | 0.417 (full-universe number) | **6.85** |
| naive_exceed \|Z\|>1.96 (nominal 0.05) | 0.058 | 0.053 | 0.050 |

**Key finding**: mean/std/exceedance-rate are all comparably well-calibrated between the two. But **shape** is not -- OUTRIDER's raw Z is strongly left-skewed with heavy tails (skew -2.0, kurtosis 6.9, vs our -0.23/0.56 even before our own SHASH correction). OUTRIDER has no equivalent second-stage shape correction (our SHASH step exists precisely to fix this kind of residual non-normality before FDR control is trusted) -- so a naive |Z|>1.96 or normal-theory p-value on OUTRIDER's output inherits that skew directly into any downstream FDR calculation, right where our engine has already spent a dedicated calibration stage on exactly this problem.

## 2. Posterior predictive check (obs vs. model-implied NB moments)

| metric | our engine (12,305-gene subset) | OUTRIDER |
|---|---|---|
| obs_zero_frac | 0.059 | 0.059 (same data, same genes -- sanity check, must match) |
| pred_zero_frac | 0.037 | 0.023 |
| zero_diff | -0.018 | -0.024 |
| pearson_chi2 (mean 1.0 = well-specified) | 0.961 | 0.875 |
| mean_rel_err (median) | 0.015 | 0.017 |
| var_rel_err (median) | -0.062 | -0.440 |

Median-level fit quality is broadly comparable, OUTRIDER slightly under-predicts variance more (-0.44 vs -0.06 median relative error) but isn't dramatically worse on the genes it can fit.

**Outlier-gene robustness (mean, not median -- shows tail behavior)**:

| metric (mean, not median) | our engine | OUTRIDER |
|---|---|---|
| mean_rel_err | ~500 (few pathological genes) | ~5.3e7 |
| var_rel_err | ~1.2e6 | ~7.4e21 |

Both sides have a handful of pathological outlier genes that blow up the mean-based statistic (this is expected with per-gene NB fits at this scale) -- but OUTRIDER's worst-case blowups are ~15 orders of magnitude larger. Our engine clips fitted mu to [1e-6, 1e8] (`model_engine_mixed.py:score()`) and applies EB dispersion shrinkage specifically to prevent single-gene MLE degeneracy; OUTRIDER's `normalizationFactors`/`theta` have no equivalent safety bound in this pipeline.

## 3. Gene coverage, now including BOTH engines' own modeling failures

| | universe | production/full-fit failures | CV-fold-level failures (any of 5 folds) | net coverable |
|---|---|---|---|---|
| Our engine | 20,097 protein-coding | 239 (route=excluded, both cascade stages failed) | 455/88,990 gene-fold rows (0.51%) -- affects some genes in only some folds | 19,858 in production; ~19,600-19,858 robustly stable across all 5 CV folds (need exact per-gene fold-consistency count, see caveat) |
| OUTRIDER | 19,858 (already our engine's own successful set) | 7,194 fail FPKM>=1 filter (upstream, deterministic, same genes every fold) | 7 genes fail NB optimizer, **identical 7 genes dropped in all 5 folds** | 12,305, stable across folds |

**This directly answers the fairness question raised earlier**: our engine is *not* 100% -- 239/20,097 (1.2%) fail outright in production, plus a further 0.51% of gene-fold cells fail only within specific CV folds (sample-subset-dependent instability, not a full exclusion). OUTRIDER's failure structure is different in kind: its 7 CV-failing genes are **exactly the same 7 in every fold** (deterministic optimizer divergence on those specific genes, not stochastic fold-sampling sensitivity), and its dominant loss (7,194 genes, 36.9%) happens entirely upstream of any fold-level fitting, via the FPKM filter alone. Net: our engine still covers 98.8% vs OUTRIDER's 62.0%, and the failure *mechanisms* differ (our residual 1.7% combined failure rate is spread thin and fold-sensitive; OUTRIDER's 38.0% is concentrated and deterministic, dominated by the sparse/low-expression tail this project's low-expression modeling strategy was built to rescue).

Rolled up to distinct genes (model-route only, 17,798 genes): 17,473 (98.2%) succeed in **all 5** CV folds; 325 (1.8%) fail in at least one fold (218 fail in 4/5 folds, 87 in 3/5, 17 in 2/5, 3 in 1/5) -- a graded, fold-sampling-sensitive failure mode, not a hard binary exclusion.
