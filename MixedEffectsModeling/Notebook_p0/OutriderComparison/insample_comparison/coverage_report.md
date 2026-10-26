# Gene Coverage: Our Engine vs OUTRIDER (HC cohort, n=676)

## 1. Overall coverage

| | genes attempted | genes modeled | coverage |
|---|---|---|---|
| Our engine (model + pool route) | 20,097 protein-coding | 19,858 | 98.8% |
| OUTRIDER (official `filterExpression`, FPKM>=1 in >=95th percentile sample) | 19,858 | ~12,530 pass filter, 12,305 in final fit (7 dropped: NB optimizer divergence) | 62.0% |

OUTRIDER's own hard minimum (`checkCountRequirements`: total reads > 0 and >= n_samples/100) is lenient (98.2% would pass) -- the real bottleneck is its *recommended* default filter (FPKM>=1), which is far stricter and is what actually ships in its documented workflow.

## 2. Coverage by detection frequency (nz = # HC samples with count > 0, out of 676)

| nz bin | n genes | our engine: modeled | OUTRIDER: FPKM>=1 pass |
|---|---|---|---|
| 1-10 (near-absent) | 873 | 100% (pool route) | 0.0% |
| 11-30 (routed to pool, `NZ_A_MAX`=31) | 1,227 | 100% (pool route, 40 genes with nz in [11,30] individually convergent enough to stay on model route) | 0.0% |
| 31-100 | 2,206 | 100% (model route) | 1.3% |
| 101-300 | 2,972 | 100% (model route) | 20.1% |
| 301-676 (near-universal detection) | 12,580 | 100% (model route) | 94.7% |

**Interpretation**: our engine's coverage is flat at 100% across every detection-frequency stratum -- the pool-route cascade (`NZ_A_MAX` routing, `core/glmm_fit_pool.R`) exists specifically to rescue genes too sparse for individual GLMM convergence. OUTRIDER's FPKM filter is essentially a step function of detection frequency: it structurally cannot model the sparse/rare-detection tail of the transcriptome (nz<100, ~4,300 genes, 21.7% of the universe) at all, and even in the "near-universal" bin loses 5.3%. This is a genuine, literature-consistent difference (not a bug on either side): OUTRIDER's FPKM-based filter was designed for solid-tissue RNA-seq (GTEx-scale detection rates), and cfRNA's intrinsically sparser per-gene detection runs into that filter directly.

## 3. Fit-failure caveat (see also fairness discussion in conversation)
7 genes that pass OUTRIDER's FPKM filter still crash its per-gene NB optimizer (L-BFGS-B non-finite) inside `controlForConfounders`/`fit()`, which has no partial-success mode -- confirmed to fit cleanly in our engine (`route=model`, `ok=True`, `nbi_full_eb`, all with nz in [53,77]). These were excluded via a manual retry-drop patch not part of OUTRIDER's native workflow; a typical user hitting this gets a fully dead pipeline, not a silently-reduced gene set.

*(CV-level Z-moment/PPC comparison to follow once the corrected 5-fold OUTRIDER CV run completes.)*
