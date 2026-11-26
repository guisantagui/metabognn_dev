> **2026-11-25 (later) — the producer-average null and the residual both come back POSITIVE for MetChange.**
> Four arms, one design (HEPG2, Attene-Ramos MMP, n=144/475, all recomputed in the same run): on `tca` and
> `pyruvate` in MMP-down, MetChange sits at the random-set floor (AUROC 0.362–0.395, rand 0.005) while the
> **producer-average** — the same genes, aggregated to the same metabolites, with no LP and no stoichiometry —
> is **null** (0.422–0.490, rand 0.10–0.92), and gene-level expression is null too (0.405, rand 0.12). The
> **global-residual** arm keeps the signal (0.352–0.382, rand 0.010–0.040), so §8's one-dimensionality does
> not invalidate the report 09 §4 result. The producer-average is not inert — it is the strongest arm on the
> MMP-**up** contrast (nucleotide 0.259 at the floor) — so this is a real comparison.
> **BUT THE DIRECTION CONFLICTS WITH BOTH EXPRESSION ARMS AND IS UNEXPLAINED.** No sign flip is applied
> anywhere: AUROC 0.362 means MMP-down drugs have more NEGATIVE TCA log2FC, i.e. a GAIN of transcriptome
> consistency. Yet both expression arms say TCA genes go DOWN (producer-avg 0.422, expr 0.405), which should
> make TCA metabolites more expensive, not less. Two candidates — the LP's min-over-routes, and rank-relativity
> (penalties use within-sample ranks of ABSOLUTE expression, the expression arm reports log2FC, and the two
> can move oppositely by construction) — each with a diagnostic in report 11 §10.3b. **Do not build a
> mechanistic reading on the TCA result until one is settled.**
> **The statistical claim is the first result isolating the stoichiometric integration itself rather than
> MetChange vs expression**, and the disagreement arguably strengthens it (a proxy for TCA gene expression
> could not disagree with TCA gene expression) — but what MetChange tracks instead is currently unknown. Also established: the LP is what collapses the
> dimensionality (producer-average retains 5.18 effective dimensions, MetChange 1.10). On the MMP-**up**
> contrast the arms' sign inversion does reconcile them (producer-avg 0.259 vs MetChange 0.630 are the same
> biological statement); on MMP-**down** it does not, which is the conflict above. Gates remaining for the
> FNN: resolve §10.3b, and a grouped cross-validated target-prediction task.
> Launchers: `reproduce_report11.sh 9 10 11`. Full account: **report 11 §9–§11**.

> **2026-11-25 — the MetChange label matrix is effectively one-dimensional.** PC1 = 95.5% (HEPG2) and 95.0%
> (A549) of the feature-correlation variance; participation ratio 1.10 and 1.11; **one** eigenvalue carries
> 90% of the variance; random metabolite pairs correlate at r = 0.97. The same statistic on the *expression*
> arm of the same cells and perturbations gives 52–58 eigenvalues and r = 0.045–0.066, so the collapse happens
> in the expression→MetChange step, not in the data. The global component is **not** response magnitude
> (Spearman +0.14 against mean |expr log2FC|), so the `n_sig` covariate matching used in reports 08–10 does
> **not** target it, though the `z` statistic itself partially does (report 11 §8.3). Curated sets do nonetheless carry set-specific signal above size-matched random sets
> (`dntp_endpoints` at the 99.7th percentile), so the existing results are not void — but every set-level
> claim should be read against its `rand=` column, not its parametric p. Full account, the correction of a
> wrong intermediate figure, and the consequences for the FNN: **report 11 §8**. Script:
> `scripts/mc_phenotype_assoc/feature_dimensionality.py`.
>
> Also: folate is a dead MetChange feature (0th-percentile SD in both cells), so the NCI-60 design is **25
> metabolites, not 26**. And the medium does *not* explain the amino-acid null — alanine, the one amino acid
> absent from `rpmi_medium.json`, behaves like the 19 that are in it (report 11 §8.5).

# Project state — 2026-11-19

Handoff note for continuing in a new session. Everything described here is committed up to
HEAD `d11744b`, except the files listed in section 8.

> **2026-11-19 — validation datasets compiled and mapped (report 07).** `map_smiles.R` rewritten (salt
> stripping, all pairs kept, stereochemistry labelled not assumed). **Attene-Ramos (mito membrane potential,
> HepG2, ~10k compounds) is the best-powered external set: 904 LINCS HepG2 compounds mapped**, HepG2 MetChange
> labels already on scratch. **Hallinger (Seahorse) is limited to its 270 named tier-2 actives** — the 1,042-
> compound tier-1 screen was never registered in invitrodb (checked: analytical_qc, assay_annotations, pubchem
> per-aeid, sc1_sc2 from the 7.5 GB summary; 0/1,051 spids), so only the authors have the plate maps. Currency
> metabolites (ATP, NAD+/NADH, O2) **are** MetChange features but **not** network nodes: label-level tests can
> use them, model queries cannot. See section 6.2c.

> **2026-11-18 — recipe state (report 06, not frozen).** Best complete run: **expr `_qsigtop` 10.5 points**
> (q < 0.05 then top 100 by |log2FC|, 1,358 drugs, one seed). `_qsigtop` matches `_min15` for all three signature
> types (±2 points, p ≥ 0.12) and beats `_qsig` by 5–8, so **requiring significance is free**; figure:
> `results/predict_targets/A549/performance_figure/` (`scripts/predict_targets/plot_performance_panels.py`); **SD scaling is
> dropped** (`_sdscale` 6.7, `_qsigtop_sdscale` ≤ `_qsigtop`). MetChange `_qsigtop` 10.3 is provisional (fold 0
> complete since 2026-11-18, 10.2 points) and covers only the 758 drugs with significant metabolites, so **`_min15`
> stays the primary MetChange run**. One seed for every `_qsigtop` run. See section 6.0.

> **2026-11-17 (evening) — selection settled, Schuhknecht null explained** (report 05, sections 7–8).
> `_qsigtop` (q < 0.05, then top 100 by |log2FC|) equals the unfiltered top 100 for expr (10.5 points),
> TF (7.8) and MetChange (10.3, folds 1–4; fold 0 missing): significance can be required at no cost.
> SD scaling does not help (`_qsigtop_sdscale` 8.1, `_sdscale` 7.1 provisional). The Schuhknecht
> benchmark is null in both seeds and both controls. Cause: on the limma labels MetChange does not agree
> with measured Schuhknecht metabolomics for the same drug (drug-specific ρ −0.013, below the
> expression-producers baseline +0.033, p = 7e-6; specificity permutation p = 1.0).
>
> Earlier today: corrected MetChange matches expr and TF (report 05, section 5): per-drug gap 8.3 / 7.7
> (seeds 0 / 1), vs expr 7.8 / 10.4 and TF 9.0 / 9.3, not significantly different.

**2026-11-25 — NEXT ANALYSIS ASSESSED: NCI-60 exchange fluxes, see `reports/2026-11-25_nci60_exchange_flux_feasibility.md`
(report 11).** Feasible and the FIRST dataset that measures FLUX — the ground truth report 05 9.3 says
MetChange's claim requires and that pools cannot supply. **11 cell lines x 26 exchange reactions** (NOT 111
metabolites), already mapped to Human-GEM reaction ids. All 26 reachable in the MetChange feature space, all
with an EXTRACELLULAR species and no hub flags there. All 11 lines have CellMiner RNA-seq covering **2,482 of
2,522 GEM genes (98.4%) — better than the LINCS space the existing labels use (97.9%)**, so this would run
MetChange on direct RNA-seq rather than L1000-inferred values.
**THREE THINGS TO HANDLE: (1) `786-O` is `RE:786-0` in CellMiner (zero, not letter O) — normalised matching
silently drops it. (2) 1,103 duplicated symbols need aggregating; ranks make the log scale harmless but SUMMING
duplicates is only valid in linear space — convert first or use max. (3) CellMiner omits ALL 13 mtDNA-encoded
OXPHOS subunits (LINCS has all 13, so this is a NEW gap). It hits 5 reactions, all `Oxidative phosphorylation`
with AND-dominated complex GPRs (14-44 ANDs, 1-7 subunits missing of 10-45) — so it sets the penalty for
complexes I/III/IV/V. Impute at the median of the nuclear subunits, and run BOTH ways to report sensitivity.**
**DESIGN FIXED BEFORE RUNNING (report 11, 6): do NOT correlate across metabolites within a cell — fluxes span
SEVEN orders of magnitude, so that is report 10 8's confounder in extreme form. Use the between-cell
differential: 11 cells = 55 PAIRS vs PRISM's 10, so the saturated Wilcoxon stops binding. Split by DIRECTION
and predict OPPOSITE signs (8 released: more capacity -> more release; 18 consumed: more capacity -> LESS
uptake) — a crossed prediction like the GCP rationale.**
**COMPARATOR CORRECTED (report 11, 6.5-6.6): an FBA baseline (E-Flux/GIMME/iMAT/ec-GEM FVA) is the WRONG
comparator here and is WITHDRAWN.** That asks whether MetChange is a good flux-prediction method — the
MetChange paper's own question, already published. **The question here is whether the integration justifies
MetChange as the FNN's metabolic SIGNATURE, so the comparator must be whatever would otherwise be fed to the
FNN: the transcriptional signature.** E-Flux/GIMME/iMAT also fail on output space — they are REACTION-indexed
while the bipartite network propagates from METABOLITE nodes, so they are not substitutable inputs.
**THE RIGHT COMPARATOR IS A METABOLITE-CENTRIC EXPRESSION NULL: per metabolite, the mean expression of enzymes
catalysing reactions that PRODUCE it.** That is what MetChange reduces to if the LP and network-wide
stoichiometric coupling add nothing over local gene averaging. Same output space; **far tighter than the
current expression arm, which uses SUBSYSTEM gene sets** (`nucleotide gene` spans a 192-metabolite subsystem);
and it isolates the actual claim. `build_target_relevance.py` already derives producing routes, so the hop-1
producers per metabolite are a by-product of existing code. **It also retro-sharpens report 10 section 12: the
104% retention was measured against a subsystem gene set that was itself flat (-0.008); against a
per-metabolite producer average the test is much harder and correspondingly more convincing if it survives.**
**CAUTION: 16 of the 26 are AMINO ACIDS, so this is largely a test of the arm that FAILED in report 10 section
13. Either that failure was the contaminated drug group (13.4: ~95% non-amino-acid drugs) or amino-acid
capacity genuinely is not predicted — this distinguishes them, which is the dataset's main value.**
**The other two candidates: To et al. 2019 K562 CRISPR deposited NO expression data (verified — GEO returns 0
datasets for PMID 31730859; K562 also has only 136 compounds in LINCS) — DEPRIORITISE. Kang et al. 2024
(PMC11310444) is NOT reprogramming but contact-inhibition QUIESCENCE in MOUSE fibroblasts, one P-vs-Q contrast
— worth a cheap pre-registered directional case study (glycolysis down, TCA up, MPC up, ME1 route down, MDH1
route reversed) but needs Mouse-GEM and is n=1 contrast.**

**Reproducing any reported number:** see `scripts/mc_phenotype_assoc/REPRODUCING.md` — maps each report
section to the script that produces it, and names the three known gaps (`confounder_checks.py` /
`dose_exposure_check.py` for report 08's van Elteren MoA strata, potency and dose analyses;
`signature_conservation.py` for report 09 §6.3; `build_xu_phenotype.py` for report 10 §5).

**Read first:** `reports/2026-11-24_prism_gcp_cellpainting.md` **section 12 — the pre-registered conditional
test is CONFIRMED: on 5 NON-HepG2 cell lines, MetChange's cell-specific prediction tracks cell-specific drug
sensitivity (mean r +0.291, 10/10 pairs) and the matching GENE set carries NONE of it (-0.008, 4/10);
MetChange-given-expression +0.302, retention 104%. First measurement of INCREMENT in the project. The 28-drug
set is still post hoc (section 11.4), and section 12.3 shows 93 more PRISM cell lines are in LINCS without
labels — 20 cells would give 190 pairs instead of 10.** Then section 10 — **the first positive PRISM result:
with the biomass sets and TARGET-based drug selection, the NUCLEOTIDE arm behaves as predicted
(`dntp_endpoints` 8/10 cell pairs positive, matched-vs-complement paired Wilcoxon p = 0.002, complement ~ -0.05)
while the AMINO-ACID arm contradicts (2/10 positive, mean r -0.10). Not significant after honest correction and
partly a forking path — read 10.3 before quoting it.** Then the rest of report 10 — PRISM enables the first
well-powered cell-matched test, via a within-drug between-cell DIFFERENTIAL design that cancels the viability
confounder instead of matching it; GCP has the best rationale and a library that defeats it; the Way data is
not actually downloaded), then `reports/2026-11-23_cross_cell_specificity.md` (report 09 — the association is HepG2-specific
across 8 cell lines on identical drugs, which rules out the chemotype explanation; the expression arm is NOT
cell-specific, which is the real argument for the stoichiometric layer), then
`reports/2026-11-20_phenotype_confounders.md` (report 08 — are the report 07 phenotype
associations chemistry or metabolism: the ToxCast-subset finding that explains Hallinger's null, mechanism-class
stratification, the absent potency gradient, and what may/may not be claimed), then
`reports/2026-11-19_validation_datasets.md` (report 07 — external validation sets, structure
mapping, what each can test), then `reports/2026-11-18_recipe_state.md` (report 06 — where the recipe stands, what is settled and
what is pending), then `reports/2026-11-17_signature_selection_scale.md` (report 05 — q-values unused; MetChange's
top-100 is chosen by metabolite scale; `_qsig` / `_sdscale`; section 5: corrected MetChange results), then `reports/2026-11-16_labels_network_coverage.md` (report 04 — the label bug in the
MetChange runs, expr on `combined`, full coverage check), then `reports/2026-11-14_seed_tf_expr.md` (report 03 — second seed, MetChange vs expr
vs TF, early stopping and lr schedule; sections 9–10 have the current recipe, the two-seed
result and the metabolic-coverage check), then `reports/2026-11-14_lr_noglobals_matched.md` (report 02 —
the current recipe), then `reports/2026-11-13_A549_matched_negatives.md` (report 01 — metric
definitions, the dose/time leak, how to read a comparison). `scripts/predict_targets/README.md`
documents every launcher and how to submit runs with their controls.

## Goal and standing constraints

Adapt the Menacho et al. 2026 GNN drug-target framework to LINCS signatures. A549 first,
then the 9 core cell lines (A375, A549, HA1E, HCC515, HEPG2, HT29, MCF7, PC3, VCAP).

- Do not commit; the user commits.
- Keep gene symbols in the prediction pipeline.
- KEGG-derived edges must not be redistributed.
- Uncharacterized drugs are not negatives (no annotation is not evidence of no target).
- **Priority is MetChange**; expr is a welcome byproduct. TF reproduces the published Menacho
  setting and is the **positive control** (the level to reach), not a result.
- Compare every run against **its own shuffled-signature control**, never against chance, and
  test **per drug**, not per perturbation (doses of one drug are not independent).
- Long write-ups go in `reports/`, not into chat (see section 9).

## 1. Data

LINCS L1000, A549, from the Ma'ayan lab **CycleGAN** release: transcriptome-wide, but most
genes are **inferred** from the 978 landmarks. Consequence to keep in mind: the inferred genes
give the coverage MetChange needs over the GEM's GPR genes, but add no independent
measurement — the information ceiling is set by the 978 landmarks, and imputation smooths
toward typical profiles. Two cheap, unrun checks that would size that ceiling (section 6.3).

## 2. Labels

**limma reference-level sign bug — fixed.** `model.matrix(~ perturbation)` took whichever of
the perturbation name and `DMSO_nan_<time>` sorted first as the reference, so for ~55% of
A549 perturbations limma returned `control - treated`. Correction: `log2FC -> -log2FC`,
`basal -> basal + log2FC` (`fix_limma_reference_sign.py`), with R's own level order
(`limma_reference_order.R`) deciding, not the data heuristic. Corrected expression labels are
in `expr/`; PC3 was already correct; MetChange has `rank_signfixed`. `limma_DE.R` now sets the
control reference explicitly and writes `reference_level.txt`, so new runs need no correction.

**But the home MetChange copy was never corrected** (report 04, section 1; checked 2026-11-16 with
`compute_labels/compare_label_dirs.py`). `<home>/results/labels/A549/metchange/rank/` vs the scratch
rerun: 19,884 perturbations identical, 24,512 exactly negated (r = −1, all features), and the
negated set is exactly `flipped_R = TRUE`; q / p-values and significant sets identical. The
correction was only ever copied into scratch (`move_corr_files_mc.sh`). Scratch is correct by
construction (fixed `limma_DE.R`); **confirmed 2026-11-16**: its `reference_level.txt` reads
"reference level: ctrl (explicit)" (written only by `limma_DE.R` ≥ `85526fe`, which includes the fix). The expr pair of that
comparison wrote no output (skipped or failed; `logs/cmp_labels.out` not yet looked at).

**MetChange limma rerun over the 9 cell lines — in progress, verify.** The jobs segfaulted in
`fread` (data.table crashes on the ~150k-column MetChange matrix). Fixed in
`launch_limma_DE_generic_mc.sh` (commit `de76198`): `mc` mode now passes `--reader lines`, and
output goes to `labels/<cell>/metchange/<penalty>/` instead of `labels/<cell>/metchange/`.
Resubmit with `launch_limma_mc_celllist.sh` and check the first `.out` for the
"N features read" progress lines.

**MetChange concatenation.** `mc_res_concat_list_launcher.sh` is a Slurm array, one task per
cell line (`--array=1-9`, `-c 8`, `POLARS_MAX_THREADS`, skips cell lines that already have a
`*_concat.csv` unless `overwrite=true`). Note: `concat_mc_res.py`'s docstring claims it
computes deltas against matched controls and averages per perturbation; the code only
concatenates batches and drops all-missing metabolites. Only the description is wrong.

## 3. Where the model stands

Two problems were found and are now understood:

1. **Per-protein prior.** The model learned which proteins tend to be drug targets, which does
   not depend on the signature. `spec` / `rank_probc` cancel it (report 01, section 3).
2. **Dose/time shortcut.** In LINCS the experimental design tracks the screen, and so the drug
   class. The FNN exploited it through the `rel_log_dose` / `log_time` inputs and their
   products with protein features, which inflated `spec` even in shuffled runs (report 01,
   sections 5.3–5.5). Fixes: `--no-globals`, and the dose/time-stratified `spec` in
   `compare_runs.py`.

**Current recipe** (report 03, section 9): `--no-globals --lr 0.001 --min-epochs 15
--patience 15 --ffn-epochs 80` (StepLR), unmatched, epochs selected on AUPRC. The `_min15` runs.

**Gap convention** (report 03, section 1): per-drug gap = `median_diff` of `paired.tsv`, the
median of the per-drug real − control differences, in percentile points; it is what the
Wilcoxon p tests. Report 02's "5.6 points" was a difference of medians (paired: 4.2).

Per-drug gap (points), p — real vs own control, current recipe, seeds 0 / 1. MetChange on the
corrected labels (report 05, section 5); expr / TF from report 03, 9–10:

| metric | MetChange | expr | TF (positive control) |
|---|---|---|---|
| `spec_prob_timedose` (leak-free) | 8.3 / 7.7 (6e-12 / 4e-11) | 7.8 / 10.4 (2e-11 / 1e-14) | 9.0 / 9.3 (8e-15 / 7e-15) |
| `rank_probc` | 9.3 / 7.7 | 8.0 / 11.5 | 10.4 / 10.3 |
| `rank_prob` | 2.5 / 2.5 (1e-11 / 5e-9) | 3.3 / 3.6 | 3.3 / 3.3 |

Paired by drug, MetChange vs expr: p = 0.79 / 0.085; vs TF: 0.21 / 0.12; expr vs TF: 0.36 / 0.68.
**All three signatures are indistinguishable.** (Sign-buggy MetChange was 4.2 / 2.7, 6.4 / 4.4,
1.3 / 0.7; corrected vs buggy, seed 1: p = 2.5e-5.)

- **The warm-up (`--min-epochs 15`) fixed the collapse**: no FNN keeps epoch 1–2 any more
  (before: 7–8/15 for TF and expr). Expr was training-limited and now matches the positive
  control, all five folds at 0.30–0.32 (`spec_prob_timedose` per perturbation); TF also gained
  (6.8 → 9.0).
- **Schedule variants** (patience 25, `ReduceLROnPlateau`) made no difference on sign-buggy
  MetChange, nor on expr / TF; keep StepLR. Not retested on corrected MetChange.
- **Resolved: "MetChange reads the signature half as well as expr"** (reports 03–04) was the
  label sign bug. Network ruled out independently (below).
- **Expr on the `combined` network performs as well as on `gene`** (report 04, section 2; seed 0,
  same folds, `expr_combined_noGlobals_lr1e-3_min15`). Per-drug gap over own control:
  `spec_prob_timedose` 7.7 (p = 3e-12) vs 7.8 on `gene`, `rank_probc` 8.2 vs 8.0, `rank_prob` 3.0
  vs 3.3; real better in 5/5 folds. The MetChange network is not why MetChange lags expr.
- **Robust across two seeds: MetChange ≈ expr ≈ TF** (corrected labels). Per-drug gaps agree
  only moderately between seeds (MetChange Spearman 0.32): read effects in aggregate, not drug by drug.
- **Metabolic coverage** (full data, report 04 section 3): in DE terms it is not the issue
  (every active perturbation has ~14–15% of its DE genes in Human-GEM, the background share).
  In MetChange's own output it is all-or-nothing: 57% of active perturbations have no
  significant feature among the 100 the model reads, 35% have all 100. Corrected MetChange
  (report 05, 5): 14.6 / 14.7 points on fully covered perturbations (expr on the same: 12.1 /
  10.5), 5.4 / 5.6 on uncovered ones (expr 7.3 / 11.8). MetChange beats expr where it has a
  response; uncovered perturbations are what pull its average down to the expr level.
- The epoch cap does not bind for MetChange; for TF / expr 8–12/15 members reach epoch 80, but
  at an lr ≤ 1.6e-5, where more epochs change little.

**Matched negatives are superseded.** With dose/time inputs gone and training to convergence,
matched sampling is worse on both metrics (report 02, section 4). What looked like
matched-negative specificity in report 01 was mostly the dose/time shortcut plus `spec`-based
epoch selection. Version (b) (grouped specificity loss) is demoted accordingly.

**Plain ranking:** with corrected labels the MetChange signature moves per-drug `rank_prob` by
2.5 points (0.381 vs 0.415 control), close to expr / TF (3.3–3.6). The prior still dominates the
absolute ranking for all three (real 0.38–0.40 per drug vs 0.5 chance).

## 4. Code in `scripts/predict_targets/`

`metabolic_target_prediction.py`:

- Prior-free metrics in `metrics.json` / `cv_metrics.json` / rank tables: `spec_prob_heldout`,
  `spec_logreg_heldout`, `rank_probc_heldout`, `rank_logregc_heldout`; `test_scores_heldout.npz`
  per fold (every test perturbation's score for every held-out drug).
- `--no-static` (drop per-protein features) / `--no-globals` (drop dose, time, signature size).
  **Different ablations** — `noStatic` keeps dose and time.
- `--score-active-only`, `--min-de-genes` with `--de-counts`.
- `--matched-negatives K`, `--matched-weight`, `--select-metric {auprc,macro_f1,spec,auprc_matched}`;
  `epoch_selection` and `per_fold_{spec_val,auprc_matched,best_epoch}` in the metrics.
- `--min-epochs N` (default 1 = plain early stopping, previous behaviour; verified by replaying
  the 150 recorded member histories): epochs < N are trained and logged but can neither be kept
  nor stop training; `--patience` counts from N. Recorded as `min_epochs` in `epoch_selection`.
- `--lr-schedule {step,plateau}` (default `step` = `StepLR`, `--lr-step-size 10 --lr-gamma 0.5`,
  as before); `plateau` = `ReduceLROnPlateau` on the validation selection score from epoch 1,
  `--plateau-patience` (3, must be < `--patience`), `--min-lr` (1e-6). Each epoch's lr is in
  `history`. Checked with a stub-torch test of `train_ffn` only; first real run is on the
  cluster (report 03, 8.2).
- `--fold i` with `--n-folds 5` for array jobs; the last task pools (`--aggregate-only` to pool
  by hand).

Tooling:

- `compare_runs.py` — summary / per fold / paired per-drug Wilcoxon / ECDFs, plus `spec`
  **stratified** by time and by time × dose, with mid-rank percentiles (unbiased for small
  strata; the plain formula averages 0.58 at n = 5 — report 01, section 5.4).
- `check_metabolic_coverage.py` — per-perturbation metabolic coverage (expression DE in
  Human-GEM; MetChange's own significant features and significant features among the top-n
  read) and the run's spec split by it, vs its control (report 03, 10.2–10.3).
- `check_dose_time_leak.py` — `spec` of the model vs `spec` of a dose/time-only predictor, and
  the score variance explained by dose/time groups.
- **Labels (2026-11-15):** every training launcher reads its labels from scratch,
  `${project_dir_scratch}/results/labels/<cell>/{expr,metchange/rank}` (MetChange used the home
  copy until then). The scratch copy is the rerun of the fixed `limma_DE.R`; the home copy was
  ~~corrected with `fix_limma_reference_sign.py`~~ — **wrong for MetChange**: the home
  `metchange/rank` is the uncorrected output (report 04; the earlier "spot check confirms it" was
  not true). The expr launchers always read scratch.
  Label pipeline, all writing to scratch: `parse_lincs_data/launch_get_lincs_celllist.sh`
  (expression limma) → `metchange/metchange_list_launcher.sh` → `metchange/mc_res_concat_list_launcher.sh`
  → `compute_labels/launch_limma_mc_celllist.sh`.
- **Output location (2026-11-15):** every launcher writes to
  `/scratch/users/gsantamaria/test_large_files/metabognn_dev/results/predict_targets/A549/<out_dir>`
  (`${project_dir_scratch}`), and `submit_with_shuffle_control.sh` /
  `compare_runs_mc_vs_shuffled_A549.sh` read from there. Earlier runs are in the same folders
  under `/home/users/gsantamaria/projects/metabognn_dev/results/predict_targets/A549/`.
- `submit_with_shuffle_control.sh` — submits a variant's real run, its shuffled control and a
  dependent comparison; takes several variants in one call, validates names first, `--real-job`
  for an already-queued real run. Needs bash ≥ 4 (cluster yes, macOS no).
- Wrappers: `submit_train_w_shuff_mc_matched.sh`,
  `submit_train_w_shuff_mc_noGlobals_matched_auprc.sh`,
  `submit_train_w_shuff_mc_matched_noGlobals.sh`,
  `submit_train_w_shuff_mc_noGlobals_lr1e-3.sh`,
  `submit_train_w_shuff_mc_noGlobals_lr1e-3_ep100_seed1.sh` (done, report 03),
  `submit_train_w_shuff_tf_expr_noGlobals_lr1e-3.sh` (done, report 03),
  `submit_train_w_shuff_noGlobals_lr1e-3_min15.sh`,
  `submit_train_w_shuff_mc_noGlobals_lr1e-3_pat25_plateau.sh`,
  `submit_train_w_shuff_tf_expr_noGlobals_lr1e-3_plateau.sh` (done, report 03),
  `submit_train_w_shuff_noGlobals_lr1e-3_min15_seed1.sh` (done, report 03 section 10),
  `submit_train_w_shuff_expr_combined_noGlobals_lr1e-3_min15.sh` (6.1).
- Variants registered in the submit script: `tf`, `expr`, `mc`, `mc_matched`,
  `mc_matched_auprc`, `mc_matched_noGlobals`, `mc_noGlobals`, `mc_noGlobals_lr1e-3`,
  `mc_matched_noGlobals_lr1e-3`, `mc_noGlobals_lr1e-3_ep100`, `mc_noGlobals_lr1e-3_seed1`,
  `tf_noGlobals_lr1e-3`, `expr_noGlobals_lr1e-3`, `mc_noGlobals_lr1e-3_min15`,
  `tf_noGlobals_lr1e-3_min15`, `expr_noGlobals_lr1e-3_min15`, `mc_noGlobals_lr1e-3_min15_pat25`,
  `{mc,tf,expr}_noGlobals_lr1e-3_plateau`, `{mc,tf,expr}_noGlobals_lr1e-3_min15_seed1`,
  `expr_combined_noGlobals_lr1e-3_min15`, `{mc,tf,expr}_noGlobals_lr1e-3_min15_qsig`.
  `mc_noStatic` is not registered
  (launch both by hand).
- The `_lr1e-3` launchers are **array jobs** (`--array=0-4`; 4 h per task, 6 h for `_ep100` and
  the `tf` / `expr` ones, 8 h for `_min15`, `_min15_pat25`, `_plateau`). `sbatch --test-only` prints the estimated start time without
  submitting.

## 5. Runs in `results/predict_targets/A549/`

Analysed (in reports): `metchange_rank{,_shuff}`, `metchange_rank_matched{,_shuff}`,
`metchange_rank_noStatic{,_shuff}`, `metchange_rank_noGlobals{,_shuff}`,
`metchange_rank_matched_auprc`, `metchange_rank_noGlobals_lr1e-3{,_shuff}`,
`metchange_rank_matched_noGlobals_lr1e-3{,_shuff}`, and (report 03)
`metchange_rank_noGlobals_lr1e-3_{ep100,seed1}{,_shuff}`, `tf_active200_noGlobals_lr1e-3{,_shuff}`,
`expr_noGlobals_lr1e-3{,_shuff}`. Comparisons live in each real run's
`compare_real_vs_shuffled/`; cross-run comparisons in `mc_matched_vs_unmatched/`,
`mc_real_vs_suff/`, `stratified_spec/`.

Not on the laptop: `metchange_rank_matched_auprc_shuff` (its comparison ran on the cluster),
and `metchange_rank_matched_noGlobals{,_shuff}` (plain, lr 0.005) — launched 2026-11-13, not
yet synced or analysed. It is now of low value: report 02 supersedes it.

Analysed in report 03, section 9 (seed 0, same splits as `*_noGlobals_lr1e-3`):
`metchange_rank_noGlobals_lr1e-3_{min15,min15_pat25,plateau}{,_shuff}`,
`expr_noGlobals_lr1e-3_{min15,plateau}{,_shuff}`,
`tf_active200_noGlobals_lr1e-3_{min15,plateau}{,_shuff}`. Report 03, section 10 (seed 1):
`{metchange_rank,expr,tf_active200}_noGlobals_lr1e-3_min15_seed1{,_shuff}`. Report 04:
`expr_combined_noGlobals_lr1e-3_min15{,_shuff}` (seed 0), `metabolic_coverage_full{,_topnbins}/`.
All MetChange runs listed here used the sign-buggy home labels (report 04).
Corrected (report 05, section 5): `metchange_rank_noGlobals_lr1e-3_min15{,_seed1}{,_shuff}` (the laptop
copies are now the corrected runs; the seed-0 buggy pair exists only on scratch as `_signbug`),
`metabolic_coverage_fixedlabels/`. Report 05, 5.2: `metchange_rank_noGlobals_lr1e-3_min15_qsig{,_shuff}`. On scratch, `metchange_rank_noGlobals_lr1e-3_min15{,_seed1}{,_shuff}` were renamed to
`..._signbug{,_shuff}` (2026-11-16), and `..._min15_pat25{,_shuff}` were deleted 2026-11-17 (laptop
copies remain).

Pre-dating the current metrics, cannot be re-checked (no `test_scores_heldout.npz`, no FNN
shuffle control): `expr`, `tf`, `tf_active200`, the `*_linOnly_shuff` runs.

## 6. Next steps, in priority order

**6.0 Running the phenotype association (refactored 2026-11-20, report 08 section 8).**
`sbatch scripts/mc_phenotype_assoc/launch_assoc_HEPG2_slurm.sh` (or `..._A549_slurm.sh`) -> exports `cell` ->
`launch_assoc_slurm_generic.sh` (2-task array: 0 MetChange, 1 expression) -> `launch_assoc_metchange.sh` /
`launch_assoc_expr.sh` (no cell line in either) -> `mc_phenotype_assoc_common.sh`. One-off per cell:
`cell=<CELL> bash scripts/mc_phenotype_assoc/prepare_assoc_inputs.sh`, which needs
`lincs_unique_perts_<CELL>.tsv` plus `cell=<CELL> bash scripts/benchmarking/map_{attene,hallinger}_lincs_perts.sh`.
Paths (all under scratch on the cluster, repo on the laptop): phenotypes **per cell** at
`results/mc_phenotype_assoc/<cell>/phenotypes/`, sets **per project** at
`results/mc_phenotype_assoc/metabolite_sets/`, outputs at `results/mc_phenotype_assoc/<cell>/`. Sets are now
built WITHOUT `--features` so they are genuinely cell-independent; absent members are dropped at scoring time
and reported. Verified equivalent on HEPG2 (same 40 sets, same numbers). **An A549 run is the cross-cell-line
CONTROL, not a replication** - the phenotypes are measured in HepG2 regardless of which cell line is mapped.
The slurm wrappers are **submittable from anywhere in the repo** (they walk up from `SLURM_SUBMIT_DIR` to find
the root and export `repo_dir`; `$0` is unusable there because Slurm spools the script) and log to
`./mc_pheno_<cell>_%A_%a.out` in the submit dir. **Bug found and fixed 2026-11-20 while testing the A549
path:** `build_phenotype_tables.py` fell back to the legacy unsuffixed mapping for ANY cell, so `cell=A549`
silently built A549 tables from the HEPG2 mapping (identical 849/147/480 counts). The fallback is now allowed
only for HEPG2; other cells require `mapping_<cell>/` and the run stops with the commands to build it.

**Running ALL cell lines (report 08, 8.4-8.5).** `bash scripts/mc_phenotype_assoc/submit_assoc.sh` discovers
which cell lines actually have limma labels (`check_label_availability.py` over `core_cell_lines.tsv`; checks
`log2FC.csv` required, `qvals.csv` / `fit_info_nsig.csv` optional, and warns when neither exists because then
`--mode sig` AND the responsiveness stratification are unavailable), writes a `(cell, arm)` manifest, sizes
the Slurm array to it and submits `launch_assoc_manifest_slurm.sh`. `cells="HEPG2 A549"` restricts,
`dry_run=true` previews. Cell lines missing labels or phenotypes are skipped, not failed. **Logs land in
`scripts/mc_phenotype_assoc/logs/`** (kept in git: Slurm opens `--output` before the job runs, so a mkdir
inside the script is too late); `submit_assoc.sh` passes that path absolutely, the single-cell wrappers use a
relative `logs/` so submit those from `scripts/mc_phenotype_assoc`.
**Project root is configured in ONE place:** `scripts/mc_phenotype_assoc/project_paths.conf` (candidate
`scratch_proj_dir` lines, first existing wins, repo as fallback), sourced by common and submit_assoc and
parsed by check_label_availability.py. It exists because the path had been hardcoded in all three and the
checker's copy defaulted to the repo when `$labels_base` was unset - so on the cluster it read `$HOME` and
reported "no directory" for every cell line whose labels are on scratch. Fixed 2026-11-20.
**One mapping instead of nine:** `cell=ALL bash prepare_assoc_inputs.sh` maps against
`lincs_unique_perts.tsv` (the union; all 5,284 HEPG2 pertnames verified inside it) and writes
`<assoc_base>/phenotypes_ALL/`, which is cell-independent because the association script keeps only the drugs
with labels in the cell line being run. Common prefers a per-cell dir and falls back to the union, printing
which. Union mapping: Hallinger 176 compounds vs 115, Attene 1,865 vs 918; union tables hold 1,763 Attene and 177
Hallinger drugs. **NOT a drop-in replacement** - restricted to drugs with HEPG2 labels the union table gives
144 vs 147 MMP-down, 475 vs 480 inactive, 27 vs 29 Hallinger-active, because mapping against 34,419 perts
instead of 5,494 lets compounds match further pert_ids and a few aggregated labels resolve differently.
**So for ANY cross-cell comparison use `use_union=true bash submit_assoc.sh`**, which forces one identical
table onto every cell line; otherwise cells mix per-cell and union tables and a 2-3% assignment difference
would look like a cell-line effect. The per-cell HEPG2 table stays the one reproducing reports 07-08.
**Data note:** `lincs_unique_perts.tsv` originally had 4 header fields for 5 data columns (no `V1`), which R
tolerates as row names but silently misaligns every Python `DictReader`. Regenerated with `V1` 2026-11-20 -
check other files written by that code path.

**6.0a-arch (2026-11-20, report 08 section 9.1) — THE DOWNSTREAM PASS IS DEAD FOR METCHANGE.** `catalysis`
runs protein -> reaction with **zero** reaction -> protein edges (23,098 vs 0; the only metabolite -> protein
route is `second_messenger`, 59 edges). Measured on the combined network, a 100-metabolite signature reaches
**0 of 18,760 proteins at downstream hops 1-4 and 4 (0.02%) at hop 5**, against upstream 1.0% / 4.4% / 38.9% /
**64.8%** at hops 2-5. Consequences: (i) consumer edges (metabolite -> reaction, 20,373 of 40,786 metabolic
edges) carry NOTHING to proteins; (ii) **every `down_*` dyn/diff/slog feature plus `up_*_hop1` is identically
zero in MetChange mode** - about half the dynamic block - while expr/TF inject at protein nodes so both passes
work from hop 1, i.e. the two signature types are not using the same model capacity; (iii) the informative
signal sits at upstream hops 2-3 where only 1-4% of proteins are reachable, a plausible mechanical reason why
the no-network kNN beat the model (report 06, 6.1). **Cheap untried experiment: drop the dead features for
MetChange runs and re-test whether the network adds anything.**

**6.0a DO THIS FIRST — finish cleaning the lumped / placeholder species out of the network.**
The existing filter (a hand-made string-match list of aggregated species, plus the currency metabolites)
catches the obvious cases but has gaps. Found on 2026-11-20 while chasing the O-glycan result
(report 08, 3.3): **`Protein-Linked Serine Or Threonine Residue (O-Glycosylation Site)` (MAM03395g) is still a
metabolite node**, as are the **lumped glycan classes** `core 2`-`core 8`, `T-antigen`, `Tn-antigen`,
`sialyl-T`/`sialyl-Tn`/`disialyl-T antigen`, `glcnac-alpha-1,4-core 1/2`, `F1alpha`, `G00031`, `G00032`.
The three true `[protein]` entities ARE correctly absent, so the filter works, it is just incomplete.

Why it matters: these are not single chemical species, so a MetChange "production capacity" for them has no
measurable referent; they are low-degree, so their value is close to a readout of a few local enzymes; and
they can be selected as sources when significant and large. In report 08, 3.3 a 21-member set of exactly
these carried an association as strong as the TCA one (per-member median |AUROC-0.5| 0.097 vs 0.100).

A second, name-based candidate list is in **`data/networks/combined/lumped_candidates_for_review.tsv`**
(560 distinct names over 772 metabolite nodes, 9.8%, by category) — useful only for catching aggregates that
*do* have a clean formula, e.g. the lipoprotein particles. **Its patterns over-match on purpose — it is a
review list, not a deletion list** (it flags things like a tetrahydrobenzo-phenanthro-oxirene-diol purely for
containing a bracket). Categories in it where the aggregate reading is clear:

| category | nodes | examples |
|---|---|---|
| glycan core / antigen classes | 66 | `core 2`-`core 8`, `T-antigen`, `Tn-antigen`, `I-Antigen`, `F1alpha` |
| KEGG glycan ids | 37 | `G00005`, `G00020`, `G00031`, `G00032` |
| lipoprotein particles | 15 | `HDL`, `LDL`, `VLDL`, `IDL`, `chylomicron`, `* remnant` |
| protein-linked residues / sites | a few of the 49 | the O-glycosylation-site residue above |
| explicit pool, polynucleotide | 5 | `Gm4-Pool`, `polynucleotide` |

Needs judgement, probably KEEP (chemically defined despite the pattern match): `[ACP]` thioesters
(`(2E)-decenoyl-[ACP]`...), the GAG and glycolipid series with specified linkages
(`III3Fuc-nLc4Cer`, `chondroitin sulfate A (GalNAc4S-GlcA)`...), chain-specified lipids
(`cholesterol-ester-10,13,16-docosa`).

**Criterion to use — verified, not just suggested: the Human-GEM molecular formula.** Keep a node if it
denotes a single chemical species with a fully specified formula; drop it if the formula is empty, contains an
`R` (unspecified side chain) or an `X` (attachment point), or is a `(...)n` repeat unit. The formula is in
`Human-GEM.yml` per metabolite (NOT in `metabolites.tsv`, which has no formula column).

**This catches every culprit found on 2026-11-20**, all via the `X` atom: `T-antigen` C14H24NO10**X**,
`Tn-antigen` C8H14NO5**X**, `core 2` / `G00032` / `F1alpha` / `glcnac-alpha-1,4-core 1` C22H37N2O15**X**,
`core 6` / `G00031` C16H27N2O10**X**, `sialyl-T` C25H40N2O18**X**, `disialyl-T` C36H56N3O26**X**, and the
O-glycosylation-site residue whose formula is simply **`XH`**. 11 of 11 flagged.

**IMPLEMENTED 2026-11-20 in `build_network.py` (report 08, 9.2-9.3); THE NETWORK IS NOT YET REBUILT.**
`formula_pseudo_species()` applies exactly this criterion (`--no-formula-filter` restores the old behaviour):
238 curated names -> 642, i.e. 599 -> 1,264 metabolite entries (7.1% -> 14.9%), and **all 11 report-08
culprits are caught where the curated lists caught 0**. Open judgement call: 410 of the 665 new exclusions are
`R`-flagged generic lipids (`1-alkyl-2-acylglycerol`, phosphatidylinositols) - dropping the `R` branch still
catches all 11 culprits. Also added `CURRENCY_PRODUCER_ONLY = {ATP, NAD+, NADPH}`
(`--currency-producer-only`): kept as nodes with producer edges ONLY - verified 103 / 223 / 266 in-edges and
**0 out-edges** each, graph 19,894/38,564 -> 19,916/39,156. Safe because a one-directional node can never be a
path intermediary (no `R1 -> C -> R2` in either pass), and consumer edges are inert anyway per 6.0a-arch.
Chosen as the cheapest: 87 / 148 / 158 producing reactions, against ADP 792, NADP+ 797, NADH 456, so the
partners stay excluded. **PARP and the sirtuins remain unreachable from an NAD+ query** - that follows from
6.0a-arch, not from this policy, and my earlier "split each currency node into producer and consumer halves"
suggestion is withdrawn because the consumer half would be inert.

**`data/networks/combined/lumped_by_formula.tsv`** is that list applied to the network: **659 of 7,841
metabolite nodes (8.4%), 398 distinct names** — 410 with an `R` atom, 246 with an `X`, 3 with no formula.
Prefer this file over the regex list below; it is principled and complete where the string match was neither.

Decide separately what to do with the `R`-atom group, which is a different kind of case: generic lipids such
as `1-alkyl-2-acylglycerol`, `1-phosphatidyl-1D-myo-inositol-3-phosphate`, `1-(1-alkenyl)-sn-glycerol`. They
are unspecified in composition for the same reason, so the same argument applies, but there are 410 of them
and dropping them changes the lipid coverage of the network substantially.

**Knock-on effects to plan for before running it:** the label matrices do NOT need recomputing (MetChange is
computed per metabolite and limma ran over the whole matrix, which is why these appear in `log2FC.csv` at
all), but **`select_sources` output changes, so model runs and any analysis keyed on node identity would need
redoing**. Decide the criterion first, apply once, and note the network version in the run metadata.


**6.1 Corrected MetChange: done** (report 05, section 5: MetChange ≈ expr ≈ TF, both seeds).
Leftovers:
1. Read `compute_labels/logs/cmp_labels.out` for why the expr pair of the label comparison wrote
   nothing (fine if skipped as same provenance; the expr runs read scratch).
2. Check other consumers of the home MetChange copy (similarity checks, `metpath`,
   `mc_schuhk_assoc`) and the other cell lines' MetChange labels (`compare_label_dirs.py`).
3. The recipe parts established only on buggy MetChange (lr 0.001 vs 0.005, matched negatives
   "worse") are low priority now that MetChange matches the positive control.

**6.1b Signature selection** (no MetChange-vs-expr gap left to explain; these may still improve it).
- Network: **done, ruled out** (report 04, section 2: expr on `combined` 7.7 vs 7.8 points on `gene`).
- Coverage: **done, real** (reports 04–05): 57% of active perturbations have no significant
  MetChange feature among the 100 read. Corrected MetChange: 14.6–14.7 points on fully covered
  perturbations (beats expr there), 5.4–5.6 on uncovered. `_qsig` is in effect the MetChange
  activity filter (keeps 1,791 of 3,970 perturbations, 758 drugs).
- **MetChange `_qsig` done: worse** (report 05, 5.2): 6.1 points vs 11.9–14.2 for `_min15` on the same
  perturbations (p ≤ 3e-3). Half the training data dropped + random 100 of thousands of significant
  features. Magnitude-based top-100 stays; `_sdscale` is the test of scale.
- **Expr `_qsig` done: signal gone** (report 05, 5.3): 1.6 points (p = 0.22) vs 7.8 / 10.4 for `_min15`,
  with no perturbation dropped (all have ≥ 200 significant genes). Random 100 of the significant
  features loses the signal; the largest changes carry it. TF `_qsig` still running.
- **Significance (launchers ready, 2026-11-17):** every run so far, MetChange, expr and TF, used
  `--ignore-qvalue --top-n 100` (top 100 by |log2FC|, q-values unused; added 9 Nov in `5a28a0f`
  with the activity filter, no recorded reason). The published setting is significant features,
  capped at 100. `bash submit_train_w_shuff_noGlobals_lr1e-3_min15_qsig.sh` runs `_min15_qsig`
  (`--qvalue 0.05 --max-sources 100`) for all three, with controls. Perturbations with < 3
  significant features drop out, so compare signature types on shared drugs / perturbations.
  `_qsig` samples the significant features at random and lost the signal (report 05, 5.2–5.3).
  **`_qsigtop` (ready, 2026-11-17):** q < 0.05, then the top 100 of those by |log2FC|
  (`--qvalue-top-n`, new). Reason: of the top 100 by |log2FC|, a median of 1 is significant for MetChange
  and 84 for expr (report 05, 5.4). `bash submit_train_w_shuff_noGlobals_lr1e-3_min15_qsigtop.sh` (also
  MetChange `_qsigtop_sdscale`: plus `--signature-scale sd`).
  **`_qsigtop` done (report 05, section 7): no cost.** expr 10.5, TF 7.8, MetChange 10.3 points (folds 1–4;
  **fold 0 of the MetChange real run missing — rerun it**;
  `..._min15_qsigtop/` has fold1–4 only and no `cv_metrics.json`, while its `_shuff` control and the expr / TF
  `_qsigtop` runs have all 5, then re-run `compare_runs.py` for it); paired vs `_min15` within ±2 points (p ≥ 0.12), +5–8 over `_qsig`.
  `_qsigtop_sdscale` 8.1, not better than `_qsigtop` (−2.0, p = 0.12). MetChange `_qsigtop` covers only
  the 758 drugs with significant metabolites.
- **Feature scale (report 05):** MetChange log2FC is not standardised per metabolite (the MetChange
  paper z-scores against controls), so the current top 100 is mostly a fixed set of high-variance
  metabolites (selection frequency vs metabolite SD: Spearman 0.89; 33 metabolites in > 50% of
  signatures; different drugs share 22% of features vs 2–3% ranked by p). Expr is affected much
  less (5%). Fixes: `_qsig` (scale-free via significance) and **`--signature-scale sd`**
  (implemented: per-feature SD across perturbations, not centred, floored; z-scoring against
  controls rejected, many controls are constant in many features).
  `bash submit_train_w_shuff_mc_noGlobals_lr1e-3_min15_sdscale.sh` (MetChange, seed 0, with control).
  **`_sdscale` done** (fold 3 rerun 2026-11-17, compare re-run on all 5 folds, 1,358 drugs): 6.7 points, vs
  `_min15` −1.6 / +0.4 (p = 0.46 / 0.66), vs `_qsigtop_sdscale` −1.3 (p = 0.62). **SD scaling is neutral: drop
  it from the recipe.**
- Representation (code + GPU): continuous MetChange values or a larger `--top-n` instead of the
  top-100 signs (`--signature-weights`, section 7).

**6.2 Prototype path** (no GPU; this is what "working prototype" needs):
1. **Query path** — nothing yet turns a user's metabolic query into a ranked drug list. With
   `--no-globals` the model can score a query directly (a query has no dose or time), and its
   predicted targets map to drugs that hit them.
2. **Signature-similarity baseline** — rank perturbations by similarity of their MetChange
   signature to the query, no network. If this matches the model, the network is not yet
   earning its place. Needed for the paper either way.
3. **Schuhknecht A549 benchmark** — query with the measured metabolic change of a drug and
   check whether that drug (or its mechanism class) ranks high. This is the only validation on
   the actual task; `spec` only measures whether the signature is read at all. Partly in
   `scripts/mc_schuhk_assoc`.
   **Compartments (to handle when running it):** training signatures are compartment-specific
   (`MAM00003c`; 3,747 of 4,086 MetChange features are exact network nodes, weight ±1 each, the
   ±1/n split in `select_sources` never triggers). Intracellular metabolomics gives one value per
   metabolite, not per compartment, so a Schuhknecht feature must be spread ±1/n over its
   **intracellular** copies only (all compartments except `e`). Today `feature_node_map` maps a
   compartment-free id to *every* copy, including `e`: needs an option to exclude
   extracellular. Train/query mismatch: training sources are single nodes, query sources are
   spread; if the benchmark is weak, consider also training on compartment-collapsed signatures.
   Inputs now in `data/validation/schuhknecht/` (maps from `scripts/schuhknecht_map/`); counts and
   design in report 05, section 6: 473 non-ambiguous intracellular ions → 959 nodes (412 ions /
   801 nodes are MetChange training sources, cluster labels); 646 of 750 matched drugs have targets, 293 were
   CV-evaluated in A549, ~353 never entered training. Use `cpd_info_wSchuhknecht.csv` (all
   compartments), not `metabolite_match_table_non_ambig_match.csv` (`c` only).
   **Script ready (2026-11-17):** `benchmark_schuhknecht.py`, `sbatch launch_benchmark_schuhknecht_A549.sh`
   (report 05, 6.1). Must run on the cluster: locally the rebuilt features do not reproduce the
   training scores (local network / baseline / targets differ from the cluster copies), and the
   script refuses to score in that case. Copy `data/validation/schuhknecht/` to the cluster first.
   Cluster reproduction check passed (r = 1.000000, 8/8). Local MetChange labels are an older
   version (4,086 features vs 5,707 on the cluster): local counts based on them are approximate.
   **Result: chance in both seeds, real and control** (report 05, 6.2 and 8.1).
   Next: does MetChange agree with Schuhknecht at all, on the limma labels:
   `scripts/schuhknecht_mc_assoc/` (`sbatch launch_schuhknecht_mc_assoc_A549.sh`: prep, overlap
   enrichment with MetChange q-values, the R benchmark vs the expression baseline; R ran on the laptop,
   no R on the cluster). **Done (report 05, 8.2): no agreement.** MetChange drug-specific ρ −0.013 vs
   expression producers +0.033 (p = 7e-6); MetChange significant for only 85 of 750 matched drugs, and
   then for most metabolites; overlap enrichment equal to the size-matched expression baseline; sign
   agreement 57% vs 54% mismatched on Schuhknecht-significant metabolites. The benchmark cannot work
   with MetChange-trained models; options: train / calibrate on measured metabolomics, or drop the
   Schuhknecht query as validation of the MetChange path.

**6.0 Recipe as it stands (report 06, 2026-11-18). Not frozen** — kept here to take over from.

| signature | `_min15` (top 100 by \|log2FC\|) | `_qsigtop` (q < 0.05, then top 100) | `_qsig` | `_sdscale` |
|---|---|---|---|---|
| MetChange | 8.3 / 7.7 (seeds 0 / 1; 1,358 drugs) | **10.2** (758 drugs, complete) | 6.1 (758) | 6.7 (1,358) |
| expr | 7.8 / 10.4 | **10.5** (1,358) | 1.6 | — |
| TF | 9.0 / 9.3 | 7.8 (1,346) | 4.2 | — |

(per-drug gap over own shuffled control, `spec_prob_timedose`, points; MetChange `_qsigtop_sdscale` 8.1 / 758.)

Flags, as in the `_qsigtop` launchers (`launch_predict_targets_A549_{mc,tf,expr}_noGlobals_lr1e-3_min15_qsigtop{,_shuff}.sh`,
wrapper `submit_train_w_shuff_noGlobals_lr1e-3_min15_qsigtop.sh`):
`--signature-sign -1 --mode metchange --qvalue 0.05 --qvalue-top-n --top-n 100 --de-counts <expr>/fit_info_nsig.csv
--min-de-genes 200 --score-active-only --lr 0.001 --patience 15 --ffn-epochs 80 --min-epochs 15 --no-vae
--no-globals --n-folds 5`, **no `--signature-scale`**; `networks/gene` for expr / TF, `combined` for MetChange,
corrected scratch labels.

**Settled:** (1) requiring significance is free — `_qsigtop` ≈ `_min15` (±2 points, p ≥ 0.12) and ≫ `_qsig`
(+5–8, p ≤ 5e-5), and every feature read is then significant; expr and TF switch at no cost. (2) **SD scaling
out** (6.7 on full coverage, p ≥ 0.46 vs `_min15`; no gain on top of significance). (3) **Threshold-then-sample
out** (`_qsig`: expr 1.6 points, p = 0.22) — significance must rank, not sample. (4) **MetChange trade-off:**
`_qsigtop` ties `_min15` on its 758 covered drugs but cannot score the other 600, so **`_min15` remains the
primary MetChange run** and `_qsigtop` the significance-backed variant.

**Performance figure (report 06, section 5):** `scripts/predict_targets/plot_performance_panels.py` →
ECDF real vs own control per signature, forest of per-drug gaps with bootstrap CIs, top-k recovery — all on
`_qsigtop` (no other recipe to explain; `--extra sig=run:label` adds comparison rows to the forest when
wanted). Presentation format: 11 x 5.95 in (2:1), runs labelled only MetChange / expression / TF, no panel letters,
transparent background (`--background surface` for opaque). Light ink by default; `--theme dark` writes a
`..._dark` variant with the dark steps of the same hues for a dark slide. Output
`results/predict_targets/A549/performance_figure/performance_panels_<metric>.{png,pdf,tsv}`. The coverage-split
panel was dropped: it only exists for `_min15` (in `_qsigtop` every feature read is significant), so it is the
motivation for the recipe, not a result of it.
**Metric:** headline `rank_probc_heldout` (the deployment task, prior-controlled; MetChange / expr / TF
`_qsigtop` = 10.0 / 8.5 / 8.8 points), `spec_prob_timedose_heldout` beside it as the dose/time-leakage-controlled
check (10.2 / 10.5 / 7.8; the metric of reports 01–05), `rank_prob_heldout` in the supplement (5.8 / 3.4 / 4.0).
**Watch out:** `rank_*_heldout` in `cv_drug_ranks.tsv` are **raw ranks** — divide by `n_candidates_heldout` for a
percentile (only the stratified `spec_*` are percentiles already). Report 05 section 7's `rank_probc` numbers are
in rank units ×100, so their magnitudes are not comparable with the `spec` points; signs and p-values hold.
Re-run when a second `_qsigtop` seed exists.

**Pending / caveats:**
- **MetChange `_qsigtop` fold 0: done** (landed 2026-11-18 00:40, `compare_runs.py` re-run). 10.2 points,
  758 drugs, p < 1e-4; paired vs `_min15` s0 / s1 +0.3 / +1.1 (p = 0.76 / 0.82), vs `_qsig` +6.0 (p = 0.0023),
  vs `_qsigtop_sdscale` +1.8 (p = 0.12). All runs in the table are now complete (5 folds).
- **One seed only** for all `_qsigtop` runs; needed before calling it the recipe (the expr seed spread,
  7.8 vs 10.4, is larger than most differences being compared).
- `_sdscale` is complete (6.7) and can be dropped from future comparisons.
- Open, not selection: representation (continuous MetChange values or a larger `--top-n` instead of 100 signs,
  `--signature-weights`), and validation (6.2b).

**6.0b How to improve performance (2026-11-18).** The current recipe beats its shuffled control by 8-10 points
(median percentile ~0.36-0.40 vs ~0.47-0.50; true drug in the top 10% for ~12% of drugs vs ~5% control): real,
reproducible, not strong. **Selection is exhausted** — `_qsigtop` ≈ `_min15` ≈ `_sdscale` and only `_qsig` is
worse, so *which* 100 features are read barely matters; do not spend more there. In order:

1. **Diagnostic — a predictive similarity baseline scored on the model's metric (no GPU).**
   *Not* the model-free check, which is **already done**: `check_signature_target_similarity.py` →
   `results/signature_similarity/A549/` (expr) and `.../A549/metchange/` ask whether drugs sharing a target
   have similar signatures, and the answer is clearly yes:

   | | expr (3,970 perts, 20,210 genes) | MetChange (2,022 perts, 5,707 features) |
   |---|---|---|
   | same drug, other dose/time (reproducibility ceiling) | AUROC 0.848 | 0.736 |
   | different drugs sharing a target (same time point) | 0.672 (0.693) | 0.561 (0.568) |
   | on drug consensus signatures | 0.590 | 0.529 |
   | 10-NN: neighbours of other drugs sharing a target | 26.6% vs 2.7% expected (**9.9x**) | 15.1% vs 2.8% (**5.5x**) |

   (permuted-label nulls all ~0.50 / 1.0x.) So the signatures do carry target information, and expr carries more
   of it than MetChange — consistent with the model results.
   **RESULT (2026-11-18, report 06 section 6.1): the kNN baseline ranks the true drug BETTER than the model on
   `rank_probc`** (median percentile 0.29-0.31 vs 0.38-0.41, controls all ~0.50; the kNN wins for 55-61% of drugs,
   p ≤ 0.012 MetChange, ≤ 7e-14 expr / TF). Full profile ≈ same input (only MetChange +1.7). Model and kNN succeed
   on different drugs (per-drug Spearman 0.18-0.31) and a rank-average of the two gives 14-16 points. **Caveat:**
   with k = 25, 68-78% of candidates score exactly 0, so `rank_prob` and `spec` are not comparable and the gap
   statistics are noisy; **next is a dense baseline** (every training drug votes, soft weight) before this goes in
   any presentation. **Dense version written and launched 2026-11-18** (`--beta 2 5 10`, launcher default; report 06
   6.2): no zeros except token-level orphans (~7-12% of candidates, flagged per query); check the conclusion
   holds across B and on `rank_prob` / `spec` too. If it holds, the network does not earn its place with the current data, and the priorities
   below change (stacking the kNN vote into the model, or rethinking the architecture, before more data).
   **The predictive version was written (2026-11-18) and run the same day:** `scripts/predict_targets/
   baseline_signature_knn.py` + `sbatch launch_baseline_signature_knn_A549.sh` (cluster only: each fold is checked
   against the run's held-out drugs and test perturbations, and the laptop inputs differ — older labels, no `gene`
   network). For each held-out perturbation, the k = 25 most similar **training** drugs vote their targets, and the
   protein scores go through the model's own `DrugRanker` / `signature_specificity` / `centered_ranks` on the
   run's own folds, queries and candidates, against the run's `_shuff` (shuffle replayed). Two inputs:
   - `sources` = exactly the model's input (selected features, signs, on nodes; cosine) → **network vs no
     network on the same input**;
   - `profile` = the full continuous signature (Spearman) → **is there more in the data than the model reads?**
     (directly informs item 3).
   One call per run scores the run and its `_shuff`, both representations and every k from a single read of the
   labels; the launcher is a 3-task array (one per signature type), so wall time is one expr read, ~30-50 min.
   Output is run-shaped, so `compare_runs.py` and the performance figure read it unchanged; the launcher draws
   `performance_figure/vs_baseline/` with the baseline rows under each model row. Smoke-tested end to end on the
   laptop with `--allow-mismatch` (code path only; numbers meaningless there).
   **Reading:** baseline ≈ model on `sources` → the network adds nothing over similarity of the same input;
   `profile` > `sources` → the 100-sign input is the bottleneck (go to item 3); model ≫ both → the network earns
   its place and the effort goes into data (item 4). Not a baseline: ranking a held-out drug by similarity to its
   *own* other doses (that is the same-drug ceiling, AUROC 0.85). Also not this: `rank_baseline_heldout` in the
   run outputs is the protein-prior baseline (`prior_scores`), signature-independent, identical in real and
   shuffled runs (0.508 / 0.508).
2. **Diagnostic — split the 10 points into prior vs signature.** The model learns a protein prior (validation
   AUPRC 0.348 real vs 0.280 shuffled, base rate 0.027; full proteome 0.0116 vs 0.0074). Retrain with
   `--no-static` (per-protein features removed, so signature only) and, the other way round, with the signature
   removed, and see how the gap splits. If signature-only is near chance, the weak link is signature → target
   and **more data will not fix it**; if it holds up, the prior is just adding a constant and the signature path
   is worth investing in.
3. **Representation: stop throwing the magnitudes away.** `metabolic_target_prediction.py:372` is
   `sigs[p] = (f, np.sign(v[keep].to_numpy()))` — the signature the model sees is **100 features, ±1 each**. It
   never sees how much anything changed, nor feature 101. The `--top-n 100` default arrived with the activity
   filter in `5a28a0f` with no recorded reason. Two small changes, the best effort-to-payoff ratio on this list
   (and the natural reading of why selection makes no difference: the bottleneck is that there are only 100
   signs):
   - **`--signature-weights`** (open in section 7, not implemented): continuous weights instead of signs —
     scaled log2FC, or signed −log10 p. Keep `--signature-scale sd` available here: scaling was neutral for
     *selection* but matters once the value itself is the input.
   - **`--top-n` sweep**: 100 → 300 → 1000 → all significant, each with its `_shuff` control.
4. **Biggest data lever: train across cell lines.** Everything so far is A549 alone (~1,358 drugs); LINCS has
   hundreds of lines and a GNN this size is probably data-starved. Pool 3-5 well-sampled lines with cell
   identity as a covariate: a drug's targets are the same everywhere, so the label is shared while the
   signature varies, which is exactly the contrast the model needs. Costs a label run per line
   (`scripts/compute_labels`), so it is the expensive item — but the one most likely to move the number
   substantially. Ties into 6.3 (data ceiling).
5. **Objective, if 2 points at the signature path:** the FFN trains with per-protein cross-entropy while we
   report a ranking metric; a pairwise / listwise ranking loss optimises what is actually measured.

**Already tried, do not repeat:** matched negatives (`--matched-negatives`, report 01) — once the dose/time
inputs were removed and training ran to convergence, they no longer help (report 02); signature scaling
(`_sdscale`, neutral, 6.0); threshold-then-sample selection (`_qsig`, worse).

**6.2b Validation strategy after the Schuhknecht null (report 05, sections 9–9.3).** Inferred labels do not
agree with measured metabolomics of the same drug (8.2). **Key distinction (9.3): pool metabolomics does not
validate the MetChange labels at all.** A MetChange label claims production / consumption *capacity*; the ground
truth for that is flux. An enzyme inhibitor acts post-translationally, so a correct MetChange need not mirror the
substrate-up / product-down pattern in pools — disagreement is not evidence against it, agreement is a bonus.
Schuhknecht / Sarvin tests are **characterisation**, not validation. Schuhknecht training is deprioritised (one
cancer cell line, ~1,500 drugs vs thousands of LINCS perturbations in hundreds of cell lines; pools cannot
separate more production from less consumption). The expr model is not restricted to metabolic genes: the expr
launchers use all DE genes on `networks/gene` (no GEM layer), so it is a generic transcriptome → target model and
the stoichiometric reading is MetChange's specific contribution.

The chain and what tests each link (report 05, 9.3): (a) signature → target = **done** (CV target recovery,
8–10 points over control); (b) transcriptome → MetChange = the LINCS genetic-perturbation test below; (c) enzyme
transcript → enzyme protein = drug proteomics atlas, overlap not checked; (d) capacity → flux = the claim itself,
prospective tracing only.

**To do, in order:**
1. **LINCS genetic perturbations of metabolic enzymes (link b, the one real label test).** `trt_sh` / `trt_xpr`
   (knockdown / knockout) and `trt_oe` (overexpression) signatures of GEM enzymes, in the cell lines with the
   most coverage. Two readouts: (i) **label-level** — does knockdown of an enzyme read as less production of its
   products (and overexpression as more), including 1–2 reactions away from the perturbed gene? (ii)
   **model-level** — query the trained model with those MetChange signatures: does the perturbed enzyme rank near
   the top? Exact ground truth, no polypharmacology, no new data modality. Partly circular (MetChange reads the
   perturbed transcript itself), so report the 1–2-reaction-away version separately.
   *First step (no GPU):* count how many GEM enzymes have `trt_sh` / `trt_xpr` / `trt_oe` signatures per cell
   line, and whether A549 or another line has the best coverage. Then compute MetChange labels for them
   (`scripts/compute_labels`, same recipe as the drug labels) on the cluster.
2. **MetChange around the target reaction for known metabolic inhibitors** (characterisation, cluster: needs the
   **corrected** labels, the laptop copy is the old uncorrected one). LINCS drugs whose target is a GEM enzyme
   (MTX / DHFR, MPA / IMPDH, leflunomide-teriflunomide / DHODH, statins / HMGCR, pemetrexed / TYMS, …): restricted
   to the target reaction's substrates and products, is MetChange specific to the drug vs other drugs, and in
   which direction (block vs compensation)? Expression read the same way as comparison. Tests the compensation
   hypothesis for the negative slope in 8.2.
3. **Local Schuhknecht test** (characterisation, laptop, runs on the corrected matched matrices in
   `results/mc_schuhk_assoc/limma_rank_*`): drugs whose annotated target **is a GEM enzyme** (see the warning
   below: *not* the `validated_metabolic` column), only the target reaction's substrates / products. Does
   Schuhknecht show substrate up / product down at 48 h, and what does MetChange say at the same metabolites?
   Caveat: few drugs will have an annotated, measured substrate or product.
4. **Knowledge-based enzyme benchmark** (link a, interpretable): expected substrate-up / product-down query
   around the enzymes of drugs with simple documented mechanisms (DHFR, DHODH, IMPDH, HMGCR, TYMS, PAICS, UMPS,
   complex I) — do the known inhibitors rank at the top? Self-built from the GEM, so partly circular with the
   network; cheap and easy to read.
5. **Inborn errors of metabolism** (link a, external ground truth): deficiency → biomarker metabolites
   (Shlomi et al., Mol Syst Biol 2009; IEMbase, OMIM / HMDB). Query = biomarker profile, answer = the deficient
   enzyme. Caveat: plasma / urine pools.
6. **Prospective**: isotope tracing / exometabolomics for a few top predictions (link d, the only direct test of
   the claim).
7. **Seahorse (OCR / ECAR) as a query** (design in report 06 section 7, nothing run). Real flux readout, so it
   tests the capacity claim that pools cannot; coarse (whole-cell respiration / glycolysis only). Query: not
   "every TCA + ETC metabolite up" (at higher flux intermediates often stay flat or fall) but the **endpoints**
   (O2 / NADH / acetyl-CoA consumed, CO2 / ATP produced; glucose consumed, lactate produced), and the two rates
   **together** — OXPHOS up (OCR up, ECAR down), glycolytic shift (the reverse), both down (toxicity, a negative
   control). Two benchmarks: **(b) enrichment of each drug's significant MetChange metabolites in glycolysis /
   TCA-ETC vs its measured OCR / ECAR direction** (one value per drug, well powered, reuses
   `scripts/schuhknecht_mc_assoc/`, tests the labels) **first**, then **(a) query -> drug rank** (the tool's task,
   but only 2-3 queries, so it is a rank-enrichment test of the known positives against a permutation null).
   **Binding constraint: data.** Public Seahorse at screen scale is scarce and usually only in summary tables; for
   (a) curated known ETC / OXPHOS / glycolysis inhibitors and mitochondrial-toxicity screens (Tox21,
   glucose-galactose) supply the positives without raw rates. Check what exists for a LINCS cell line first.

**0. (Do first, cheap) Protein-level target rank in the Schuhknecht benchmark.** `scores.npz` only holds the
drug x drug matrices, but `benchmark_schuhknecht.py` already computes the ~924 candidate-protein scores
internally. **To do: save them and report, per drug, the rank / percentile of its true target protein(s)** among
the candidates (and the mean over the 5 folds for never-trained drugs), with the same permutation null and the
same group split. Small code change + one rerun (~2 min per run on the GPU, runs are on scratch). Why it matters:
the model's output *is* protein scores; ranking drugs is a downstream aggregation, and "lands near the right
enzyme" is both the real task and far better powered than drug retrieval. Report it for all 646 drugs and
separately for the GEM-enzyme-target subset. If the model does not get near the right enzyme here, no new dataset
will rescue it.

**Warning: `validated_metabolic` is LLM-generated and is not "targets a metabolic enzyme".** It flags drugs
*known to have metabolic effects*, which may be indirect (signalling, transport, toxicity), not drugs that
directly inhibit an enzyme. Do not use it as the metabolic-target subset in any of the tests above. **Derive the
subset properly instead:** intersect each drug's annotated targets (the target table the training / benchmark
scripts already use: `--targets`, gene symbols) with the **GEM enzyme genes** (`data/networks/combined`
`gpr_ensembl` / the GEM gene list used by `check_metabolic_coverage.py --gem-genes`), and keep drugs whose target
is a GEM enzyme *and* whose mechanism is inhibition (ChEMBL `mechanism_of_action` / DrugBank actions where
available). Report `validated_metabolic` only as a secondary, clearly-labelled annotation. The same correction
applies to report 05 sections 6, 6.2, 9.1 and 9.2, where the 70 `validated_metabolic` drugs were used as the
metabolic subset, and to the seed-1 `validated_metabolic` spec 0.418 (p = 0.085), which is therefore not even
about enzyme-targeting drugs.

**Sarvin et al. (Nat Commun 2026, report 05, 9.2).** 76 known drugs, 30 metabolic enzymes, HCT116 (a LINCS cell
line), 6 doses × 6 times (10 min–24 h) — the best-matched public dataset (Schuhknecht is A549 at 48 h only).
**But figshare 29268887 holds raw data only:** `Data_mzxml.zip` (8.9 GB mzXML) + `Sample description.xlsx`
(downloaded to `data/validation/sarvin/`: 76 drugs by catalogue number / name, control coding). No processed
intensity matrix, no annotation, no drug → target table (that is Supplementary Data 1 of the paper). Full
untargeted processing (XCMS / MZmine + HMDB annotation + control normalisation) is days of work and would not
reproduce their values exactly; a **targeted extraction of the ~30 target enzymes' substrates / products** from
the mzXML is the lighter option, and enough for test 2–3 style comparisons. Supplementary Data 1 is needed either
way.
**Usable as a Schuhknecht-style benchmark? Only weakly, and not first.** Crude name match of the 76 drugs:
33 in LINCS, **22 evaluated in our A549 CV** (methotrexate, brequinar, leflunomide, mycophenolate mofetil,
pemetrexed, merimepodib, 3 statins, raltitrexed, pralatrexate, TEPP-46, AGI-6780, cytarabine, gemcitabine,
capecitabine, clofarabine, metformin, lonidamine, fomepizole, disulfiram, indoximod). SMILES matching should add
5-FU, mycophenolic acid, 2-deoxy-D-glucose, telaglenastat and the other statins: ~45 in LINCS, ~30 scorable.
Drug-ranking with n ~ 30 has far too little power (Schuhknecht had 646 and was flat), and it is HCT116, so it
needs either cross-cell-line transfer or HCT116 labels + a second training run. **Its value is the enzyme ground
truth** (30 known target enzymes), i.e. the protein-level metric of item 0, which works at n ~ 30 because each
drug is one trial against ~924 candidates. Processing only the useful conditions (8 h / 24 h at ~10 uM, ~300
mzXML files instead of 2,736) keeps it bounded. **Order: do item 0 on Schuhknecht first; process Sarvin only if
that shows signal.** Either way it is a benchmark of the tool, not validation of the labels (9.3): a positive
result counts, a null is uninformative.
Their model predicts the target *pathway* (373 candidates: top 1 29%, top 3 47%, AUC 0.87) from a
time × dose metabolomics matrix, and the enzyme was pinned down by rescue experiments — complementary to ours
(proteins from transcriptomes, free for thousands of LINCS compounds).

R on the cluster: conda env `r-4.3.1` (`launch_schuhknecht_mc_assoc_A549.sh` activates it for the R step).

**6.2c External validation datasets (report 07, 2026-11-19).**

| dataset | readout | cell | mapped to LINCS | status |
|---|---|---|---|---|
| Schuhknecht | metabolomics (pools) | A549 | 750 | done, null (report 05, 8.2) |
| **Attene-Ramos** | **mito membrane potential** | **HepG2** | **890** | **next: section 6 of report 07** |
| Hallinger | Seahorse OCR (flux) | HepG2 | **115 (29 active / 86 inactive)** | **all 1,042 tier-1 chemicals recovered, report 07 sections 8-9** |
| Sarvin | metabolomics, dose x time | HCT116 | ~22-30 | raw mzXML only; enzyme ground truth |
| Lackner | isotope labelling | breast | 1 drug class | case study |
| LINCS shRNA / CRISPR / OE | genetic perturbation | many | — | still the one test of the labels (6.2b item 1) |

> **2026-11-19 (later) — Hallinger solved, two mapping bugs fixed (report 07, sections 8-9).** The paper's
> Supplemental Table 2 lists **all 1,042 tier-1 chemicals with activity calls** (243 active / 799 inactive) and
> the 243 actives' **mechanism and inhibited ETC complex**, so the EPAPLT plate maps are not needed; parsed by
> `scripts/benchmarking/parse_hallinger_supp.py`. **Hallinger -> LINCS HepG2 is 115 compounds, 29 active / 86
> inactive** (19 ETC inhibitors, 3 uncouplers, 4 ATP synthase; rotenone = complex I), up from 42 actives-only.
> **`map_smiles.R` had two bugs** (9.2): `.I` under a filter gave subset positions used as table rows, so right
> structures carried wrong ids; and explicit `[H]` in DSSTox SMILES broke canonicalisation, losing true matches
> (rotenone, estradiol). Both fixed. **Attene-Ramos was re-mapped: 890 compounds, and 0 of the old 1,509 pairs
> survive** - anything built on the old `mapped_drugs.csv` must be redone (old file kept as
> `mapped_drugs_BUGGY_2026-11-19.csv.bak`). **Schuhknecht checked and is clean** (report 07, section 10): that
> mapping is separate code which matches on values, so the report 05 null stands - 998 of its 1,040 pairs
> reproduce and the only drug lost is Gramicidin, a `NA`-matched-`NA` artefact. Re-mapping it would *gain*
> 212 drugs (1,251 total), since explicit `[H]` (228 of 1,520 Schuhknecht SMILES) and salts blocked them.

> **2026-11-19 (later) — first non-null external association (report 07, section 11).** Label-level test of the
> HepG2 MetChange labels against measured mitochondrial phenotypes, `scripts/mc_phenotype_assoc/` (40 Human-GEM
> metabolite sets, per-drug `z` / rank statistic, Mann-Whitney + bootstrap + permutation + **size-matched
> random-set control**, no imposed sign, model activity filter deliberately NOT applied).
> **Attene-Ramos: MMP-down drugs have lower MetChange TCA capacity** (`tca mito_nohub` AUROC 0.364, q = 7e-6,
> random-set p = 0.010; non-hub variants strongest, glycolysis weak at 0.43, **oxphos / ETC flat**), and
> **MMP-up drugs have higher energy-currency capacity** (0.663, q = 2e-4). The 24 cytotoxic drugs are flat, so
> it is not a cytotoxicity readout. **Hallinger is null** (29/86; 19 ETCi vs 86 — same direction, CIs cross 0.5,
> uninformative at that n). **Confound, carry it with the result:** MMP-active drugs are 3.4x more likely to have
> a responsive MetChange signature (31% vs 9%), and within the *controls alone* responsive drugs already show the
> low-TCA pattern (0.355, p = 0.0015). Stratified, TCA holds at ~0.38 in both strata (non-responsive p = 2e-4) and
> MMP-up/currency is clean (0.698 non-responsive, confounder absent). **Decisive next control: the same test on the
> HepG2 expr labels with GPR gene sets** - needs `log2FC.csv` for HEPG2/expr (the laptop copy has only
> `qvals.csv` / `basal.csv`), then `launch_assoc_expr_HEPG2.sh` runs it; sets via
> `build_metabolite_sets.py --emit genes` (Ensembl both sides, no mapping).

> **2026-11-19 (later still) — three modalities (report 07, section 12).** `--mode value|sign|sig`;
> refactor verified (all 80 shared statistics reproduce exactly). **The directional reading confirms TCA and is
> the cleanest evidence**: `sign_dev tca cyto` AUROC 0.415 (q = 0.041, rand = 0.020), non-responsive stratum
> 0.432 (p = 0.03) and **the responsiveness confounder is absent** (0.453, p = 0.30) where it was present for
> the magnitude statistic. `sign_mean` null, so TCA is down *relative to each drug's own mean sign*.
> **The `sig` modality is a responsiveness artefact - do not use it pooled**: `tca mito_nohub` q = 4e-10 pooled
> but AUROC exactly 0.500 (p = 1, all ties) among non-responsive drugs and 0.091 (p = 2e-71) for responsiveness
> within controls alone; restricted to responsive drugs nothing survives (best p = 0.079). Always pair
> `--mode sig` with `--min-nsig 1`. **The apparent cytotoxic signal is the same artefact** - the clean cytotox
> control of section 11 holds only for value/sign. Expr baseline started but only in `sig` (no `log2FC.csv`):
> no gene set beats its random-set control (rand >= 0.25), suggestive but not like-for-like.
> **Launchers:** `scripts/mc_phenotype_assoc/launch_assoc_{metchange,expr}_HEPG2.sh` +
> `launch_assoc_HEPG2_slurm.sh` (2-task array); config in `mc_phenotype_assoc_common.sh`, which auto-detects
> the scratch base on the cluster.

> **2026-11-24 — LABEL GENERATIONS: the 2026-11-12 files ARE the corrected ones (user). The superseded buggy
> generation was dated 2026-11-02.** `results/labels/A549/metchange/rank/log2FC.csv` (4,427,631,022 bytes,
> mtime Nov 12 13:17) is correct and is what the local analyses used, so **the cross-cell signature-conservation
> result STANDS** — a withdrawal issued earlier on 2026-11-24 was unnecessary and has been reversed. Nothing in
> a label file's name or content identifies its generation, so **state the matrix mtime in any write-up** and
> read the `labels` path in a downloaded result's `run_info.json` (a `/scratch/...` value means it was computed
> on the HPC, not from the local copy). `results/labels/A549/compare_corrected_vs_rerun/` exists from that
> comparison. `scripts/mc_phenotype_assoc/signature_conservation.py` (new) carries the provenance warning.
> **VERIFIED (report 10, 6.1): the section-3 A549 null reproduces EXACTLY from the local Nov-12 labels** —
> all 40 shared set x variant combinations identical (median AUROC diff +0.0000, max |diff| 0.0000; tca
> mito_nohub 0.477 p=0.34 both), same 179/803 contrast. So the HPC cross-cell runs used this same generation,
> which means the set-level conservation table is on corrected labels too. (Minor: the local matrix has 9,830
> perturbations / 1,085 drugs vs the HPC run's 10,150 / 1,118 — the extra columns are in OTHER contrasts
> (`MMP_up`, `cytotox`, `incon`), so don't compare those across the two.)
> **CONSERVATION MEASURED (report 10, 6.2-6.4), answering the correct objection that shared MoA + shared core
> metabolism make SOME transfer the EXPECTED BASELINE, not a finding.** Full signatures HEPG2 vs A549, 250
> drugs x 5,707 metabolites: same-drug median r +0.297 vs different-drug +0.031, but **same-vs-different AUROC
> only 0.551 (p=0.0055)** and just 10% of same-drug r above the different-drug 95th pct — **conservation is
> REAL but WEAK; signatures are dominated by cell context, not compound identity.** Set-level, responsive-only:
> energy_currency +0.399, glycolysis +0.474, oxphos +0.337, **tca +0.113 to +0.145 (lowest)**. **Set size does
> NOT explain it: Spearman(members, conservation) = -0.058 (p=0.72).**
> **CONSEQUENCE: section 3 is no longer the anomaly — weak conservation means A549 labels carrying little
> information about a HepG2 TCA phenotype needs NO appeal to HepG2 being special.** Reframes 9.3 ("refutation
> of the strong reading") and 10.5 ("the damaging direction") as the expected baseline being observed.
> **WITHDRAWN OVERSTATEMENT: "TCA is not conserved at all" and the clean "what transfers doesn't predict"
> complementarity.** On the best-powered pair TCA conservation is POSITIVE (+0.11..+0.32), merely lowest; the
> negatives that prompted it (HCC515 -0.371, A375 -0.154) rest on 45-61 responsive drugs. Directional only.
> Not addressed: only HEPG2-vs-X pairs, responsive n 31-68, and `z` is already within-perturbation standardised.

> **2026-11-25 — REPRODUCIBILITY GAP FOUND AND CLOSED (report 10, 16).** Asked whether the reported results had
> a launcher; the answer was **no**. The SCORING steps always did (`launch_prism_cell/slurm.sh`,
> `prepare_assoc_inputs.sh`, `run_naia_assoc.sh`), but **every number in report 10 sections 8-15 came from an
> ad-hoc command line** — reproducible only by remembering which of `prism_differential.py`'s **33 flags** were
> set. Nominal reproducibility: script in git, result not, and a changed default would alter a published number
> silently. **Section 7 identified this exact failure for the cross-cell work, fixed it with
> `cross_cell_compare.py`, and it was then repeated a dozen times**; section 13's subsystem loop existed only in
> a session temp dir.
> **ALSO ADDED `reproduce_report09.sh`** — the HepG2 MMP association itself was always scripted
> (`launch_assoc_HEPG2_slurm.sh` runs BOTH arms with stratified + matched checks), **but report 09's cross-cell
> table (section 3) and the MetChange-vs-expression HEAD-TO-HEAD (section 4) were NOT**, because both apply the
> common-drug restriction the per-cell launchers skip. **Section 4 holds the project's most load-bearing number
> for the thesis and had no launcher.** Now verified to return HEPG2 MetChange 0.350 pooled / 0.365 matched and
> expression 0.418 / 0.500 matched, exactly as reported. Sections 9/9.6 use `run_naia_assoc.sh` (always existed).
> **NEW CONTEXT the launcher surfaced: the RESPONSIVE SHARES differ hugely between label types — MMP-down in
> HEPG2 is 28%/7% (MetChange) vs 62%/50% (expression). Expression calls far more drugs responsive with a much
> smaller case-control gap, which is why the matched designs behave so differently — bear it in mind before
> reading the 0.500 as a clean null.**
> **NOW: `reproduce_report10.sh` runs every section 8-15 analysis with exact flags** (`bash ... reproduce_report10.sh`
> for all, or `... 12 13` for selected); `prism_cell_coverage.py` reproduces section 15 and is reusable for any
> cell-line question; `build_target_relevance.py` for section 14. Only prerequisite is the per-cell scores from
> the cluster. **Verified: section 12 returns +0.291 (10/10), expression -0.008 (4/10), MC|expr +0.302 (10/10),
> retention 104%, Wilcoxon at the n=10 floor — exactly as reported.**
> **NOT FIXED, and the right next housekeeping: `prism_differential.py` is 811 lines / 33 flags and the headline
> depends on a 6-flag combination. A launcher PINS that but does not make it LEGIBLE. A load-bearing claim wants
> a ~40-line script with the specification HARDCODED and no arguments**, with the flexible script kept for
> exploration. Also: of ~2,600 lines written for report 10, ~600 serve no live claim (`build_xu_phenotype.py`
> underpowered test, `signature_conservation.py` one use, `build_target_relevance.py` withdrawn as a selector
> within the hour). **Volume is not progress; the remedy is a decision criterion, not more code.**

> **2026-11-25 — WHICH CELL LINES ARE WORTH COMPUTING: TWO, NOT 93 (report 10, 15). CORRECTS 12.3.** Measured
> from the new `data/LINCS/lincs_unique_perts_bycell.tsv` (227 cells, 201,295 rows) against PRISM, BEFORE
> committing to label runs.
> **89 cells are in both. Coverage of the 31 nucleotide-targeting drugs: 29-30 for the five cores (labels
> done); 24 for YAPC; 18 for MDAMB231; 11 for U2OS; <=6 for the remaining 81.**
> | cell | LINCS panel | PRISM ovl | nucleotide | labels |
> | PC3/A375/MCF7/A549/HT29 | 12,156-13,892 | 861-927 | 29-30 | yes |
> | **YAPC** | 1,716 | 605 | **24/31** | no |
> | MDAMB231 | 1,357 | 480 | 18/31 | no |
> | U2OS | 20,836 | 438 | 11/31 | no |
> | HEC108/CAL29/BEN/... | ~723 | ~220 | 5/31 | no |
> **CAUSE: panel size and COMPOSITION, not cell availability.** Non-core LINCS panels are modally 120 perts
> (13 cells) or 723-724 (17 cells) — a small fixed plate — vs 12,000+ for the cores. **U2OS is the instructive
> case: 20,836 perturbations and still only 11 of 31, because its panel is a DIFFERENT chemical space.**
> **WHAT IT BUYS (fixes the saturated Wilcoxon of 11.6, floor = 1/2^pairs): +YAPC -> 6 cells, 15 pairs, floor
> 3.1e-05 — BREAKS THE SATURATION on its own; +MDAMB231 -> 21 pairs, 4.8e-07; +U2OS -> 28 pairs, 3.7e-09.**
> **RECOMMENDATION: compute YAPC, then MDAMB231 if the first is informative. DO NOT compute the other 86** — a
> cell with 5 of 31 drugs gives a correlation on five points. **YAPC is also one of the seven Panorama GCP cell
> lines (section 2), so one label run serves two analyses.**
> **FIFTH INSTANCE OF THE SAME ERROR (15.3): 12.3's "93 available" came from `cellinfo` membership without
> checking per-cell drug coverage** — as with Hallinger, naia (683 -> 10), xu (53 -> 11/12) and the PRISM
> subsystem comparators (1-4 drugs each). **General form of the rule: COUNT THE UNITS THAT CAN ACTUALLY CARRY
> THE PREDICTION, for the specific comparison proposed, BEFORE quoting an opportunity.**
> `lincs_unique_perts_bycell.tsv` makes this a one-command check — consult it first for any cell-line question.

> **2026-11-25 — PRINCIPLED TARGET FILTER IMPLEMENTED, deliberately NOT run-and-quoted on PRISM (report 10,
> 14).** `scripts/mc_phenotype_assoc/build_target_relevance.py` (new) emits `gene, set, hops` — shortest
> PRODUCING route from a gene to a set member, currency metabolites excluded. Three new flags on
> `prism_differential.py`: `--relevance/--max-hops`, `--min-target-sources`, `--exclude-target-classes`.
> **IT REDERIVES THE HAND-MADE LIST: for `dntp_endpoints` it rules out 1,509 of 2,522 GEM genes (60%) as NEVER
> REACHED, and the never-reached set includes PDE, PI3K and EZH2** — exactly the contaminations found by
> inspection in 10.1, 11.3 and 13.4.
> **TWO LIMITS, MEASURED: (a) hop distance does NOT separate biosynthetic from incidental routes** — to a dNTP,
> TYMS/TK1 sit at 4 alongside EHMT2, DHFR/IMPDH1 at 5 alongside CYP2C19, DHODH at 6 alongside MAOA (gene counts
> 19/43/117/402/837/1415 at hops 1-6). **(b) it only discriminates for SMALL endpoint sets**: `dntp_endpoints`
> gives 8/14/42/116 drugs at H<=2/3/4/6, but `nucleotide` (192 metabolites) and glycolysis/tca/oxphos/ppp give
> ~235 at ANY depth, because near-currency members put every gene a few hops away.
> **TWO FAILURE MODES DISTINGUISHED, and a CORRECTION to 11.3: GCK, HK1 and EEF2 should NOT have been in the
> exclusion list — they are weak ANNOTATIONS (all n_sources=2 from CLUE/RepurposingHub), not misplaced enzymes;
> a real hexokinase inhibitor WOULD be a legitimate glycolysis drug.** `n_sources >= 3` removes all three while
> every genuine antimetabolite annotation survives (3-5, incl. ChEMBL + DrugCentral). 11.3's result is
> unaffected; the reasoning was wrong.
> **CLASS LIST (report 10, 14.3), principle = present in the GEM for catabolic / cofactor-consuming /
> signalling / xenobiotic reasons rather than as a biosynthetic route:** PDE + lipid kinases (already automatic);
> then CYP/UGT/SULT/NAT/FMO, SAM methyltransferases (EHMT/EZH/DNMT/SETD/KMT/PRMT/COMT/NNMT), MAO/ALDH/ADH/AOX1/DBH,
> SIRT/PARP/CD38, TXNRD/GSR/GPX/GST/NQO1/SOD/CAT, PTGS/ALOX/PTGES/HPGDS/TBXAS1, NOS, CA, HSD/SRD5A/STS.
> **DELIBERATELY ABSENT: NAMPT (NAD synthesis IS biosynthetic) and TRANSPORTERS (SLC7A11 cystine, SLC1A5
> glutamine genuinely limit precursor supply — flag, don't blanket-exclude).**
> **COST: the full stack leaves 25 drugs overall and 12 for `dntp_endpoints` (vs 28 hand-cleaned) — it HALVES
> the power on the one set where the effect was visible, because hop<=4 drops the DHFR/IMPDH/DHODH drugs (hops
> 5-6). And A BETTER FILTER CANNOT CREATE DRUGS: PRISM still has 1 PPP drug, 3 TCA, 4 glycolysis, 4 OXPHOS.**
> **THEREFORE NOT RUN ON PRISM — it is specified now so it can be applied BEFORE looking at results on a library
> where the drugs exist (Tox21/ToxCast). NEXT STEP: check usable n under this stack on those libraries FIRST.**
> **CORRECTED 2026-11-25 (report 10, 14.5) — USE THE CLASS EXCLUSION ALONE. The producer filter is a
> DIAGNOSTIC, not a selector.** The objection that halving the power on drugs that genuinely do affect dNTP
> production means the filter is doing harm is correct. **DHFR, IMPDH and DHODH ARE on producing routes to
> dNTPs** (DHFR makes the THF that TYMS needs; IMPDH runs IMP->XMP->GMP->dGTP; DHODH feeds orotate->UMP->
> dCTP/dTTP) — they sit at hops 5-6 because REAL BIOSYNTHETIC ROUTES ARE LONG. Hop distance measures path
> length, which is not relevance.
> **DECISIVE COMPARISON (with PDE + lipid kinases moved INTO the class list): subsystem + CLASS = 33 drugs
> retaining ALL 28 of 11.3. Adding producer-reachability = 31 but LOSES ABT-702 (a real adenosine kinase
> inhibitor). `n_sources>=3` = 25 but LOSES carmofur, cladribine, AVN-944 (real antimetabolites).**
> **REVISED RULE, simpler than 14.3-14.4: (1) target in a matching GEM subsystem, effect=inh; (2) exclude the
> 14.3 classes WITH `^PDE`, `^PIK3|^PIP[45]K|^SPHK|^DGK` and `^XDH` ADDED to the list (XDH is purine catabolism
> — how curcumin and diphenyleneiodonium get in); (3) nothing else — NO hop cutoff, NO source threshold in the
> primary rule.** Gives 31 drugs, retains all 28, ~10% residue (3 bad annotations) which
> `--min-target-sources 3` can quantify as a SENSITIVITY CHECK without being baked in.
> **The producer filter's remaining role: it DERIVED the class list** — running it is what established that
> PDE, PI3K and EZH2 are NEVER on a producing route to a dNTP, replacing three rounds of inspection (10.1,
> 11.3, 13.4) with a stoichiometric criterion. Keep `build_target_relevance.py` / `--relevance` for that and
> for vetting a new class before adding it. **Do not use them to choose drugs.**

> **2026-11-25 — SPECIFICITY TESTED PROPERLY: PRISM CANNOT ANSWER IT (report 10, 13).** Ten runs, one per
> subsystem, drugs selected by inhibited GEM target and scored on the SAME subsystem's metabolite set.
> **nucleotide 28 drugs +0.254 (10/10, p=0.00098); pyrimidine 20 +0.306 (10/10); purine 13 +0.152 (8/10);
> amino_acid 48 -0.061 (4/10, p=0.78); folate 4 / glycolysis 4 / oxphos 4 / pyruvate 4 / tca 3 / ppp 1 = NOT
> SCORABLE.**
> **(1) The three positives are ONE finding, not three: purine (13) and pyrimidine (20) are ENTIRELY inside the
> nucleotide 28, and 27 of the 28 fall in one or the other.** Same drug set re-sliced.
> **(2) Six comparators have 1-4 drugs and are contaminated with the SAME drugs:** `ppp` has ONE drug; `tca`'s 3
> are a NOX inhibitor, an ALDH inhibitor and an uncoupler; `oxphos`'s 4 include brequinar + teriflunomide
> (DHODH inhibitors already in the 28 — DHODH is in OXPHOS because it feeds ubiquinone); `folate`'s 4 are the
> antifolates, also in the 28. A repurposing library has dozens of antimetabolites and essentially no approved
> drugs against GAPDH, citrate synthase or complex I.
> **(3) `amino_acid` looked like THE comparator and is NOT one.** 48 drugs, larger than the 28, nearly
> independent of it (3 shared), and null — so not a power artefact. But ~95% contamination: **13 via CYP2C19**
> (imatinib, gefitinib, vemurafenib, axitinib, rifampicin, amoxicillin... drugs METABOLISED by it), **6 via
> MAO** (monoamine catabolism sits in tyrosine/tryptophan), **2 SAM-dependent methyltransferase inhibitors**
> (BIX-01294, E-7438 — in cysteine/methionine because their targets consume SAM), 2 sirtuin, ~20 via
> AKR1C1/AOX1/TYR/NOS/MIF/ALDH/COMT. **Genuinely amino-acid-metabolism: ~1 (nitisinone).** So the null means
> "no amino-acid-directed drugs here", not "amino-acid capacity is not limiting". **AND THE COMPARISON WAS
> ASYMMETRIC IN MY FAVOUR: the nucleotide arm was cleaned 44 -> 28, this one was left uncleaned at 48.**
> **(4) METHODOLOGICAL CONCLUSION THAT GENERALISES: target-based selection on GEM subsystems needs a CURATED
> ENZYME WHITELIST, not a subsystem regex.** Human-GEM subsystems correctly absorb drug metabolism (every CYP),
> monoamine catabolism (MAO), SAM-dependent methylation (every methyltransferase), NAD/cNMP second messengers
> (sirtuins, PDEs) and phosphoinositide signalling (PI3K). **The same failure has now produced a
> plausible-looking drug group THREE times: PI3K under inositol phosphate (10.1), PDE under purine (11.3),
> CYP/MAO/methyltransferase under amino acid (13.4). The subsystem a gene's reactions belong to does not tell
> you whether a drug acting on it perturbs that pathway's flux.**
> **(5) WHERE IT STANDS: the 12.1 result holds as "for the ONE metabolic drug class this library contains,
> MetChange's cell-specific prediction works and the matching gene expression does not" — and CANNOT be
> sharpened to "specifically nucleotide" on this evidence.** Routes out: (a) more cell lines, same drugs (12.3,
> 93 available) — fixes the saturated p-values but not specificity; (b) a library with broader metabolic
> coverage (Tox21/ToxCast), which brings back the chemical-space mismatch that produced the Hallinger and Xu
> nulls — **check usable n for a CLEANED, WHITELISTED target selection there BEFORE designing anything around
> it.**

> **2026-11-24 — EXPRESSION ARM RUN: THE PRE-REGISTERED CONDITIONAL PREDICTION IS CONFIRMED (report 10, 12).**
> All 5 expr arms completed and downloaded (20,209-20,214 genes, 679-760 drugs).
> **IMPLEMENTATION GAP FOUND AND FIXED FIRST: the conditional columns computed the partial over ALL ~680 drugs,
> not within the matched subgroup** — a different question from the one registered in 10.7. Added
> `r_partial_moa` / `r_expr_moa` / `r_partial_other`, computed within the matched group and the complement
> separately. The numbers below are the registered statistic.
> **`dntp_endpoints`, 28 matched drugs, per cell pair: MetChange mean +0.291 (10/10 pairs positive);
> EXPRESSION (`nucleotide` gene set) mean -0.008 (4/10); MetChange|expression +0.302 (10/10). RETENTION 104%**
> (registered prediction: "well above 50%"). Paired Wilcoxon of the partial, matched vs complement: 10/10, p at
> the n=10 floor.
> **STRONGEST RESULT IN THE PROJECT FOR THE THESIS, and the FIRST measurement of INCREMENT rather than a
> side-by-side comparison** (10.6) — on five NON-HepG2 cell lines, with a mechanism-matched drug set.
> **BUT READ WHY RETENTION IS 104%: not that expression carries a strong signal MetChange survives, but that
> the GENE-SET DIFFERENTIAL IS FLAT (-0.008), so conditioning removes nothing because there is nothing to
> remove.** Cleanest possible form of "the transcriptome alone does not see this", and a WEAKER statement than
> "MetChange adds to a strong transcriptional signal". The comparator is fair (`nucleotide` gene set = GPR genes
> of the same subsystems whose metabolites `dntp_endpoints` summarises) but the finding is an ABSENCE on the
> expression side.
> **THE FORKING PATH IS UNTOUCHED BY THIS: the conditional statistic was registered, the 28-DRUG SET WAS NOT**
> (five successive narrowings, 11.4). Correct description: *within a post-hoc-selected drug set, a
> pre-registered conditional test behaves exactly as the thesis requires.* All p-values still saturated at
> n = 10 — quote effect sizes.
> **SCALING — 12.3's "93 cell lines" is WRONG, see report 10 section 15: only TWO are worth computing.**
> (Original claim, retained for context: PRISM covers 568 cell lines, 99 are in `cellinfo_beta.txt`,
> 6 have labels and 93 do not — incl. HCT116, HUH7, JHH5, JHH7, DU145, BT474,
> 22RV1, HCC827, HCC44 + deep endometrial/ovarian/bladder/colon panels. **With 5 cells the test cannot resolve
> below p = 0.000977 (10 pairs); 20 cells gives 190 pairs, 50 gives 1,225. The confirmatory cohort does not need
> a new dataset — it needs more LABEL RUNS on the one in hand.** GATING CHECK is LINCS DEPTH not `cellinfo`
> membership — and that check has now been done, see below.)

> **2026-11-24 — JACCARD BOTH WAYS + A CONTAMINATION; THE SECTION-10 RESULT IS NOW FROZEN, NOT CONFIRMED
> (report 10, 11).**
> **(1) Jaccard as AGGREGATION substantially WEAKENS section 10.** `dntp_endpoints`, nucleotide-targeting:
> ungrouped +0.136 / 80% of pairs / Wilcoxon p = 0.0020 -> **Jaccard 0.5: +0.055 / 60% / p = 0.080**;
> Jaccard 0.7: +0.056 / 60% / p = 0.053. `nucleotide`/`purine`/`pyrimidine` FLIP NEGATIVE. The effect-size
> halving is the estimate moving, not lost power. **Cause: of 44 matched drugs, the largest Jaccard group is 8
> FLUOROPYRIMIDINES (`{DPYD,TYMS}`: capecitabine, carmofur, doxifluridine, floxuridine, ftorafur,
> trifluridine...) — one chemotype family counted 8 times.** 44 drugs -> 26 groups.
> **(2) Jaccard as STRATIFICATION (the per-mechanism version) CANNOT test the hypothesis with this library: NO
> matched group reaches 10 drugs** (8, 4, 4, 3, 3, 2, then 20 singletons). The 9 groups that do are ALL
> non-metabolic (MET 35, ATR/MTOR 30, FGFR 19, HDAC 16, AR/ESR1 16, MEK 13, HSPA1A 13, tubulin 10, Aurora 10)
> and behave as a **negative-control panel**: median of group medians **+0.016**, individual values -0.217..+0.236,
> sign fractions 20-70% — noise. Weak reassurance that the design does not manufacture associations where the
> hypothesis does not apply.
> **(3) CONTAMINATION: 16 of the 44 "nucleotide-targeting" drugs are NOT antimetabolites — 9 are
> PHOSPHODIESTERASE inhibitors** (papaverine, anagrelide, enoximone, cilostamide, zardaverine, YM-976,
> trequinsin, avanafil, arctigenin), because cAMP/cGMP hydrolysis is filed under purine metabolism in
> Human-GEM; plus GCK, HK1, EEF2, TXNRD1/2, XDH, aspirin. **GEM subsystem membership is BOOKKEEPING, NOT
> MECHANISM** — same failure as PI3K-in-inositol-phosphate, one level subtler. `--exclude-target-genes` added.
> Clean 28 drugs: **ungrouped +0.291, 10/10 pairs positive; Jaccard 0.7 +0.240, 90%; Jaccard 0.5 +0.187, 70%
> (p = 0.019).**
> **EXACT VALUES (`cross_cell/prism_clean28_*.tsv`): ungrouped +0.291, 10/10 positive, 10/10 diffs favouring,
> p = 0.000977; Jaccard 0.7 +0.240, 9/10, 9/10, p = 0.001953; Jaccard 0.5 +0.187, 7/10, 8/10, p = 0.018555.**
> **ONLY the ungrouped p is SATURATED: 0.000977 = 1/2^10, the floor at n = 10 pairs**, so it says only "all 10
> pairs favoured matched" and carries no effect-size information (would read 1e-06 with 20 pairs).
> **CORRECTION: an earlier claim that the Jaccard-0.7 row ALSO sat on the floor and was "indistinguishable by
> p-value" was WRONG — its p is 0.001953, exactly 2x the floor (9/10 diffs, not 10/10). It came from a value
> typed into a demonstration instead of read from the output, and the conclusion was then written into the
> report.** Also corrected: the illustration that `frac r>0` and the Wilcoxon can disagree rested on that same
> bad number — in all three rows they track closely (10/10 & 10/10, 9/10 & 9/10, 7/10 & 8/10). They are still
> DIFFERENT statistics (`r_moa > 0` vs `r_moa > r_other`) but the disagreement was not observed here. **The clean 28 SURVIVES pseudo-replication control where the
> contaminated 44 did NOT (p = 0.080 at Jaccard 0.5) — the one piece of internal evidence here that is not just
> a consequence of narrowing.**
> `amino_acid` stays negative ungrouped (-0.084) but moves to +0.032 p = 0.023 grouped — **incoherent, a direct
> warning that the grouped n = 14 estimates are too noisy to read precisely.**
> **(3b) WHAT THE TEST IS, LEVEL BY LEVEL (report 10, 11.6).** 28 drugs = 821 in the phenotype table -> 747
> with any annotated target -> 44 with an INHIBITED target in a `nucleotide|purine|pyrimidine|folate` GEM
> subsystem -> 28 after the exclusion regex. **Level 1:** per cell pair, Spearman across drugs — one rho on the
> 28 matched, one on the ~650 complement. Matched ranges 0.094-0.502; complement -0.097..+0.030. **Level 2:**
> Wilcoxon SIGNED-RANK on the 10 differences, one-sided. **So the test's n is 10 CELL PAIRS — not 28 drugs, not
> 820 combinations.** W = 55 = the maximum at n=10, hence p = 1/2^10.
> **DEPENDENCY not stated earlier: the 10 pairs come from 5 CELLS (each cell in 4 pairs) and the SAME 28 DRUGS
> appear in every pair — so the 10 differences are non-independent in both directions while the Wilcoxon assumes
> independence. Its p is optimistic on top of being saturated.**
> **LEAVE-ONE-CELL-OUT (6 pairs each): mean difference +0.288..+0.362, every p = 0.0156 = ITS floor (6/6
> positive). Per-cell mean rho_matched A375 +0.273, A549 +0.271, HT29 +0.285, MCF7 +0.378, PC3 +0.246 — NO
> SINGLE CELL DRIVES IT.**
> **EVERY p-value in this analysis is saturated; they report sign consistency only. QUOTE THE EFFECT SIZE AND
> ITS STABILITY (matched rho ~0.25-0.38 per cell vs complement ~ -0.05), NOT A P-VALUE. This is also why the
> confirmatory cohort needs more CELL LINES, not just more drugs: with 5 cells the aggregate test cannot
> resolve below 0.000977 however strong the effect.**
> **(4) THIS SEQUENCE CAN NO LONGER ESTABLISH ANYTHING. The effect grew at EVERY refinement: all drugs -0.036 ->
> MoA (incoherent) -> any energy+biomass GEM target +0.099 -> nucleotide target +0.136 -> contaminants removed
> +0.291. That is the textbook forking-path signature, and each step being individually defensible does not
> change it** — every narrowing was chosen after seeing the previous number and every one moved toward a larger
> effect. **So nothing in 10.1-11.3 is confirmatory evidence.**
> **WHAT IT DID PRODUCE: a fully specified hypothesis with no free parameters, FROZEN in report 10, 11.4** —
> the 28-drug list (nucleotide/purine/pyrimidine/folate GEM target, `inh`, minus
> `^PDE|^TXNRD|^XDH|^GCK$|^HK1$|^EEF2$`), `dntp_endpoints` variant `all`, per-pair Spearman aggregated by
> paired Wilcoxon vs complement, Jaccard 0.5 group medians as the conservative form, positive direction, with
> `amino_acid` predicted null/negative as the internal contrast. **The next number quoted for this hypothesis
> must come from data none of this touched: PRISM primary, CTRP (usable n unchecked), or GDSC — run once.**
> **(5) TARGET-PROXIMAL SETS: not runnable locally** — per_drug_scores holds the 87 predefined set x variants,
> while a target-proximal set is PER DRUG, so it needs a scoring pass against the label matrices on the
> cluster. **Now the most valuable remaining analysis:** with a shared `dntp_endpoints` a DHFR, a DHODH and a
> TYMS inhibitor are all judged on dNTP pools; per-drug neighbourhoods would score each on its own
> (folate/one-carbon, orotate, dTMP/dUMP) AND **collapse the family from 820 combinations to ~10**. Needs a
> `--drug-sets FILE` mode in `mc_phenotype_assoc.py` (columns `drug, node_id`), built by a 1-2 hop walk from
> each target gene in `nodes.tsv`.

> **2026-11-24 — THE BIGGEST GAP: BOTH "MetChange BEATS EXPRESSION" INSTANCES USE **HEPG2 LABELS**, and the
> PRISM runs have NO EXPRESSION ARM (report 10, 10.5 / report 09, 9.6 note).**
> **There are NO SHSY5Y labels and none are possible** (224 LINCS perturbations in SHSY5Y, only 10 naia
> compounds among them); `results/labels/` has no SHSY5Y dir. Both `naia_atp` and `naia_atp_expr` have
> `labels = results/labels/HEPG2/{metchange/rank,expr}` — **SH-SY5Y is the cell the PHENOTYPE was measured in,
> nothing more.** Report 09, 9.6's column header "(ATP, SH-SY5Y)" invited the opposite reading.
> So the two instances (report 09 section 4 MMP: MetChange 0.365 vs expression 0.500; 9.6 ATP: 0.173 vs 0.329)
> vary phenotype, phenotype cell line, compound library and confounder control — **but NOT the label cell line.
> Anything peculiar to HEPG2's label generation (reference level, limma fit, its LOWEST responsive share of the
> eight lines) is common to both.** Calling them simply "independent" overstated it.
> **FIX, and it is the cheapest test of the project's central thesis: add an expression arm to
> `launch_prism_cell.sh`** (follows the `launch_assoc_expr.sh` pattern; `gene_sets.tsv` already has the
> nucleotide/amino-acid GENE sets from the 9.5 rebuild). Every Attene/Hallinger run has an `_expr` counterpart;
> `prism_viability` has none. **It would be the FIRST MetChange-vs-expression comparison on NON-HEPG2 labels —
> five cell lines at once (A375, A549, HT29, MCF7, PC3).**
> **THE INCREMENT QUESTION HAS NEVER BEEN ASKED IN THIS PROJECT (report 10, 10.6 note).** Every result so far
> uses an UNCONDITIONAL statistic: sections 8/9/10 are plain Spearman of the differences, and the two existing
> MetChange-vs-expression comparisons (report 09 sections 4 and 9.6) are **SIDE-BY-SIDE** — greedy covariate
> matching then Mann-Whitney, run separately on each label type. **Partial Spearman has NEVER been run on real
> data** (added 2026-11-24, tested only on a synthetic arm). MetChange surviving a control that expression fails
> is suggestive of increment but does NOT measure it — two statistics can both be non-null while sharing all
> their information. **So the present evidence for the central thesis is weaker IN KIND, not just in quantity,
> than "MetChange 0.365 vs expression 0.500" suggests.**
> **IMPLEMENTED 2026-11-24 (report 10, 10.6-10.7), with the statistic REFRAMED.** The claim is NOT that
> MetChange beats expression in a horse race — expression being informative too is expected, since MetChange is
> computed FROM expression. The claim is that **integrating stoichiometry and the network adds information
> expression alone lacks**, which is a question about INCREMENT, so the statistic must be CONDITIONAL.
> `launch_prism_cell.sh` now runs **two arms** (metchange + expr with the GPR gene sets, same drugs/pairs);
> `arm=metchange|expr|both`, default both, passed through by `launch_prism_slurm.sh`; the expr arm self-skips
> with a message if the cell has no `expr` dir or gene sets are missing.
> `prism_differential.py --expr-run prism_viability_expr` adds per row: `r_expr` (same differential on the
> matching GENE set) and **`r_partial_mc_given_expr`** (MetChange differential with the gene set PARTIALLED
> OUT — the statistic that matters) plus `p_partial`. Only the 10 subsystem sets have same-named gene
> counterparts, so `--expr-set-map` pairs the endpoint sets to their parent subsystem (dntp/ntp/nmp ->
> nucleotide, amino_acid_endpoints -> amino_acid, tca_endpoints -> tca, etc.). A summary line reports what % of
> the MetChange magnitude survives conditioning. Verified on a SYNTHETIC expr arm (plumbing only, not a result).
> **PRE-REGISTERED BEFORE RUNNING (report 10, 10.7): primary = `r_partial_mc_given_expr` on `dntp_endpoints`,
> nucleotide-targeting drugs, paired Wilcoxon vs complement over the 10 cell pairs. PREDICTION: retention stays
> well above 50% and the paired test stays positive. REFUTED BY: retention near/below 50%, or the partial
> losing significance while the raw keeps it — either means the association is transcriptional information
> wearing a metabolic label and section 10 stops supporting the thesis WHATEVER its own p-value. NOT claimed
> either way: `r_expr` being large (expression sharing the signal is fine; expression explaining ALL of it is
> not). Secondary: the same columns for `nucleotide`/`purine`/`pyrimidine`/`folate`, whose gene counterparts are
> EXACT rather than mapped — the cleanest like-for-like increment test.**

> **2026-11-24 — FIRST POSITIVE PRISM RESULT: the NUCLEOTIDE arm behaves as predicted once drugs are
> MECHANISM-MATCHED; the AMINO-ACID arm contradicts (report 10, 10).** Rescore with the 9 new sets landed same
> day (82 sets/cell). Drug selection by ANNOTATED TARGET (user's suggestion, made before any of these results).
> **Mean r_diff, prediction = POSITIVE:**
> `dntp_endpoints`: all 821 drugs **-0.036** -> energy+biomass targets (~80/pair) **+0.099** (60% of pairs +)
> -> **nucleotide-pathway targets only (~41/pair) +0.136 (80% +)**; complement **-0.046**.
> Same shape for `ntp_endpoints` (+0.154, 80% +), `nucleotide` (+0.105), `pyrimidine` (+0.091), `purine`
> (+0.081). `folate` null (+0.018). **`amino_acid` -0.100 (20% +), `amino_acid_endpoints` -0.066 — the WRONG
> way, and worse as the match tightens.**
> **Three features: (1) the SIGN FLIPS with mechanism matching — pooled over all drugs every nucleotide set is
> negative; (2) the effect GROWS as the match tightens while n falls 684 -> 41, the opposite of a chance
> finding; (3) the COMPLEMENT behaves as a control (~ -0.05).**
> **ONLY `dntp_endpoints` and `ntp_endpoints` have a single variant, so only their 10 rows are one-per-cell-pair.
> The other sets have 6 variants each — their "60 pairs" are NESTED and must NOT be read as 60 tests.**
> Clean aggregates (nucleotide-targeting drugs, 10 pairs): 8/10 positive (sign p = 0.055), mean Fisher z +0.140,
> **matched-vs-complement PAIRED WILCOXON p = 0.0020** (the right statistic: tests the crossed prediction and
> removes pair-level effects, which is what defeated section 8). Strongest pair HT29/MCF7 +0.374 (p = 0.018).
> Amino acid fails the same test at p = 1.000 against the prediction.
> **CORRECTED ACCOUNTING (report 10, 10.3). All per-pair p-values quoted are UNCORRECTED; the "6e-05" was
> alpha/m = 0.05/820, the BH threshold a p must beat AT RANK 1 — not a q-value.** Actual BH q: **0.62** over the
> 820-test per-pair family, and still **0.176** over the pre-registered `dntp_endpoints`-only family of 10. So
> **no individual cell pair survives, and not narrowly.** (For p=0.018 to pass BH over 820 it needs rank >= 296;
> only 12 tests are at p <= 0.018.)
> **But the AGGREGATE statistic is a different, much smaller family — one test PER SET, m = 9 — and it DOES
> survive: `ntp_endpoints` q = 0.0018, `dntp_endpoints` q = 0.0029; `folate` q = 0.043; `amino_acid` q = 1.000.**
> **DISCARD, do not celebrate, the four sets showing q ~ 1e-09 (`pyrimidine`, `nucleotide`, `nmp_endpoints`,
> `purine`): their Wilcoxon runs over 60 rows = 10 pairs x 6 NESTED compartment/hub variants, so
> pseudo-replicates are treated as independent and the p is inflated by the variant multiplicity.** Only
> `dntp_endpoints` and `ntp_endpoints` have one row per cell pair.
> **The surviving aggregate still rests on 10 pairs sharing only 5 CELL LINES, so effective n < 10 — halving the
> df puts it near p = 0.03, q ~ 0.05. A real pre-specified result at the EDGE of significance, not a comfortable
> one.**
> **PRE-REGISTRATION, apportioned honestly — and the phrase "pre-specified aggregate" is WITHDRAWN.**
> REGISTERED with a written record: the positive direction (7.4, in the script docstring before the first run),
> `dntp_endpoints` as the set to read first (9.5/STATUS, written before the rescore landed),
> complement-as-control (in `--moa-regex` help before any stratified run).
> **NOT registered: narrowing the drug set 90 -> 44 (chosen after seeing the 90-drug result), and — the one that
> matters most — THE AGGREGATE STATISTIC ITSELF, chosen after it became clear the per-pair tests would not
> survive.** So it is an aggregate test on a pre-specified SET in a pre-specified DIRECTION with a statistic
> picked AFTERWARDS.
> **Three defensible aggregates span a factor of 25: sign test (8/10) p = 0.055; mean Fisher z p = 0.005;
> paired Wilcoxon vs complement p = 0.0020. I led with 0.0020, the most favourable.** The Wilcoxon is still the
> best BY DESIGN (only one using the control the crossed prediction needs) but also the one least likely to be
> chosen had it come out worse. **Quote 0.055 as the conservative figure and 0.0020 as the design-appropriate
> one; the gap is the cost of picking the statistic afterwards.**
> The sequence all drugs -> MoA -> any GEM target -> nucleotide targets is FOUR LOOKS at one dataset.
> **NEXT (report 10, 10.4): (1) a HELD-OUT test — PRISM primary screen or CTRP/GDSC, with direction and
> drug-selection rule fixed in advance; (2) DRUG-SPECIFIC target-proximal sets (score each drug against
> metabolites 1-2 reactions from ITS OWN target) — collapses 820 tests to 10, the biggest available power gain;
> (3) `--group-jaccard` (implemented, unused here) — 237 of 747 drugs share a target vector, so n=41 contains
> pseudo-replicates, run it BEFORE quoting any number as final; (4) A375's pairs are the two negative ones for
> `dntp_endpoints`, consistent with the 8.3 coupling anomaly — but be suspicious of a result that improves by
> dropping a cell line.**

> **2026-11-24 — BIOMASS-PRECURSOR SETS ADDED; they were missing from OUR set file, not from Human-GEM
> (report 10, 9.5). CORRECTION: I described building them as needing curation and a scope decision. It needed
> NEITHER.** `SUBSYSTEMS` in `build_metabolite_sets.py` maps a set name to Human-GEM subsystem names and the
> builder takes every metabolite in those subsystems' reactions — and **Human-GEM has 19 nucleotide /
> amino-acid subsystems out of 147** (Phe/Tyr/Trp biosynthesis 157 reactions, Nucleotide 119, Purine 108,
> Arg/Pro 99, Tyr 78, Pyrimidine 75, Folate 64, Gly/Ser/Thr 61, Val/Leu/Ile 53, Cys/Met 44, Ala/Asp/Glu 41,
> Lys 32, + 6 more). It was a lookup.
> **ADDED — subsystem sets:** `nucleotide` (192), `purine` (150), `pyrimidine` (140), `folate` (97), and ONE
> combined `amino_acid` (823) rather than 14 separate ones to stop the BH family doubling. `Amino sugar and
> nucleotide sugar metabolism` deliberately EXCLUDED (glycan-adjacent, belongs nearer NEG_CONTROL).
> **ADDED — curated endpoint sets:** **`dntp_endpoints` (14) is the one to read first** — dATP/dGTP/dCTP/dTTP,
> direct products of RNR and TS, the best-represented metabolic targets in PRISM (RNR 4 drugs, TS 4, DHFR 3),
> and NOT currency hubs so `all_nohub` is also 14. Plus `ntp_endpoints` (20; **no nohub variant exists because
> all four NTPs are hubs**), `nmp_endpoints` (39), `amino_acid_endpoints` (89).
> **Set file now 20 sets / 87 set x variant** (was 11 / 40), `metabolite_sets.tsv` = 62 kB.
> **TWO NAMING TRAPS HIT WHILE DOING THIS: (1) amino acids carry NO `L-` prefix in Human-GEM — `alanine`, not
> `L-alanine`; probing the prefixed form returned 1 of 20 and looked like the model lacked them. (2) DO NOT
> check metabolite names in `data/networks/combined/nodes.tsv` — that is the BUILT network and currency
> metabolites/hubs are filtered out during construction, so ATP, GTP, CTP, UTP and every amino acid are absent
> BY DESIGN; checking there suggested the GEM lacked ATP. Authoritative source is
> `data/GEMs/Human-GEM/model/Human-GEM.yml` parsed as the builder parses it; `model/metabolites.tsv` has only
> cross-reference IDs, NO metName column.** The builder does warn on unmatched subsystems and curated names, so
> typos are loud — but only if checked against the right file.
> **NEXT: rescore on the cluster.** Copy the regenerated `metabolite_sets.tsv` (62 kB) across, or rebuild with
> `rebuild_sets=true cell=ALL bash scripts/mc_phenotype_assoc/prepare_assoc_inputs.sh`, then re-run
> `launch_prism_slurm.sh`. The PRISM phenotype table is unchanged. **When reading the result: 9 more sets means
> a stricter BH correction, and `amino_acid` at 823 metabolites is a seventh of the matrix — close to a global
> measure. Read `dntp_endpoints` first: small, mechanistically specific, matched to the drugs present.**

> **2026-11-24 — PRISM IS A WEAK INSTRUMENT FOR THIS HYPOTHESIS; the section-8 null is NOT a refutation
> (report 10, 9).** Objection raised and correct: growth inhibition has many mechanisms, and a hypothesis about
> metabolite production applies only to drugs acting THROUGH metabolism.
> **Only 45 of 821 PRISM drugs (5%) have a metabolic MoA.** The library is EGFR (22), tubulin (17), HDAC (17),
> MEK (14), topoisomerase (13), mTOR (13), CDK (13), glucocorticoid (13)...
> **SIMULATION — pooling across mechanisms destroys the test: 45 drugs with TRUE r embedded in 776 nulls, over
> all 821 -> true r=0.4 gives pooled r=0.022 (8% power); true r=0.8 gives pooled r=0.045 (24% power). THE RUN
> REPORTED MEDIAN |r_diff| = 0.030 — indistinguishable from a true r of 0.6-0.8 among applicable drugs.** So
> section 8's pooled null cannot separate "no effect" from "a very large mechanism-specific effect".
> **`--moa-regex` added** to `prism_differential.py`: reports each combination for matching drugs AND the
> complement (the built-in control — the hypothesis predicts an effect in the former, not the latter).
> **Stratified result: only 30 matching drugs per pair; 0 of 400 reach BH q<0.05 (min q 0.995), and BH over 400
> needs |r|~0.65 at n=30, so individual rows are meaningless** — and the big values are sign-incoherent across
> pairs for the same set. **The well-powered AGGREGATE comparison is flat: predicted-positive sign 46% matching
> vs 42% complement; median r -0.015 vs -0.010. The hypothesis predicts matching > complement; it does not.**
> (Caveat: the 400 tests are non-independent — 10 pairs over 5 cells, nested set variants — so read it as
> absence of the predicted separation, not an estimate.)
> **The amino-acid arm is UNTESTED for lack of drugs** — the matching compounds are nearly all nucleotide/
> one-carbon directed.
> **REVISED VERDICT: PRISM does not refute the hypothesis, it is structurally a weak instrument for it.** Four
> compounding problems: wrong sets (8.5), cell-level confounding not cancelled (8.3), 95% of the library
> irrelevant (9.2), stratified test underpowered with a flat aggregate (9.3).
> **1.2 called PRISM "a far better-powered version of the cross-cell test" on 821 drugs x 5 cells — counting
> compounds without asking how many were RELEVANT. Same error as counting whole-LINCS overlap without the
> cell-line panel (5.1). RULE: usable n is the number of drugs for which the hypothesis makes a prediction, not
> the number of drugs in the file.**
> **PRISM remains good for:** a well-powered cell-line-SELECTIVITY readout for questions about bioactive
> compounds generally, and the 8.3 diagnostic of how much cell-specific score-to-sensitivity coupling each
> label set carries.

> **2026-11-24 — PRISM DIFFERENTIAL RUN: PREDICTION FAILS, and the design is WEAKER than advertised
> (report 10, 8).** All 5 covered cores scored on the cluster against their own labels (679-760 drugs each),
> 40 set x variant x 10 cell pairs = 400 tests.
> **Of 400, 171 (43%) carry the predicted POSITIVE sign (50% = no signal), and of the 13 at BH q<0.05, ZERO
> have the predicted sign — all negative.** Top: A375/MCF7 `glycolysis cyto_nohub` r = -0.187 (BH q 3.3e-04).
> Per the pre-registration this CONTRADICTS the prediction rather than supporting it with a flipped sign.
> |r| <= 0.19 throughout.
> **ALL 13 significant results are pairs CONTAINING A375**, and mean r_diff per pair sorts by A375 membership
> (A375/* : -0.026..-0.061, 5-35% positive; non-A375 pairs: -0.009..+0.035, 30-80% positive).
> **Not variance dominance:** A375 does have the largest score dispersion (SD 0.328 vs 0.263-0.303) and
> corr(r_diff, naive_r_in_x) = +0.40 (+0.45 for A375), so `--standardise-within-cell` was added (z-score each
> cell per set before differencing) — it changes almost nothing (9 vs 13 significant, still all negative, still
> all A375).
> **THE CAUSE, and a CORRECTION to report 10, 1.2: differencing cancels confounding SHARED between the two
> cells (a drug toxic everywhere gives dP~0) but NOT a CELL-SPECIFIC score-vs-sensitivity coupling.** Those
> within-cell couplings differ systematically: mean Spearman(score, AUC) over 40 sets = A375 +0.071,
> A549 +0.036, HT29 +0.126, MCF7 +0.126, PC3 +0.120. A375 `glycolysis cyto_nohub` is **-0.155** vs +0.059 in
> MCF7, which is exactly the headline r_diff = -0.180. **So r_diff measures the DIFFERENCE IN CONFOUNDING, not
> the biology.** 1.2's "cancelled by construction" claim is now marked corrected in the report.
> **The one thing pointing the predicted way CANNOT be claimed:** within-cell correlations are mostly positive
> (the predicted direction) in 4 of 5 lines, BUT **26-28 of 40 set x variant combinations are significant at
> p<0.01 in HT29/MCF7/PC3** — near-uniform across glycolysis, TCA, OXPHOS, PPP, adenylates and energy currency
> alike. That is general bioactivity, exactly what 1.1 predicted would happen. **Negative-control sets were NOT
> in this run and must be next time** (one extra `--sets`): if `negctrl_o_glycan` also sits at +0.12 the matter
> is settled.
> **THE STATED HYPOTHESIS WAS NOT THE ONE TESTED (report 10, 8.5).** `prism/notes.md` says "Biomass precursors
> (nucleotides, amino acids) -> proliferation", but the 11 sets in `metabolite_sets.tsv` are ALL central
> carbon/energy metabolism — **there is NO nucleotide set and NO amino-acid set** (grep for
> nucleotide|amino|purine|pyrimidine returns 0). `ppp` and `adenylates` are the whole overlap. So section 8
> tests an ENERGY-METABOLISM proxy. 8.1-8.4 stand (the 8.3 design limitation is set-independent) but the verdict
> applies to energy metabolism only. **Second mismatch: `z` is standardised WITHIN each perturbation, so it
> measures pathway SELECTIVITY within a drug's response — a drug that uniformly collapses everything scores ~0
> on every set. For "collapse of biomass precursors stops proliferation" the relevant quantity is plausibly
> TOTAL disruption, which `z` removes by construction.**
> **FIRST FIX, before anything else: build nucleotide and amino-acid sets** (`build_metabolite_sets.py` needs
> Human-GEM subsystem definitions for them). Without them the stated hypothesis cannot be tested at all.
> **VERDICT: PRISM does not support an ENERGY-METABOLISM version of the hypothesis, and does not test the
> biomass-precursor version at all.** It was the best-powered cell-matched test available
> (821 drugs, 5 lines), so this negative carries more weight than the small-n nulls of report 09 8.1 / report
> 10 5. **It does NOT bear on the HepG2 MMP or SH-SY5Y ATP results** — different phenotypes, and crucially
> phenotypes that are not themselves the confounder.
> **For a usable version: (1) include negative-control sets; (2) handle CELL-level confounding explicitly
> (regress each cell's score on its own AUC and difference the RESIDUALS, or match drugs on potency within
> cell) — differencing alone is provably insufficient; (3) a phenotype that is not viability. The 1.1 problem
> was not solved by the differential, only moved.**

> **2026-11-24 — MICROMAMBA ENV FIX + a SILENT-SUCCESS bug (report 10, 7.7).** Cluster runs failed with
> **"No module named numpy"**: deps are in the micromamba env `metabolic_drug_repurposing` and **a Slurm job does
> NOT inherit an activated environment**, so `python3` was the system interpreter. Invisible locally because the
> laptop's system python3 has numpy/scipy.
> **Fixed in ONE place: `scripts/mc_phenotype_assoc/python_env.sh` (new)**, sourced by
> `mc_phenotype_assoc_common.sh`, `submit_assoc.sh` and `run_naia_assoc.sh` — the three that each had their own
> `python3` default, which is why it bit only some.
> **THE ROOT PREFIX IS THE CRUX: the env is under a SHARED PROJECT DIR, not `~/micromamba`, so the working
> procedure is `export MAMBA_ROOT_PREFIX=/work/projects/car_t_methylome/software/micromamba` then
> `eval "$(micromamba shell hook --shell bash)"` then `micromamba activate`.** Without that export micromamba
> looks in the wrong place and says no such env. The resolver reproduces it, with the root recorded ONCE as
> `mamba_root_prefix` in `project_paths.conf` (repeatable; only existing dirs tried, same convention as
> `scratch_proj_dir`). Order: exported `$py` -> current `python3` -> **micromamba asked about itself**
> (`run -n`, then `info --base` for a candidate root) -> each candidate root with MAMBA_ROOT_PREFIX exported,
> hook-then-activate first, else `<root>/envs/<env>/bin/python3` directly (needs neither hook nor micromamba on
> PATH) -> explicit `mamba_env_prefix`.
> **The self-query step was added because the FIRST version FAILED ITS OWN TEST** — micromamba on PATH and
> configured but its root not in the candidate list meant nothing was tried. Verified against a stub micromamba
> in 3 configurations (on PATH w/ unlisted root; MAMBA_ROOT_PREFIX exported w/o micromamba on PATH; neither):
> resolve, resolve, refuse. **`provenance` now prints the interpreter and how it was
> found**, so logs record which Python produced the numbers. `require_python` refuses BEFORE doing work and
> prints the three ways to fix it.
> **THE FIX EXPOSED A WORSE BUG: with no usable interpreter `launch_assoc_slurm_generic.sh` EXITED 0** — a
> green Slurm array task with no output. The arm scripts run several associations in sequence and never
> aggregate exit status, so every `run_assoc` failed in turn and the script still finished cleanly.
> `require_python || exit 1` now sits at the top of `launch_assoc_metchange.sh` and `launch_assoc_expr.sh`.
> Verified every entry point now exits 1 on an unsatisfiable requirement (launch_prism_cell, prepare_assoc_inputs,
> submit_assoc, run_naia_assoc, both generic-slurm tasks) and the laptop path is unchanged.
> **Fourth member of the family in report 10, 7.6: the failure mode that matters is not the error, it is the
> confident-looking success.**

> **2026-11-24 — THREE MORE DATASETS: PRISM, Panorama GCP, Way (report 10).**
> **PRISM (`data/validation/prism/`) is the valuable one — 1,449 BRD compounds x 482 cell lines, 8-dose,
> AUC/IC50/EC50; `broad_id` + `ccle_name` so mapping is EXACT, no structure work. 1,103 of 1,449 in LINCS.
> Covers 5 of 9 cores: A375 (988 LINCS drugs w/ AUC), A549 (1,053), HT29 (1,045), MCF7 (1,043), PC3 (1,008).
> HEPG2, VCAP, HCC515 ABSENT (verified: PRISM has 14 liver lines, none is HEPG2); HA1E absent by construction.
> HEPG2's absence is the headline limit — every positive in reports 07-09 is in HepG2, so PRISM cannot EXTEND
> them, only test other lines.**
> **DO NOT run the obvious design: the phenotype IS viability, which is THE confounder of reports 07-09. You
> cannot match on cytotoxicity when cytotoxicity is the phenotype** — a naive test would likely light up every
> set including the glycan negative controls.
> **RUN INSTEAD the within-drug between-cell DIFFERENTIAL (report 10, 1.2): per drug, correlate
> AUC(X)-AUC(Y) against setscore(X labels)-setscore(Y labels). Differencing CANCELS general potency by
> construction**, leaving only the selectivity question.
> **IMPLEMENTED AND RUNNABLE ON THE HPC — NO LARGE FILE MOVES (report 10, 7).** `build_prism_phenotype.py`
> runs on the LAPTOP (where the ~800 MB PRISM download is) and writes ONE 86 kB table, 821 drugs x 5 cells;
> scp that single file to `<scratch>/results/mc_phenotype_assoc/phenotypes_ALL/`. Then on the cluster
> `cd scripts/mc_phenotype_assoc && sbatch launch_prism_slurm.sh` (array 0-4, one task per covered core;
> `cells="A549" sbatch --array=0-0 ...` for one; `cell=A549 bash launch_prism_cell.sh` without slurm).
> **Bring home only `<CELL>/prism_viability/per_drug_scores.tsv`, ~3 MB each** — the GB label matrices never
> move. Then `prism_differential.py` (runs either side) does the actual test.
> **IGNORE each run's `set_results.tsv`: it is the sensitive-vs-resistant group contrast, i.e. the CONFOUNDED
> NAIVE TEST.** Both launchers say so in their headers and output, and per-cell resampling is deliberately
> minimal (n_perm=200, n_random_sets=10) because that contrast is not the point.
> **PREDICTION RECORDED BEFORE RUNNING (report 10, 7.4), since 9.4 criticises itself for post-hoc pathway
> choice and 10 cell pairs x ~160 combinations makes it worse here: `r_diff` POSITIVE** (lower predicted
> capacity in X -> lower AUC in X; both differences negative together), **strongest for the biomass-precursor
> sets rather than uniform**. A NEGATIVE r_diff CONTRADICTS the prediction; a positive result spread evenly
> over every set incl. negative controls = a confounder that survived differencing.
> Controls in `prism_differential.py`: `perm_p` (drug identity of the predictor differences permuted),
> `bh_q` over the whole family, `r_selective` (the **355** drugs with between-cell AUC spread > 0.1,
> pre-specified), and `naive_r_in_x` (the undifferenced, confounded value, printed as a foil only).
> **CHAIN VALIDATED (report 10, 7.6).** A549 scored end to end on the local corrected labels: 683 of 821 PRISM
> drugs have A549 labels, `per_drug_scores.tsv` = **2.1 MB**. Differential exercised on pseudo-cells in a
> throwaway dir (nothing fake written to `results/`): jittered-null case gave 0 of 40 at BH q<0.05 and 55%
> predicted sign (50% = no signal), as it should.
> **BUG THE VALIDATION CAUGHT AND FIXED: a constant difference vector gave `r_diff = nan` but
> `perm_p = 1/(1+n_perm) = 0.0196` — a spuriously SIGNIFICANT p-value for an input with no information**,
> because `|null| >= |r|` is False for every permutation when r is nan. NOT hypothetical: a set whose members
> are all absent from both matrices gives exactly that, and 23 of 41 sets have members missing from the A549
> matrix. Now such combinations are SKIPPED and listed by name, non-finite permutation values are dropped
> before the p-value, and an all-degenerate run exits with an explanation instead of a table.
> Same family as the HEPG2-labels-for-A549 fallback, the whole-LINCS-vs-panel overlap and the stale matrix:
> **each would have produced a confident number rather than an error, and each was caught by testing a case
> whose answer was known in advance.** **821 drugs with an AUC in all 5 cores and in LINCS;
> 361 of them cell-selective (between-cell AUC SD > 0.1) — the informative subset.** Needs MetChange labels
> for >=2 PRISM-covered cores; locally only A549 has `metchange/rank/log2FC.csv`, the HPC has all 8.
> **PRE-REGISTER the direction** (lower predicted nucleotide/amino-acid capacity in X vs Y -> lower AUC in X vs
> Y): report 09, 9.4 criticises itself for not doing this and with 10 cell pairs the post-hoc risk is larger.
> **Panorama GCP — best rationale in the project, library largely defeats it.** 28 plates, 49 histone-mark
> probes, 24 h, 0.1-10 uM, `pert_id` exact. **7 cell lines with THE SAME compounds in each** (MCF7, A549, A375,
> PC3 cores + NPC, YAPC, Astrocyte) — the crossed structure we have been missing. The rationale's strength is
> that the SIGNS DIFFER BY METABOLITE and are predictable in advance (acetyl-CoA + on acetylation, NAD+ - on the
> same marks, SAM + on methylation, 2-OG - / succinate+fumarate + on methylation). **BUT only 117 compounds and
> the library is built around DIRECT inhibitors of the enzymes whose marks are read** — HDACi (vorinostat, TSA,
> belinostat, MS-275, CI-994, RGFP966), HMTi (EPZ-5676/004777/5687, GSK126, CPI-169, UNC-0321/0646,
> BIX-01294/01338, UNC1215), KDMi (GSK-J4, methylstat), BETi (JQ1, GSK525762A, GSK1210151A), HAT C646, sirtuin
> (EX527, salermide), DNMTi (decitabine, zebularine). **Vorinostat raises H3K9ac by blocking HDACs, not because
> acetyl-CoA rose** — run naively the test "succeeds" for non-metabolic reasons. **Salvage: exclude them
> (~70-80 left, incl. niclosamide = protonophore, pravastatin, bafilomycin A1) and use the excluded drugs as a
> POSITIVE CONTROL for the assay (HDACi must raise acetylation or the extraction is wrong); EX527/salermide
> check the NAD+ arm's SIGN, since blocking NAD+-consuming deacetylases mimics low NAD+.**
> **Way — NOT ACTUALLY DOWNLOADED. All 20 profile files are 133-byte Git LFS POINTER STUBS, and
> `lincs-cell-painting-v1.zip` is the same commit with the same stubs.** Only `metadata/moa/repurposing_info*.tsv`,
> `platemaps/` and `clue_manifest/` are real (worth keeping). Fix: `git lfs pull` in that dir, or fetch the
> profiles from Zenodo 5008187. **Worth getting:** batch `2016_04_01_a549_48hr_batch1` = A549 48 h ~1,500
> compounds, and A549 is the one non-HepG2 core with LOCAL MetChange labels -> a properly powered cell-matched
> test. Caveat from nyffeler carries: only `Mito_*` features bear on the hypothesis, and mito-selectivity
> (mito vs the rest) is the version that controls general bioactivity internally.
> **ORDER (report 10, 4): 1 PRISM differential; 2 `git lfs pull` Way; 3 GCP with the epigenetic arm split out
> and signs pre-registered.**
> **xu — RUN 2026-11-24, and UNDERPOWERED BY CONSTRUCTION (report 10, 5). Only 11 MMP_down vs 12 inactive
> reach the HEPG2 panel** (A375 panel 14/20; full-LINCS ceiling 17/25 — not a runnable option). **At 11/12 a
> true AUROC of 0.350 would give p ~ 0.21**, so the test cannot reach significance even if the effect is as
> large as the one it was meant to check. Phenotype extraction was clean (all 61 actives verified as
> membrane-potential DECREASERS, negative curve rank and efficacy). Result: `tca cyto` 0.295 (p = 0.10,
> randset 0.16), `tca mito_nohub` 0.311, `oxphos mito` 0.470 — nothing survives any control.
> **READ THE DIRECTION: the lean is the DAMAGING one, not the reassuring one.** A positive xu result would mean
> HEPG2 labels predict MMP in a cell line they know nothing about; the point estimates lean weakly that way,
> through TCA — the SAME pathway as HepG2's own MMP result, where the SH-SY5Y ATP transfer ran through
> glycolysis. At 11/12 it is a hint only.
> **Cause: chemical-space mismatch.** Xu screened Tox21 ENVIRONMENTAL chemicals (4-nonylphenol,
> 2-ethylhexylparaben, kadethrin) and LINCS is drug-focused — the same cause as the Hallinger null (report 07).
> `map_smiles.R` cannot rescue it: a handful more compounds does not change the power.
> **MISTAKE TO REMEMBER (report 10, 5.1): xu was recommended twice on a WHOLE-LINCS overlap of 53 without
> checking the HEPG2-panel count — the same error as naia's 683 -> 10, made a day later. RULE: a validation
> dataset's usable n is its overlap with the SPECIFIC label panel, computed BEFORE the dataset is ranked,
> never its overlap with all of LINCS.** The falsification role falls to the PRISM differential, which has 821
> drugs and cancels the viability confounder; xu could not have controlled it anyway, as Xu reports NO
> viability readout.

> **2026-11-23 — STRONGEST RESULT SO FAR: MetChange beats expression under a MEASURED viability control, and it
> REPLICATES (report 09, 9.6).** Expression labels ALSO predict SH-SY5Y ATP loss, through MITOCHONDRIAL sets
> (`oxphos gene` 0.314-0.348, random-set p = 0.015; `glycolysis gene` null) — the reverse ordering to MetChange,
> so the two layers are NOT reading the same thing, and expression's ordering is the biologically expected one.
> **But matching on naia's MEASURED cytotoxicity: expression DIES (oxphos 0.348 p = 0.013 -> 0.329 p = 0.115)
> while MetChange SURVIVES (glycolysis 0.265 p = 1.1e-04 -> 0.173 p = 0.0025; ppp -> 0.156 p = 0.0014).**
> **This is section 4's pattern independently replicated across FOUR differences:** phenotype (MMP -> ATP),
> phenotype cell line (HepG2 -> SH-SY5Y), library (Tox21 -> Prestwick), and the confounder control (`n_sig`
> PROXY -> MEASURED cytotoxicity). Matched: MetChange 0.365 / 0.173, expression 0.500 / 0.329.
> **This is a much better basis for "the stoichiometric layer carries metabolic information the transcriptome
> does not" than report 08 section 13**, which compared UNMATCHED magnitudes in one cell line and found no
> advantage. **Limits: 15 matched pairs, |AUROC-0.5| = 0.15-0.33 (small), pathway post hoc in both arms. And two
> failures still unexplained and NOT to be dropped: Schuhknecht per-drug metabolite correlation null
> (rho = -0.013) and Hallinger OCR null.**

> **2026-11-23 — FALSIFICATION TEST RUN: POSITIVE, but through a DIFFERENT PATHWAY (report 09, 9).** Naia
> SH-SY5Y ATP phenotype vs **HEPG2** MetChange labels (the converse of the section-3 design: labels fixed,
> phenotype cell line varied). 25 ATP_down / 27 ATP_up / 226 inactive.
> **`glycolysis cyto` AUROC 0.265, MW p = 1.1e-04, random-set p = 0.005, BH q = 0.005**; glycolysis_endpoints
> 0.280, ppp 0.303, pyruvate 0.306. **TCA is nominally significant but FAILS the random-set control**
> (p = 0.08-0.19); `oxphos mito` null (0.428); **`energy_currency` NULL (0.482) even though the phenotype IS
> ATP**; `redox_couples` moves the OPPOSITE way (0.648). **ATP_up is null everywhere** (all BH q = 0.94).
> **NOT cytotoxicity:** naia ships a measured per-compound cytotox z-score, so viability is controlled directly
> rather than via `n_sig`. The imbalance was +6.62 z-units; matching removes it and **the effect STRENGTHENS to
> 0.173 (p = 0.0025, 15 pairs)**; joint cytotox+n_sig 0.222 (p = 0.010). Leave-one-out never moves p above
> 3.1e-04; dropping the 5 most extreme cases still 0.321 (p = 0.008); MoA heterogeneous. Case group is
> chemically plausible (hexachlorophene, triclosan = uncouplers; proscillaridin/digoxin/digoxigenin = cardiac
> glycosides).
> **CONSEQUENCE FOR REPORT 09 SECTION 3 — the strong reading is REFUTED.** HEPG2 labels are NOT specific to
> HEPG2-measured phenotypes; they transfer. The asymmetry predicted in 8.1 is observed. **But it is NOT the same
> signal:** MMP->TCA/mitochondrial vs ATP->glycolysis/cytosolic, `energy_currency` null in both, and a generic
> chemotype axis would have produced the SAME sets. **"The association is HepG2-specific" is no longer sayable
> without naming WHICH association.**
> **THE DEPTH EXPLANATION (reading 2 of section 6) IS DEAD: A375 has 13,314 profiled compounds vs HEPG2's 5,494
> — 2.4x more — and A375 is the MOST NULL line in section 3 (0.510, p = 0.73).** HEPG2 also has the LOWEST
> responsive share of the eight (28%/7%), so its labels are not "more significant" either.
> **BIGGEST REMAINING CONFOUNDER, uncontrolled:** the glycolysis hit is biologically odd — naia used
> DIFFERENTIATED (oxidative) SH-SY5Y, so mito sets were the a priori expectation. Glycolysis may mean the HEPG2
> labels report **HEPG2's own glycolytic wiring** rather than the phenotype cell's biology, a label-side
> artefact correlating with ATP loss. Also: n = 25 cases / 15 matched pairs, and **the pathway assignment is
> POST HOC over 160 set x variant x statistic combinations** — the 9.3 comparison would be far stronger if
> glycolysis had been predicted in advance.
> **ONE COMMAND REPRODUCES ALL OF SECTION 9 (report 09, 9.7): `bash scripts/mc_phenotype_assoc/run_naia_assoc.sh
> [CELL]`** (default HEPG2; chains phenotype build -> cross-tab -> MetChange -> expression -> matched controls).
> New scripts: `build_naia_phenotype.py`, `covariate_matched_check.py` (matches on any numeric phenotype column,
> generalising `matched_check.py`; `--robustness` adds leave-one-out / drop-5 / MoA / case listing),
> `compare_phenotype_tables.py`, `run_naia_assoc.sh`.

> **2026-11-23 — THREE NEW VALIDATION DATASETS assessed (report 09, 8).** `data/validation/{naia,nyffeler,xu}/`.
> LINCS overlap below is the InChIKey-first-block FLOOR; `map_smiles.R` should recover more.
> **naia — DISQUALIFIED for the symmetric test (resolved 2026-11-23 by the gating check). LINCS has only 224
> perturbations in SHSY5Y vs 5,494 in HEPG2, and only 10 naia compounds are in that panel — 1 ATP-down, 0
> ATP-up, 4 flat** (the lone ATP-down is D-cycloserine, a peptidoglycan-synthesis inhibitor, i.e. exactly the
> implausible-hit noise flagged in 8.1). The 683 whole-LINCS matches said nothing about the cell line that
> matters, and the differentiation worry is moot. **naia's remaining value: it REPLACES xu as the falsification
> test with 6x the power — 340 of its compounds are in the HEPG2 panel**, so ask whether HEPG2 labels predict an
> ATP phenotype measured in another cell line. NOTE THE DIRECTION: report 09 section 3 fixes the phenotype and
> varies the LABEL cell line; this fixes the labels and varies the PHENOTYPE cell line. A positive would be an
> odd asymmetry pointing at **reading 2 of section 6 — HEPG2 labels are just the deepest/best, and MetChange
> works where LINCS data is deepest** — the alternative report 09 cannot exclude.
> *(Original naia assessment, now moot, kept in report 09, 8.1a:* `SHSY5Y` IS
> in `cellinfo_beta`, the phenotype is ATP (a metabolite in `energy_currency`, now producer-only), and it ships
> an INDEPENDENT cytotoxicity z-score - a better confounder control than the `n_sig` proxy, and genuinely
> independent (Spearman(ATP, cytotox) = -0.083, p = 0.004, so ATP loss here is NOT cell death). Weaknesses:
> single concentration, noisy, enhancer-screen design.)
> **naia CAVEAT that could disqualify it: the screen used DIFFERENTIATED SH-SY5Y, LINCS `SHSY5Y` is
> UNDIFFERENTIATED.** Differentiation shifts SH-SY5Y from glycolytic to oxidative metabolism - exactly the axis
> an ATP phenotype measures - so "cell-matched" would mean matched on cell LINE, not on metabolic state, and a
> null could just mean the labels describe a glycolytic cell and the phenotype an oxidative one. Check the
> differentiation protocol in the methods before treating naia as the symmetric test.
> **xu — AC16 cardiomyocytes, MEMBRANE POTENTIAL (the same phenotype as Attene-Ramos), 247 compounds, 53 match.
> AC16 is NOT in LINCS**, so no cell-matched test is possible. Its value is as a **FALSIFICATION TEST of report
> 09 section 3**: if HepG2 MetChange predicts AC16 MMP about as well as it predicts HepG2 MMP, the signal
> travels with the compounds, the chemotype explanation revives and the cell-specificity claim is in trouble.
> Run it BECAUSE it can only damage the current conclusion - and note a positive is informative while a null is
> not, at n = 53.
> **nyffeler — U2OS cell painting with 10 `Mito_*` categories, CONCENTRATION-RESPONSE (BMD), 462 compounds.
> Best design, worst overlap:** only 87 match LINCS and just 19 of those are mito-selective
> (`min(Mito_* BMD) - min(other BMD) < 0`), so no binary design; and 457 of 462 have a responding `Mito_*`
> category, so there is effectively NO inactive group. Usable only as a continuous rank correlation over 87
> drugs. Its mito-selectivity contrast is attractive because it nets out general bioactivity internally.
> **GATING CHECK before any of this (report 09, 8.4): `cellinfo_beta.txt` lists every line LINCS ever touched
> and says NOTHING about depth.** Count distinct `trt_cp` `pert_id`s for SHSY5Y and U2OS in `siginfo_beta.txt`
> on scratch (awk one-liner in 8.4). If SHSY5Y has only a few dozen compounds profiled, naia is dead despite
> the 683 structure matches. Adding either line also means a new MetChange + limma label run for it.

> **2026-11-23 — CROSS-CELL: the association is HepG2-SPECIFIC (report 09).** The label-level test run over the
> 8 core cell lines with MetChange labels, all against the single union phenotype table, then restricted to the
> **641 drugs present in every cell line with a consistent group (125 MMP-down / 427 inactive)** so only the cell
> line varies. **MetChange `tca mito_nohub` z: HEPG2 0.350 (p = 3.3e-07); A549 0.493, A375 0.510, HA1E 0.434,
> HCC515 0.481, HT29 0.464, PC3 0.478, VCAP 0.485 — all null.** Responsiveness-matched per cell: HEPG2 0.365
> (p = 3.4e-04), everything else null. **This RULES OUT the compound-class / chemotype explanation of report 08
> section 2** — identical compounds cannot act in one cell line and not in seven — and it is not power (HEPG2 has
> the FEWEST drugs, 717 vs up to 1,364, and the restriction makes its effect stronger). PC3's earlier nominal
> signal (p = 0.036) was coverage and disappears on the common set (p = 0.46). **A prediction stated in advance**
> (uniform = chemistry, HepG2-only = HepG2 biology).
> **CORRECTED 2026-11-23 (report 09, 3.1 / 4 / 6):** (a) **MMP-UP is NOT HepG2-specific** - `energy_currency
> mito` survives matching in HCC515 too (0.624, p = 0.011), A549 borderline (p = 0.062); the specificity claim
> holds for **MMP-down only**, and the first version of the report generalised from the MMP-down run alone.
> (b) The expression arm is **not** "cell-general" - that came from POOLED numbers; under matching it collapses,
> **in HEPG2 to exactly 0.500 (p = 0.99)**, surviving only in PC3 (0.404, p = 0.0087) and HCC515 (0.408,
> p = 0.012). **The claim to use instead, within one cell line: in HEPG2, matched, MetChange `tca mito_nohub`
> 0.365 (p = 3.4e-04) vs expression `tca gene` 0.500 (p = 0.99)** - MetChange survives the control that removes
> expression entirely, no cross-cell argument needed. (c) HA1E's nominal p = 0.025 does NOT survive matching
> (0.454, p = 0.21); it has the largest responsive share of any line (64% vs 37%) and ~0.2 after Bonferroni over
> 8 cells.
> **The o_glycan negative control is also NOT cell-specific** (HEPG2 0.590, HA1E 0.582, VCAP 0.578), supporting
> the glycosyltransferase-expression reading of report 08 3.3 and showing the design discriminates rather than
> flagging everything.
> **Still open.** The "rules out chemistry, not transcription" caveat is **WITHDRAWN** - MetChange IS an LP over
> the transcriptional signature and the stoichiometric matrix, so "transcription projected onto TCA" describes
> the method, not an alternative. **The real limitation is that MMP was measured only in HepG2**, so two
> readings remain and nothing separates them: the phenotype is itself cell-specific and MetChange tracks it
> correctly (nulls are correct negatives, result strong), or MetChange only works in HepG2 and the phenotype is
> shared (nulls are failures). Also: HEPG2 was the discovery cell, so the evidence is the PATTERN not its own
> p-value; effect size unchanged at |AUROC-0.5| ~ 0.15.
> **Decisive next experiment: a second (phenotype, cell line) pair.** `APR_HepG2_MitoMembPot` cannot do it -
> same line. **Schuhknecht (A549 metabolomics, on disk) was proposed and WITHDRAWN**: its phenotype would be
> measured metabolite LEVELS over the same metabolites the MetChange set summarises, so the group test is the
> already-null per-drug correlation (report 05, 8.2, rho = -0.013) coarsened into bins - and report 05, 9.3
> (capacity is not validated by pools) applies in full. **Promoted instead: the LINCS GENETIC PERTURBATIONS**
> (`trt_sh`/`trt_xpr`/`trt_oe` exist PER CELL LINE, so they give (phenotype, cell) pairs across many lines with
> no new data, exact cell-independent ground truth, and no metabolite-level readout; the >=1-reaction-away
> version is the real test since MetChange reads the perturbed transcript itself). Secondary: PRISM/DepMap
> viability across all 9 lines - not metabolic, but tests the cell-specificity MECHANISM at large n. **Honest
> gap:** no mitochondrial phenotype at screen scale in a non-HepG2 core line is known to us - everything touched
> (Attene-Ramos, Hallinger, Apredica, CCTE_Simmons) is HepG2, which is itself informative. **DONE:** the
> common-drug restriction is now `scripts/mc_phenotype_assoc/cross_cell_compare.py`, which reproduces every
> section-3 number.

> **2026-11-20 — confounder analysis of the phenotype association (report 08).** **The association is carried
> entirely by the drug-like part of the library**: restricted to the 115 compounds shared with Hallinger
> (= ToxCast) the Attene effect is 0.447 (p = 0.39, n = 38/55), outside them 0.352 (p = 1.8e-06, n = 109/425).
> **That, not n, is why Hallinger is null** - its point estimates for the Attene-signal sets are 0.47-0.55, not
> 0.37, and power at 29/86 for 0.371 is 54%. **Mechanism class does NOT explain it:** within 9 coarse MoA
> families `tca mito_nohub` is 0.380 (Holm p = 0.0054) and `energy_currency mito` 0.657 (Holm p = 0.0010),
> essentially unattenuated; the negative control `oxphos mito` stays null everywhere. Narrow (fine) MoA strata
> give p = 0.13-0.57 but are 27 tiny strata in a receptor-drug population - untestable, not negative.
> **No potency gradient survives correction** (best Holm 0.076, and the negative-control set shows an
> equal-sized gradient), so there is no dose-response support. Similarity matching is uninformative (median
> Tanimoto 0.27). **An earlier within-MoA reading of mine ("fits a compound-class effect") was an artefact of
> double-counted strata, an unpaired test and no correction - withdrawn, see report 08 section 3.**
> **Dose/potency line exhausted (report 08, 5.1-5.3):** `Ratio Potency` is the qHTS AC50 (half of the compound's
> own fitted max, ~15-point titration to ~92 uM); across actives it spans 3.9 decades but the interquartile
> range is only 0.75, and **the max LINCS dose has an interquartile spread of 0.00 decades (~all 10 uM)**.
> **48% of MMP-down cases were profiled BELOW their own AC50** - a real design fact - but splitting on it gives
> identical effects (`tca mito_nohub` 0.361 exposed vs 0.367 sub-threshold, p = 0.68, rho = +0.002), and where a
> gradient appears (`tca cyto sign_dev` rho = -0.193) the negative control `oxphos mito` shows a stronger one
> (-0.217). The test is also NOT independent of the potency test (constant dose => driven by AC50 alone) and the
> MMP AC50 is not a threshold for a 6-24 h metabolic change. So this dataset cannot provide a dose-response test
> at all; no further potency analysis of it is worthwhile.
> **DECLARED NEGATIVE CONTROLS change the interpretation (report 08, 3.1).** `oxphos mito` was NEVER a valid
> negative control - it is the most on-target set for this phenotype, was designated post hoc because it came
> out null, and is 11/13 hubs; it is a POWER DIAGNOSTIC and is relabelled as such in reports 07 and 08. Five
> proper controls (glycan / glycosaminoglycan subsystems, criteria fixed before looking):
> **4 of 5 are null, but `negctrl_o_glycan` gives AUROC 0.596-0.631 (p = 4e-4 to 1e-6), beats its size-matched
> random null (rand = 0.0099, same as TCA's), and survives coarse-MoA (0.613, Holm 0.0097) and similarity
> matching.** Same magnitude as TCA (|0.13|), opposite direction. **So beating a size-matched null is NOT
> sufficient to implicate a pathway, and the defensible claim is a MULTI-PATHWAY difference between
> phenotype-active and inactive compounds (TCA down, O-glycan up, 34/40 primary sets below 0.5), not a
> TCA-specific metabolic response.** **Tested the obvious rescue and it FAILED (report 08, 3.2):**
> if o_glycan were coupled to TCA via the hexosamine pathway its nucleotide sugars (UDP-GlcNAc/GalNAc/Gal,
> CMP-NeuAc, serine) would carry the signal - they are **null** (0.488, p = 0.66) and the signal is entirely
> in the **glycan dead ends** (core 2-8, T/Tn/sialyl antigens, degree 2-12: 0.649, p = 4e-08).
> **But the same split on TCA goes the other way:** `tca_core_cycle` (citrate, AKG, succinate, fumarate,
> malate, OAA, isocitrate, cis-aconitate, succinyl-CoA; degree 8-39) gives 0.367 (p = 1.1e-06, rand = 0.005),
> stronger than `tca_peripheral` 0.393. So the two associations are **different in kind** - TCA in high-degree
> central metabolites (a capacity interpretation means something), o_glycan only in low-degree terminal
> structures (MetChange there is ~a readout of a few glycosyltransferase genes). **General diagnostic from now
> on: for any significant set, check whether the signal sits in its high-degree or its low-degree members.**
> **Refined (report 08, 3.3): per member the two effects are the SAME SIZE** (median |AUROC-0.5| 0.100 TCA vs
> 0.097 o_glycan; 8/9 vs 21/21 members at p < 0.05), so o_glycan is not the weaker finding - the difference is
> interpretability. **Important caveat on `rand_p` throughout reports 07-08:** the o_glycan cores are
> near-collinear (mean pairwise rho 0.799 vs 0.511 for TCA), i.e. ~1 effective d.o.f. measured 21 times, so a
> near-collinear set beats a size-matched random set partly BECAUSE its average is less noisy. The random-set
> control assumes comparable within-set correlation; always report per-member statistics alongside set ones.
> The o_glycan set is mucin-type O-glycan cores / Thomsen-Friedenreich antigens, mostly **lumped glycan classes**
> plus 4 **GEM placeholders** (`[protein]`, `[protein]-L-serine`, the O-glycosylation-site residue) which carry
> the same AUROC as the real glycans - so it reads as **glycosyltransferase expression**, not metabolite capacity
> (its nucleotide-sugar precursors are null, so substrate supply is unchanged).
> **Placeholders cannot reach the model** (user, 2026-11-20): the 291/5,707 pool and placeholder features are in
> the limma matrix only because MetChange is computed for every model metabolite and limma ran on the whole
> matrix; they are filtered during network construction along with the currency metabolites (hubs, no real
> biology), so they cannot be sources or queries. Verified for this set: the three `[protein]` entities are
> absent from `networks/combined/nodes.tsv` - but 18/21 remain nodes, including the O-glycosylation-site residue
> and the lumped glycan classes, so lumped pseudo-species are NOT filtered. Cheap check outstanding: how often
> pool/lumped names appear among the 100 selected sources.
> **Decisive next experiment: the LINCS genetic perturbations - a knockdown has no chemotype, no MoA and no
> dose-response requirement, so every confounder here is unavailable by construction.**

> **2026-11-19 (FINAL, supersedes the banner below) — responsiveness-matched design: MetChange survives, the
> expression baseline does not (report 07, section 16).** `matched_check.py` (1:1 greedy matching on
> log1p(n_sig), caliper 0.25, plus a residual-on-n_sig variant; balance: 140/147 matched, imbalance p = 0.64
> MetChange and 139/147, p = 0.97 expr). **MetChange `tca mito_nohub` z: pooled 0.364 -> matched 0.371
> (p = 2e-04, q = 0.0039) -> residual 0.388; `tca cyto` sign_dev 0.415 -> 0.422 (p = 0.025). Expression
> `tca gene` z: 0.413 (p = 0.0014) -> matched 0.486 (p = 0.70, q = 0.87); `oxphos gene` 0.410 -> 0.450
> (p = 0.15).** So **the expression association is essentially all responsiveness and the MetChange TCA one is
> not** - both designs agree, in both statistics, and the `oxphos mito` negative control stays null throughout
> (so it is not lost power). **This corrects section 13.5's "the baseline matches MetChange".** The pooled
> effects were similar in size but different in kind. Also: `rand_p` does NOT depend on set size (rho = +0.10,
> p = 0.52), and at matched absolute size MetChange still beats its null (42 members, rand = 0.0099) where expr
> does not (45, rand = 0.209) - so the universe-fraction worry of 13.3 is weak, and the matched design does not
> use random sets at all. **Claim as it stands (narrowed by 16.5): for MMP-DOWN (147 drugs) the MetChange TCA
> signal survives matching (0.371, p = 2e-04, q = 0.004) and the expression one does NOT (0.486, p = 0.70);
> for MMP-UP (75 drugs) BOTH survive** - MetChange via energy currency (0.651, p = 0.0014, **magnitude only**,
> sign_dev null) and expr via pyruvate genes (0.637-0.640, p = 0.003-0.004, both statistics). So the
> stoichiometric reading is distinguishable from the transcriptome on ONE phenotype, not both. MetChange
> `oxphos mito` is nominally significant for MMP-up only (0.594, p = 0.046), the hub-dominated set, consistent
> with the currency result rather than independent of it. Still not done: a second cell line for the confounder
> structure. **Null-shape comparison DONE (report 07, section 17):** at matched universe fraction the expression
> null is **1.3-1.6x wider** (median |AUROC-0.5| 0.047 vs 0.029 at ~0.22%; 0.042 vs 0.033 at ~0.74%), flat in
> size within each space. So random GENE sets separate the groups better than random metabolite sets - the
> transcriptomic group difference is **diffuse**, MetChange's is **concentrated**: MetChange `tca mito_nohub`
> |effect| 0.136 vs its null p95 0.090 (clears), expr `tca gene` 0.087 vs 0.121 and `oxphos gene` 0.090 vs 0.112
> (both inside). **Claim to use:** the stoichiometric projection CONCENTRATES phenotype-relevant signal into an
> interpretable pathway coordinate where gene space has it diffuse and responsiveness-confounded. **Stoichiometry DOES add information** (report 07, 17.4 corrects an
> earlier objection of mine in 17.3): MetChange = LP(expression, GPR, stoichiometric matrix, directionality) via
> COBRA/Gurobi, and the GEM is a prior NOT derivable from LINCS - a gene-level pathway mean cannot encode that in
> A->B->C a fall in B interrupts flux to C. The comparison is well posed for this, since the expression arm
> already uses the GEM's GPR genes and subsystems, so what differs is stoichiometry/topology/LP. Defensible form:
> *on MMP-down, expression + stoichiometry separates the phenotype (0.371 matched, q = 0.004, clears its null)
> where aggregating the SAME genes over the SAME subsystems does not (0.486, p = 0.70, inside its null)*.
> Surviving caveats: MMP-down only (MMP-up both survive); the comparator is subsystem-MEAN aggregation, so the
> claim is not "beats any expression-based method" (a topology-aware gene method such as SPIA could narrow it);
> the MoA/chemotype confounder; one cell line. **Do NOT claim MetChange "properly
> models" metabolic change** - this is a group-level shift (|AUROC-0.5| ~ 0.14), and the per-drug test (report 05,
> 8.2, rho = -0.013) points the other way. **Untested and the likeliest alternative:** MMP-active compounds may be
> enriched for a MoA/chemotype whose transcriptional response projects onto TCA; matching on n_sig equalises how
> much a drug responds, not what kind of compound it is. A MoA- or structure-matched analysis is the next step
> before presenting this.

> **2026-11-19 — the expression baseline is run; POOLED it matches MetChange (report 07, section 13; see the
> matched-design banner above, which supersedes this reading).**
> HEPG2 expr `log2FC.csv` arrived; `value`+`sign` on the GPR gene sets, same contrasts. **MMP-down: expr
> `oxphos gene` AUROC 0.410 (z, q = 0.0035) / 0.421 (sign_dev, q = 0.018) vs MetChange TCA 0.359 / 0.415;
> MMP-up: expr `pyruvate gene` 0.615 vs MetChange currency 0.660.** So **on raw separation the stoichiometric
> reading has NO advantage over reading the genes** - the kNN verdict (report 06, 6.1), now for the labels.
> MetChange looks more specific on two criteria only: it beats size-matched random sets (rand 0.005-0.02 vs
> 0.13-0.16) and its directional result is not responsiveness-confounded (0.453 ns vs expr 0.335-0.363,
> p <= 2e-7). **Both criteria are partly artefacts of comparing two feature spaces** (20,186 genes / sets of
> 23-106 vs 5,707 metabolites / sets of 5-73; and the confounder test has 265 vs 215 controls for expr but only
> 44 vs 436 for MetChange, so MetChange's "no confound" has CI ~[0.36, 0.55]). **Most interesting: the two label
> types point at different pathways** - expr at the OXPHOS *genes* (on-target, 106 genes no hubs), MetChange at
> TCA (its own oxphos set is 24/26 hubs and null everywhere). Defensible claim for a talk: the MetChange labels
> carry a TCA-capacity signal separating MMP-active from inactive drugs, at a magnitude the transcriptome also
> reaches via OXPHOS genes. Hallinger null in both arms (n = 29). Next: responsiveness-MATCHED design, random
> sets matched on universe fraction, a second cell line, and the orthogonal `APR_HepG2_MitoMembPot` readout.
> **Set construction and a multiple-testing correction (report 07, section 14):** 11 base sets x 3 compartments
> x 2 hub variants = 66 candidates, 26 below the 5-member floor (structurally: currency sets have no non-hub
> member, `oxphos` is 24/26 hubs, `ppp mito` is empty), leaving 40 - of which **39 are distinct**; the gene file
> has 10 rows but **5 distinct sets** (every `gene_nohub` equals its `gene` set). Duplicates are now collapsed,
> with an `aliases` column. **These duplicates did NOT inflate any q-value**: BH is invariant to duplicating the
> whole family (`q = p*m/rank`, both scale), so the section 13 expression q-values are correct as published; only
> the MetChange family (40 -> 39, non-uniform) shifts, in the 4th significant figure. Statistics are per set
> (features aggregated within a set within a perturbation, then per drug), never pooled across subsystems - which
> is what makes the tca/pyruvate/glycolysis/oxphos gradient readable.
> **Run inventory complete (report 07, section 15):** both phenotypes x both label types x all three modalities.
> The last gap, `hallinger_expr_sig`, is null (0.582, p = 0.25), so **every Hallinger arm is null in both label
> types and all modalities** (n = 20-29). Clarification: in the POOLED expression `sig` run `glycolysis gene`
> does survive BH (q = 0.016) - section 12.3 said otherwise, which was true only of the responsive run. It is
> still not evidence (rand = 0.144; responsiveness confounder 0.230 at p = 1.5e-26; q = 0.106 once restricted),
> so the conclusion holds but the earlier statement was incomplete. Also note `attene_expr_sig*` (702 drugs) and
> `attene_expr_valuesign` (726) differ because only the latter carried the cytotox contrast.

**Mapping tool:** `scripts/benchmarking/map_smiles.R` (rewritten 2026-11-19): largest-fragment stripping (salts),
all pairs kept, NA never matched, unique SMILES canonicalised once, and `--stereo both` = exact pass then a
stereo-insensitive retry of the unmatched, every pair labelled with `match_type` / `ambiguous_stereo` so the cost
of merging enantiomers is measured, not assumed. `data/validation/hallinger/DSSTox_CCD_dump_12092025` (CSV)
gives DTXSID / CASRN -> SMILES / InChIKey offline from a dated release. Clowder downloads:
`https://clowder.edap-cluster.com/api/files/<file-id>/blob` (**HEAD 404s, GET works**).

**Watch out: 5,785 of 34,419 LINCS compounds (564 of 5,494 in HepG2) have neither SMILES nor InChIKey** — mostly
unnamed `BRD-…` codes. They can never be mapped to an external dataset by structure, so the mappable HepG2
denominator is 4,930, not 5,494.

**Currency metabolites** (ATP, ADP, NAD+, NADH, NADPH, O2, H+, CoA) **are MetChange features but not network
nodes** (they are among the 339 dropped before `select_sources`). So an energy-metabolism phenotype can be tested
at the **label** level but cannot be built as a model **query**; queries must use the glycolysis / TCA proxies
(lactate, pyruvate, succinate, citrate, fumarate, acetyl-CoA, glucose). This is why the enrichment test comes
before the ranking test for MMP / respiration phenotypes.

**To do (Attene-Ramos, report 07 section 6):**
1. **Viability columns of `AID_720637_datatable.csv`** (laptop): is the cytotoxicity split possible? Most
   compounds that collapse MMP just kill cells; the informative contrast is **MMP-down without viability loss**
   vs inactive. Without this the test predicts cytotoxicity, not metabolism.
2. **Real n** (HPC): mapped compounds x *active* HepG2 perturbations — not 904.
3. **Enrichment test** (HPC, HepG2 MetChange labels at `<scratch>/results/labels/HEPG2/metchange/rank/`, ~1.5 GB
   per matrix so do not download): per drug, summarise MetChange over glycolytic vs TCA/ETC metabolites
   (currency metabolites included), compare MMP-active vs inactive (Mann-Whitney / AUROC), with the expression
   baseline read the same way, and **no assumed sign** (the compensation hypothesis of report 05, 9.2).
4. Aggregation rule for the 1,379 Attene-Ramos CIDs collapsing onto 904 compounds; run first with `exact` +
   unambiguous pairs, then test whether the 225 ambiguous stereo pairs change anything.
5. **Hallinger:** email Simmons (EPA/CCTE) for the `EPAPLT0033` / `EPAPLT0034` plate maps — the only remaining
   source. The 17 GB MySQL `sample` table is a long shot (tier 1 was never registered) and not worth it on spec.

**6.3 Size the data ceiling** (no GPU, ~a day):
1. PCA of the A549 expression matrix: how many components explain 95% of variance across
   perturbations. That is the budget of distinguishable perturbation states.
2. Fraction of the GPR genes MetChange uses that are landmarks vs CycleGAN-inferred. If mostly
   inferred, the metabolic labels are a smooth function of the 978 landmarks.

**6.4 Then** freeze the recipe (done 2026-11-15: `_min15`), and only after that revisit version (b), the prior + residual
model, or per-protein standardisation of the dynamic features. (Done on 2026-11-14: cap check —
does not bind; second seed — replicates, smaller; port to `tf` / `expr` — report 03.)

## 7. Open issues

- **Per-drug plain ranking is still prior-dominated** (section 3). This is the main scientific
  gap for the tool, not a bug.
- ~~MetChange reads the signature about half as well as expr~~ — **resolved**: the label sign bug
  (report 05, section 5; corrected MetChange ≈ expr ≈ TF).
- **Laptop expression / MetChange labels are an old partial set** (`results/labels/A549/expr/`
  and `metchange/rank/`, 2 Nov: 699 of the 3,970 scored perturbations). Analyses needing labels
  must run on the cluster or after syncing the scratch copy. (The MetChange `_min15` control's
  low folds 0–1 did not recur in seed 1.)
- **MetChange similarity-check anomaly — probably explained** (report 04, 1.1). TOP1 / JAK1 / XIAP
  shared-target pairs were anticorrelated and the split matched flip status exactly: consistent
  with the uncorrected home labels (camptothecin flipped, topotecan not). ATP1A1 (mixed-flip,
  r = +0.92) is not explained. Confirm by rerunning the similarity check on the scratch labels.
- **Model does not use edge signs.** `propagation_operator` gates on weight × activity only;
  the two channels propagate independently.
- **Signatures are discretized** (`np.sign`, with `--ignore-qvalue --top-n 100`), while the
  similarity check uses continuous Spearman. A `--signature-weights {sign,value,t}` option was
  proposed, not implemented.
- **`--no-globals` also drops `log_n_sig`** (signature size), which is *not* design information
  — it is computed from the signature and travels with it in the shuffled control. Splitting
  the flag (e.g. `--no-dose-time`) would keep the most direct effect-strength input. Small
  change, not done.
- **Validation `spec` used for epoch selection is unstratified**, so it can reward dose/time in
  any run that keeps those inputs. Harmless for `--no-globals` runs.

## 8. Uncommitted

- Report 04 (`reports/2026-11-16_labels_network_coverage.md`) and the correction banner on report 03.
- `scripts/predict_targets/check_metabolic_coverage.py`: `--coverage` (reuse a `coverage.tsv`),
  fixed 0 / 1–99 / 100 bins for `n_sig_topn`; README lines.
- `scripts/predict_targets/submit_train_w_shuff_mc_min15_fixedlabels.sh` (new).
- `launch_predict_targets_A549_{mc,tf,expr}_noGlobals_lr1e-3_min15_qsig{,_shuff}.sh`, wrapper
  `submit_train_w_shuff_noGlobals_lr1e-3_min15_qsig.sh`, variants registered in
  `submit_with_shuffle_control.sh`, README row.
- `--qvalue-top-n` in `metabolic_target_prediction.py`; `_qsigtop` launchers (mc, tf, expr, with `_shuff`),
  MetChange `_qsigtop_sdscale`,
  wrapper `submit_train_w_shuff_noGlobals_lr1e-3_min15_qsigtop.sh`, registered, README row.
- Report 05 (`reports/2026-11-17_signature_selection_scale.md`).
- Report 05, section 5 (corrected MetChange results); `results/.../metabolic_coverage_fixedlabels/`.
- `metabolic_target_prediction.py`: `--signature-scale {none,sd,mad}`, `--signature-scale-floor`;
  launchers `launch_predict_targets_A549_mc_noGlobals_lr1e-3_min15_sdscale{,_shuff}.sh`, wrapper
  `submit_train_w_shuff_mc_noGlobals_lr1e-3_min15_sdscale.sh`, variant registered, README row.
- This `STATUS.md` update.

## 9. Working notes

- A safety classifier has repeatedly truncated long prose answers in this project (three times
  on 2026-11-13). Code, tables, numbers and tool output were never affected; the pattern was
  broad written assessments. Keeping write-ups in `reports/` and chat answers short avoids it.
  Starting a new session helps because it drops accumulated context — hence this handoff.
- Laptop has `python3` (no `python`), no torch, bash 3.2 (no `declare -A`), and no network
  data, so only unit-level tests and result analysis are possible there.
