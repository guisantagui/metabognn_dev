# Report 12 — Why the MetChange label matrix is one-dimensional, and what to do about it

Follows report 11 §8, which established that the MetChange label matrix carries ~1.1 effective dimensions and
left the cause unknown. This report identifies the cause, quantifies it, and concludes that the remedy is a
post-hoc removal rather than a change to the penalty.

Everything here is from the diagnostic notebook `scripts/utilities/check_mmp_down_tca_sign.ipynb`, run on the
cluster against HepG2. The MMP/TCA *direction* question that the same notebook settled is written up in
report 11 §10.3b–10.3c and is not repeated.

## 1. Verdict first

**The collapse happens at the LP, and it is an averaging effect rather than information destruction.** Each
metabolite's cost is a penalty-weighted sum over many reactions. Summing averages down the independent
per-gene variation while preserving anything common to all routes, so a small shared component in the
penalties becomes the dominant component of the output.

**The shared component is the sample's overall metabolic cost level**, and 68% of it traces to one scalar: the
average rank of metabolic genes within the whole transcriptome. That is a consequence of ranking over all
~20,587 genes and then using only the ~2,783 GEM genes.

**It is not an artefact and the ranking choice was correct.** It measures how well-positioned metabolism is in
that cell's transcriptional output, which is real. The problem is dominance, not validity.

**The information survives underneath.** Removing one component takes the participation ratio from 1.10 to
**6.36**, with ~12 eigenvalues for 90% of variance and median pairwise r falling from +0.974 to −0.058.

**So: remove the component and supply it as its own feature. Do not re-rank the penalties.**

## 2. The chain test: where the collapse happens

Report 11 §9.2 attributed the collapse to the LP on the strength of the producer-average arm (5.18 dimensions
against MetChange's 1.10). That comparison was invalid, because the producer-average is built from expression
log2FC and never passes through the penalty — it eliminated metabolite-level aggregation, not the penalty. The
claim was retracted. This is the comparison that actually tests it: the three matrices of the real chain, same
samples throughout, absolute values.

| step | shape | PC1 | PR | eig@90% | median r |
|---|---|---|---|---|---|
| 1 absolute expression (genes) | 20,131 × 5,482 | 20.2% | **10.69** | 38 | +0.0149 |
| 2 rank penalties (same genes, GEM subset) | 2,783 × 5,482 | 21.7% | **9.70** | 35 | +0.0084 |
| 3 MetChange scores (metabolites) | 5,700 × 5,482 | **95.2%** | **1.10** | 1 | **+0.9733** |

Expression → penalty loses essentially nothing (10.69 → 9.70, factor 1.1). Penalty → LP destroys it
(9.70 → 1.10, factor 8.8). **The rank penalty transform is not the problem**, which also disposes of the
earlier worry that per-sample ranking might strip between-feature independence.

So report 11 §9.2's conclusion was right and its evidence was not. Both facts are now on record.

## 3. The response is smooth, not switching

A hypothesis worth killing: that the LP is insensitive because a minimum is piecewise constant, so small
expression changes never move the argmin and the cost only jumps when a route switches. That predicts a large
mass at exactly zero plus fat tails.

Measured, over the TCA endpoints and 605 MMP-down perturbations, median basal 10,455:

| | |
|---|---|
| exactly zero | **0.000%** |
| \|log2FC\| < 0.001 | 0.02% |
| below 0.01% of basal | 3.68% |
| **below 1% of basal** | **52.89%** |
| **above 1% of basal** | **47.11%** |
| kurtosis / skew | 40.6 / 3.70 |

**Nothing sits at zero.** 47% of values exceed 1% of basal and the modal fractional change is ~1%. Everything
moves, continuously, all the time.

The hypothesis was mathematically wrong, not merely unsupported: the minimum over routes of a *sum* of
penalties is **continuous and piecewise-linear**, not piecewise constant. A switching argmin produces a kink
in the derivative, not a jump in the value, and between switches the chosen route's cost still varies
continuously with its own penalties. A smooth response was always the correct expectation.

The heavy tails (kurtosis 40.6) therefore come from variation **across perturbations** — some drugs are far
more transcriptionally disruptive than others — not from route switching.

## 4. There is no fixed bottleneck, and movement is uniform

The other candidate: a small set of metabolites that always moves, dragging the rest. Tested by ranking
metabolites by |log2FC| within each perturbation and asking how often each one lands in that perturbation's
top 5%.

| | |
|---|---|
| metabolites in the top 5% for >90% of perturbations | **0** |
| in [0.50, 0.90) | 2 |
| with frequency < 0.01 | 53.3% |
| Jaccard of top-5% sets between random perturbation pairs | **0.065** (chance 0.026) |
| top 1% of metabolites carry | **1.7%** of all movement |
| top 50% carry | 67.6% |

No fixed core. The top movers rotate between drugs, and movement is almost perfectly **uniform** across
metabolites — the top 1% carrying 1.7% is near-proportional.

That also means the ranking test was the wrong instrument for the question it was built for: if everything
moves together, nothing stands out by magnitude, and the shared component shows up only in the correlation
(+0.9733) and in the 93% of each metabolite's variance explained by the global mean.

## 5. What PC1 is

An SVD of the standardised absolute matrix:

* PC1 explains **94.8%** of the subsampled matrix.
* **r(PC1, g) = −1.0000**, where `g` is the mean standardised score across all features — one number per
  sample. The sign is SVD orientation.
* Loadings are **flat**: median |loading| 0.0228, 10th percentile 0.0220, 90th 0.0229, against
  1/√2000 = 0.02236 for a perfectly uniform eigenvector. The top-25-by-loading list differs in the fourth
  decimal and reshuffles with the random subsample, so **it identifies nothing**.

Correlations with individual metabolites are therefore uninformative: ATP −0.989, NADPH −0.996, pyruvate
−0.978, glucose-6-phosphate −0.990 — but with random pairs at +0.97, *every* metabolite correlates ~0.99 with
PC1. This does **not** support attributing PC1 to particular supply metabolites; uniform loadings rule that
out.

The informative entries are the exceptions: **Pi[c] −0.207, H2O[c] −0.277, H2O[g] −0.409**, and a cluster at
exactly **−0.9008** (O2[e], Pi[e], H2O[e]) which is the medium-import degeneracy of §9.

So the matrix is essentially rank-1:

    score(m, s)  ≈  a(m)  +  b(s)

a metabolite-specific baseline plus **one sample-level global cost**. Knowing one number per sample tells you
~95% of what all 5,700 metabolite scores would.

## 6. Where the sample-level cost comes from

This was the sharpest question raised: if penalties are ranks, every sample has an identical penalty
distribution by construction, and the medium is identical — so what varies?

`create_penalties.py` answers it, in a comment:

> *"The matrix is read in row chunks into one float32 array (**all genes are kept**: the rank penalty ranks
> each sample over **every** gene). **Only the GEM genes are written out**."*

```python
n_genes = X.shape[0]                 # ALL ~20,587 genes
ranks   = rankdata(x, axis=0)        # ranked over all of them
pvals   = 1.0 - ranks / n_genes
pen[:, ...] = pvals[gem_rows]        # then only the ~2,783 GEM genes kept
```

The distribution is fixed over 20,587 genes but **not over the 2,783 the LP uses**. What varies is where
metabolic genes sit *relative to the whole transcriptome*. Measured:

| | |
|---|---|
| mean GEM penalty, median across samples | **0.3761** |
| SD across samples | 0.00665 (~1.8% relative) |
| range | 0.3459 – 0.4194 |
| **r(mean GEM penalty, g)** | **+0.8261** |
| Spearman | +0.7782 |
| **r²** | **68.3%** |

Since `g` is 95% of the matrix variance, **~65% of all variance in the MetChange output is that one scalar**.
The remaining ~32% of `g` is reassignment — which specific genes became cheap rather than the overall level.

Two corroborations already on record: the mean GEM penalty shifted by **+0.00108** between MMP-down and DMSO
samples, with **1,539 genes worse against 1,244 better** — an imbalance that a genuinely fixed distribution
could not produce. And the median of 0.376 rather than 0.5 is the direct signature of the subsetting.

## 7. Why a 1.8% penalty shift dominates the output

**The LP is an averaging operator.** Each metabolite's cost sums penalties over many reactions, weighted by
flux. Averaging treats the two kinds of variation oppositely:

* **independent** per-gene variation averages down by ~√n;
* **common** variation survives undiminished and is multiplied by route length.

So a weak shared component in the penalties becomes dominant in the output. The penalty matrix has PR 9.70 —
~10 effective dimensions, only a small shared part — and after summing ~50–100 reactions per route that
shared part is 95% of the result. The magnitudes agree: fluxes are large (LP1 fixes the demand at maximum
production, measured at 19–998 mmol/gDW/h for amino acids), so a 0.001–0.05 penalty shift becomes tens to
hundreds of cost units against a basal of 10⁴–10⁵, which is the ±100–900 range observed.

**1.10 dimensions out of 9.70 is therefore not information destruction. It is an averaging operator
concentrating a small common term.** That is a materially better diagnosis than "the LP destroys the
information", and it is why §8's remedy works.

## 8. The remedy: remove the component, do not re-rank

Removing the leading component, both matrices, `k` = 0…3:

| | log2FC (5,704 × 15,032) | absolute (5,700 × 5,482) |
|---|---|---|
| raw PR | 1.10 | 1.10 |
| **minus 1 component** | **6.36** | **5.24** |
| minus 2 | 6.08 | 5.47 |
| minus 3 | 6.17 | 5.52 |

log2FC: PC1 95.3% → 31.0%, eig@90% 1 → 12, median r +0.9739 → −0.0577.

**It plateaus immediately.** A second and third component gain nothing and even dip slightly, so this is
**one** component, not layered structure. The FNN recipe is therefore "subtract one axis and supply it as one
feature", not "project out the leading few". (The small negative median r afterwards is expected: residuals
from a common component sum to roughly zero across features.)

This also reproduces report 11 §9.2's figure of 6.21 for the same matrix, computed there by removing `g`
rather than exact PC1 — the r = −1.0000 result shows those are the same operation.

### 8.1 Re-ranking within the GEM genes: not the first move, but left open

The mechanism of §6 suggests a one-line change: `rankdata(x[gem_rows])` instead of `rankdata(x)`, which would
fix the mean penalty at ~0.5 per sample by construction and so remove the 68% level-shift component at
source. **It should not be the first move**, for three reasons — but it is worth testing, and §8.1.1 says why.

**It would destroy real information.** Enzymatic genes are genuinely better expressed than, say,
transcription factors, and the median GEM penalty of **0.376 against 0.5** is that fact measured. Ranking
within the GEM genes forces the distribution to be uniform over metabolic genes, so half of them receive
penalty > 0.5 *even when every one is highly expressed*, and a cell with uniformly excellent metabolic
transcription scores identically to one with uniformly poor transcription. The all-genes ranking is the
design that preserves the standing of metabolism within the transcriptome, and that was the stated reason for
choosing it.

**It is less complete.** It removes the level shift — 68% of the axis — and leaves the ~32% arising from
reassignment correlation. Component removal takes all of it.

**It is far more expensive.** Re-ranking means recomputing MetChange for every sample and cell line, a full
Gurobi campaign, against a matrix operation on output already in hand.

And the two are not equivalent in any case. Component removal is a **linear post-hoc operation on the
output**; re-ranking changes the **input non-linearly**, re-spreading the whole penalty distribution and
therefore changing which route is cheapest and what the LP returns. The result would not be "the same minus
PC1" but a different model that also lacks the metabolic-standing information.

**Consistency argument, also valid:** the mixture penalty in `create_penalties.py` is fitted on the pooled
values of all genes. Ranking within the GEM subset would make the two penalty modes non-comparable.

#### 8.1.1 Why it is still worth testing

The three arguments above say what re-ranking would *lose* and what it would *cost*. None of them shows it
would fail to produce a better signature, and there is a real counter-argument the other way: for a
**metabolic** signature, making every metabolic gene's penalty depend on non-metabolic transcription couples
in a dependency that may be the wrong biology for this task, however real the quantity is in itself. A drug
that shuts down ribosomal genes would lower every metabolic penalty and read as "metabolism got easier"
without a single metabolic gene changing. Whether that is signal or nuisance depends on what the signature is
for, and the project's answer — a repurposing tool for metabolic reprogramming — arguably wants it treated as
nuisance.

So it stays open as a **sensitivity test**, scoped to one cell line rather than a campaign: one-line change,
one HepG2 MetChange run, then `feature_dimensionality.py`. Prediction, so it is falsifiable: PC1 falls and PR
rises, magnitude unknown — anywhere from ~2 to ~6 is consistent, and if it stays near 1.10 then the level
shift was not the binding mechanism and the remaining 32% (reassignment correlation) is itself enough to
collapse the output. Two things to check alongside: compare the two versions on a **downstream task** rather
than on dimensionality alone, since the GEM-only version is expected to lose real information; and confirm
the `--threshold 2` ceiling does not introduce its own artefact, because it would then force a different ~56
genes to penalty 1.

### 8.2 What to do instead

1. Keep `rank` over all genes.
2. Remove the leading component from the label matrix before any feature selection, and supply it as **one
   additional feature** — it is a real quantity (overall metabolic priority) and should not be discarded, only
   separated.
3. Report 05 found that top-100 feature selection is driven by metabolite scale. Scale is this axis, so the
   current selection is choosing on the common component. Removing it first changes what gets selected, and
   that is the concrete place this finding bites.

## 9. Two incidental findings

**A set of metabolites is scored by import alone and carries no information.** The ten least-moving features
are beta-carotene, aquacob(III)alamin, cobamide-coenzyme, alpha-tocopherol, folate, thiamin, cystine, lipoic
acid and TAG-VLDL pool — **nine of them medium components** from the 60 exchanges in
`Human2-GEM/rpmi_medium.json` — sharing `mean_abs = 0.044295` to six significant figures with
`frac_perts_moved = 0.0`. Identical values across distinct metabolites means an identical LP solution: direct
import at the median penalty, expression-independent. The O2[e]/Pi[e]/H2O[e] cluster at exactly −0.9008 is the
same degeneracy.

Report 11 §8.5 tested for exactly this saturation and concluded it was not happening. That test used SD of
log2FC across perturbations, which is small but non-zero for these features; the tell is *identical values
across different metabolites*, which it did not look for. The conclusion there was wrong for the
medium-supplied set, and report 11 §8.5's folate finding was the visible tip of it.

**The top movers are coherent biology.** Margaric, pentadecylic, tridecylic, heneicosanoic, nonadecylic,
tricosanoic and heptadecenoic acids, plus Mevalonate-Lactone — long-chain fatty acids and a mevalonate
intermediate. That matches the subsystem reallocation of report 11 §10.3b exactly (fatty-acid elongation
+0.0093, cholesterol biosynthesis +0.0061, steroid metabolism +0.0109): lipid and sterol synthesis suppressed
under mitochondrial stress. So the metabolite-specific residual is tracking something real.

## 10. What this means for the FNN

**MetChange residual PR 6.36 against expression's 13.96.** MetChange contributes ~6 dimensions of
metabolite-specific information plus one global axis, where expression brings ~14 at gene level. A real
representation, not a degenerate one — but roughly half the dimensionality of the input it was derived from.

That makes the decision question concrete: **are those ~6 dimensions complementary to expression's ~14, or a
lossy summary of them?** Directly measurable — canonical correlation between the two residual matrices, or
whether MetChange-residual improves target prediction over expression alone — and it needs no new MetChange
run.

It also softens the case for abandoning the metabolite-centric approach. Report 11's framing was that a
reaction-centric method might retain dimensionality where MetChange does not. That remains worth testing, and
STATUS now counts reaction-centric reprogramming as in scope. But the collapse is now understood, is largely
reversible, and leaves ~12-dimensional structure underneath — so it is no longer a reason on its own to switch
methods.

## 11. Corrections made in reaching this

Recorded because several were asserted in reports before being tested.

| claim | status |
|---|---|
| "the LP is what collapses the dimensionality" (report 11 §9.2) | **right, but the evidence given was invalid**; retracted, then confirmed by §2 here |
| "a minimum is piecewise constant, so expect a spike at zero" | **wrong**: min of sums is continuous piecewise-linear; §3 |
| "rank-relativity explains the TCA direction" | **wrong**: corr(d_expr, d_penalty) = −0.924, 91% of genes concordant |
| "PC1 is attributable to shared supply metabolites (ATP, NADPH…)" | **not supported**: loadings uniform; the high correlations are what any metabolite gives |
| "re-rank within the GEM genes" | **demoted, not withdrawn**: not the first move (destroys real information, removes only 68%, costs a run), but a live sensitivity test; §8.1, §8.1.1 |
| "changing the penalty cannot help" | **imprecise**: the *transform* cannot (§2), but the *gene universe* the ranks are taken over drives 68% of the axis (§6). Different things, conflated |
| medium-supplied metabolites are not saturated (report 11 §8.5) | **wrong for that set**; §9 |

## 12. Open

1. **Is the residual complementary to expression, or a lossy summary?** Canonical correlation or a
   target-prediction comparison. The FNN decision rests on this and nothing else is blocking it.
2. **Confirm the averaging mechanism** by extracting the flux-carrying reactions per optimal solution and
   measuring route overlap across metabolites. Needs Gurobi. The same run would confirm min-over-routes for
   report 11 §10.3b.
3. **Decide what to do with the dead medium-supplied features** — nine or more carry no information and
   should probably be dropped from every set before scoring.
4. **Re-run the four-arm comparison with `--mode sign`**, now that `reproduce_report11.sh` includes it, to get
   the absolute-direction table for all four representations.
5. The NCI-60 run, where the absolute/cross-cell-line case for a mixture penalty is genuinely different from
   the LINCS case and should be run both ways.

## 13. Reproducing

| result | how |
|---|---|
| §2 chain test | notebook, "WHERE DOES THE DIMENSIONAL COLLAPSE HAPPEN" cell |
| §3 movement distribution | notebook, histogram cell |
| §4 bottleneck / uniformity | notebook, top-movers cell |
| §5 PC1 identity and loadings | notebook, PC1 correlation cell |
| §6 mean GEM penalty vs g | notebook, two-line check |
| §8 the k-sweep | notebook, `residual_pr()` on both matrices |
| the four-arm comparison | `bash scripts/mc_phenotype_assoc/reproduce_report11.sh`, or `sbatch launch_report11_slurm.sh` |

The notebook writes its artefacts and a `run_info.json` to
`<labels_base>/results/mc_phenotype_assoc/HEPG2/tca_sign_diagnostic/`. **The figures in this report were read
from terminal output, not from those files** — the save cell had not been run when they were produced, so
re-running the notebook end to end is what puts them on disk.
