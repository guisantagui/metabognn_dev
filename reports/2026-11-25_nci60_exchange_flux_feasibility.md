# Report 11 — NCI-60 exchange fluxes: feasibility check before computing MetChange

*2026-11-25. Assessment of `data/validation/human1_nci60/` before committing to a MetChange run. No analysis
has been done; this establishes what the data is, what is reachable, and what will go wrong.*

**Summary.** Feasible, and it would be the first dataset in the project that measures the quantity MetChange
actually claims. All 26 exchange metabolites are reachable in the MetChange feature space, all 11 cell lines
have RNA-seq, and that RNA-seq covers **98.5% of GEM genes - slightly better than the LINCS data the existing
labels are built on**. Three things to handle: the data is 26 metabolites and not 111, the flux values span
seven orders of magnitude so an across-metabolite correlation would be meaningless, and CellMiner omits all
13 mtDNA-encoded OXPHOS subunits, which breaks the GPR of five complex reactions.

## 1. What the data is

`NIHMS1590510-supplement-Data_File_S5.xlsx` (Robinson et al., Human1-GEM, *Sci Signal* 2020,
doi 10.1126/scisignal.aaz1482), 30 rows x 13 columns:

- **11 cell lines**: HT29, HOP-62, 786-O, MALME-3M, UO-31, HS 578T, HOP-92, NCI-H226, MDA-MB-231/ATCC,
  RPMI 8226, SR.
- **26 metabolite exchange reactions**, already mapped to Human-GEM reaction ids (`HMR_90xx`) - so the
  structure-mapping work that sank Hallinger and Xu is already done.
- Plus two rows that are **not** metabolite exchanges: `biomass_human` (growth rate) and
  `prot_pool_exchange` (Ptot, the enzyme-constraint protein pool from the GECKO formulation).

**It is 26 metabolites, not 111.** Sign convention is release-positive: L-lactate +0.77, glucose -0.59,
glutamine -0.12, i.e. the expected Warburg pattern.

On "only 11 had reliable nutrient exchange rates" - the likely reason, **stated as interpretation and not yet
verified in the Human1 methods**: converting the CORE consumption/release measurements (Jain et al., *Science*
2012) into absolute flux units usable as model constraints needs a measured growth rate and a dry-weight-per-cell
figure per line, and only 11 of the 60 had both. The presence of a `biomass_human` row per cell line is
consistent with that. Worth confirming before it goes in a methods section.

## 2. The 26 metabolites, and the crossed prediction they allow

| direction | metabolites |
|---|---|
| **released** (8) | L-lactate, glutamate, alanine, ornithine, proline, glycine, AKG, pyruvate |
| **consumed** (18) | glucose, glutamine, serine, lysine, leucine, isoleucine, valine, threonine, tyrosine, phenylalanine, methionine, arginine, asparagine, aspartate, tryptophan, choline, L-carnitine, folate |

**16 of the 26 are amino acids.** That makes this largely a test of the amino-acid arm - the arm that failed
in report 10, section 13 (mean r -0.06 to -0.10). Either that failure was an artefact of a contaminated drug
group, which 13.4 showed was about 95% non-amino-acid drugs, or amino-acid capacity genuinely is not
predicted. **This dataset distinguishes those two, which is its main value, and it should be entered expecting
a test of the weaker arm rather than a replication of the nucleotide result.**

The sign split gives a **crossed prediction**, the same structure as the Panorama GCP rationale and much harder
to pass by accident than a single correlation:

- for a **released** metabolite, higher predicted production capacity should mean **more** release;
- for a **consumed** metabolite the cell is not producing it, so higher internal capacity should mean **less**
  uptake - the opposite sign on the same readout.

## 3. Reachability: all 26, in the extracellular compartment

Every one of the 26 has at least one species in the MetChange feature space (5,707 metabolites), and **all 26
have the extracellular (`e`) species** - the compartment that maps directly onto an exchange flux. The hub
flags sit on the cytosolic, mitochondrial and ER species, **not** on `e`, so the extracellular species are not
hub-contaminated. 16 of the 26 are flagged as hubs in `c`, which matters for the `nohub` set variants but not
for this analysis.

## 4. Expression: all 11 lines, and better GEM coverage than LINCS

CellMiner RNA-seq composite (`nci60_RNA__RNA_seq_composite_expression.zip`), 60 cell lines x 24,911 rows.

| | value |
|---|---|
| S5 cell lines with CellMiner expression | **11 of 11** |
| distinct gene symbols | 23,808 (24,911 rows, **1,103 duplicated symbols**) |
| symbols mapping to Ensembl via `hgnc_complete_set.txt` | 22,013 -> 21,950 distinct ENSG |
| **GEM genes covered** | **2,485 of 2,522 (98.5%)** |
| for comparison, GEM genes covered by the LINCS gene space | 2,470 of 2,522 (97.9%) |

*The figure was first computed as 98.41% by calling `symbol_to_ensembl_map(hgnc, {})` with empty overrides,
which is not the project's convention. Every real consumer passes `gem_symbol_overrides(gem)`; doing so
recovers three more genes (98.53%). The overrides matter for exactly four symbols - **SOD2, LRTOMT, PPIAL4C,
SLC37A4** - where HGNC and Human-GEM assign different Ensembl ids. Small here, but SOD2 is ROS
detoxification, so the same slip in a mitochondrial analysis would not be harmless. See the conventions note
in `scripts/mc_phenotype_assoc/REPRODUCING.md`.*

Two practical notes:

- **786-O is spelled `RE:786-0`** in CellMiner, with a zero rather than the letter O. Normalised name matching
  silently drops it. Hardcode the mapping.
- **1,103 duplicated symbols must be aggregated.** The rank-based MetChange is invariant to monotone
  transforms, so the log scale of the CellMiner values is not itself a problem - but *summing* duplicates is
  only valid in linear space, which is the reason `prepare_input_matrices.py` converts. Either convert before
  summing, or aggregate with a rank-preserving function such as max.

So the input is **equal or better than what the existing labels use**: this would run MetChange on direct
RNA-seq rather than on L1000-inferred expression.

## 5. The one real defect: CellMiner has no mtDNA-encoded genes

Of the 40 GEM genes with no CellMiner expression, **13 are the mtDNA-encoded OXPHOS subunits** - ND1-6,
ND4L, COI-III, ATP6/8 (`Su6m`, `URFA6L`) and `UQCR3`. Most of the remaining 27 are USP17L and CDY paralogs in
the `Isolated` subsystem, plus QARS1, GGTLC3, TL-PST and SDR42E2.

**This is a new gap, not a pre-existing one: the LINCS gene space contains all 13.**

It touches **5 GEM reactions**, and they are the worst five it could touch - all `Oxidative phosphorylation`,
all AND-dominated complex GPRs:

| reaction | genes in GPR | ANDs | ORs | mtDNA-encoded and missing |
|---|---|---|---|---|
| MAR06914 | 19 | 14 | 4 | 3 |
| MAR06916 | 21 | 20 | 0 | 2 |
| MAR06918 | 10 | 9 | 0 | 1 |
| MAR06921 | 45 | 44 | 0 | 7 |
| MAR13081 | 19 | 14 | 4 | 3 |

With AND logic a complex needs all its subunits, so whichever way `create_penalties.py` treats an absent gene
determines the penalty for complexes I, III, IV and V. Five reactions of 12,931 is a small footprint, but
OXPHOS constrains ATP supply for the whole network, so a mis-set penalty there can propagate into every
prediction.

**Options, in order of preference:** impute the 13 at the median of the nuclear-encoded subunits of the same
complex (defensible - mtDNA transcripts are abundant and rarely the limiting subunit); or drop them from the
GPRs; or source expression for these lines from a dataset that includes `MT-` genes. Whichever is chosen, run
the analysis **both ways** and report the sensitivity, because this is the one input difference between the
NCI-60 run and every existing label set.

## 6. Design, fixed before running

1. **Do not correlate across metabolites within a cell line.** The fluxes span **seven orders of magnitude**
   (L-lactate +1.4, glucose -0.9, pyruvate 2.5e-07, L-carnitine -1.2e-05), so such a correlation would be
   almost entirely the magnitude ordering of three metabolites. This is the metabolite-identity confounder of
   report 10, section 8, in a far more extreme form.
2. **Use the between-cell differential.** Per metabolite, (flux in X - flux in Y) against
   (MetChange in X - MetChange in Y). **11 cells give 55 pairs**, against PRISM's 10, so the saturated
   Wilcoxon of report 10, 11.6 stops being the binding constraint (floor 1/2^55).
3. **Check per-cell coupling as well**, since report 10, 8.3 established that differencing cancels
   item-level but not cell-level confounding.
4. **Split by direction and predict opposite signs**, per section 2 above. The released and consumed arms are
   the test and its internal control.
5. **The comparator, corrected.** An earlier version of this section argued the baseline had to be an
   expression-constrained FBA - E-Flux, GIMME, iMAT, or the ec-GEM FVA of the Human2 paper at ~81% flux
   consistency. **That is the wrong comparator for this project's question and the recommendation is
   withdrawn.** Benchmarking MetChange against FBA variants asks whether MetChange is a good flux-prediction
   method, which is the MetChange paper's own question and already published. The question here is whether the
   integration tracks real metabolic change well enough to **justify it as the metabolic signature fed to the
   FNN** - so the comparator must be whatever would otherwise be fed to the FNN. That is the transcriptional
   signature. E-Flux, GIMME and iMAT also fail on output space: they are **reaction-indexed**, while the
   bipartite network propagates from **metabolite** nodes, so their predictions are not substitutable inputs.

6. **The right comparator is a metabolite-centric expression null**: for each metabolite, the mean expression
   of the enzymes catalysing reactions that **produce** it. This is what MetChange reduces to if the LP and the
   network-wide stoichiometric coupling add nothing over local gene averaging. It is sharper than the
   comparator used so far, in three ways:

   - **same output space** - one value per metabolite, so it drops straight into the differential design;
   - **tighter than the current expression arm**, which uses *subsystem* GPR gene sets (`nucleotide gene`
     covers a 192-metabolite subsystem). A per-metabolite producer average is the transcriptional information
     specific to *that metabolite*, which is a far more adversarial null;
   - **it isolates the claim.** "Expression integrated with stoichiometry" versus "expression without
     stoichiometry, at the same resolution".

   The machinery exists: `build_target_relevance.py` already derives producing routes from the stoichiometry,
   so the hop-1 producers of each metabolite are a by-product of code already written.

   This also retro-sharpens report 10, section 12: the 104% retention there was measured against a subsystem
   gene set that was itself flat (-0.008). Against a per-metabolite producer average the test would be much
   harder, and correspondingly more convincing if it survives.

## 7. Verdict

Worth doing, and the strongest of the three candidates discussed (the others being the To et al. 2019 K562
CRISPR screen, which deposited **no expression data** - verified, GEO returns zero datasets for PMID 31730859 -
and Kang et al. 2024, which is contact-inhibition quiescence in **mouse** fibroblasts with a single
proliferating-versus-quiescent contrast).

It is the first dataset here that measures **flux**, which report 05, 9.3 identified as the ground truth
MetChange's claim actually requires and which metabolite pools cannot supply. Expect it to be a hard test:
26 metabolites, 16 of them amino acids, and the amino-acid arm has failed once already.

**Amended by section 8, added later the same day.** Two changes to this verdict. The design is **25
metabolites, not 26** - folate is a dead feature (8.5). And the priority has moved: section 8 found the
MetChange label matrix to be effectively one-dimensional, which makes the producer-average null of section 6
the gating analysis. Run that before this dataset, because it determines whether NCI-60 would be confirming
a result or chasing one.

---

## 8. The MetChange label matrix is one-dimensional (2026-11-25, added after section 7)

This section came out of a question about the growth medium and ended somewhere more consequential. It
changes how every set-level result in reports 09 and 10 should be read, so it is recorded here rather than
in a separate file.

Reproduce with `python3 -I scripts/mc_phenotype_assoc/feature_dimensionality.py`.

### 8.1 The finding

Eigen-spectrum of the feature x feature correlation matrix of each `log2FC` label matrix:

| matrix | shape | PC1 | participation ratio | eigenvalues for 90% var | median r, random pair |
|---|---|---|---|---|---|
| HEPG2 MetChange | 5,707 x 14,632 | **95.5%** | **1.10** | **1** | **0.9748** |
| A549 MetChange | 5,707 x 44,394 | **95.0%** | **1.11** | **1** | **0.9734** |
| HEPG2 expression | 8,000 x 15,032 | 19.6% | 13.96 | 58 | 0.0655 |
| A549 expression | 8,000 x 7,633 | 15.0% | 17.83 | 52 | 0.0445 |
| iid Gaussian (calibration) | 1,500 x 5,000 | 0.2% | 1154 | 1111 | 0.0000 |

A single metabolite's variance is **93.0%** explained by the global component alone (the mean standardised
score across all 5,707 metabolites). 76.9% of randomly chosen metabolite pairs correlate above r = 0.95.

Two things make this a statement about MetChange rather than about the data:

- **It is not what the input looks like.** Same cell lines, same perturbations, same limma log2FC pipeline:
  expression needs 52-58 eigenvalues for 90% of variance and random gene pairs correlate at 0.045-0.066.
  The collapse from ~15 effective dimensions to ~1.1 happens in the expression-to-MetChange step.
- **It replicates across cell lines.** HEPG2 and A549 agree to two decimal places on the participation
  ratio despite a 3x difference in perturbation count.

### 8.2 A correction to a figure first reported in this session

The first pass at this analysis reported **9.55** effective dimensions. That was wrong. The matrix had been
centred across features (`Z - Z.mean(0)`), which subtracts each perturbation's mean over all metabolites -
and when the features share a dominant component, that mean *is* the shared component. PC1 was therefore
measured on data from which the effect had already been removed, while the pairwise correlations in the same
table were measured without that step. The two numbers in that table (PC1 20.6% beside median r 0.97) are
arithmetically incompatible, which is what exposed the error.

The true figure is more extreme, not less. `feature_dimensionality.py` works only from the correlation matrix
and carries a guard (`_check_no_double_centring`) that refuses a matrix whose shared component is already gone.

### 8.3 What the global component is NOT

The obvious candidate is response magnitude - a drug that perturbs more genes moves everything. It is not that:

| | Spearman r with the global MetChange component |
|---|---|
| mean abs(expression log2FC) per perturbation | **+0.141** |
| fraction of genes with abs(log2FC) > 1 | **+0.148** |

Both are weak, so the `n_sig` covariate matching used throughout reports 08-10 - a response-magnitude proxy -
does not target this component.

**But the reported statistics are partially protected by their own definition, which a first version of this
section got wrong.** `mc_phenotype_assoc.py` computes `z` as *subtract the column mean, divide by the column
SD* (line 346, applied at line 400), where the column is one perturbation across all features. That is a
within-perturbation centring, and the global component is close to that column mean. `auroc` is a
within-column ECDF and so invariant to any monotone column transform; `sign_dev`, `sig_dev` and
`sig_signed_dev` each subtract their own column mean too. Only `sign_mean` is uncentred.

The protection is partial, not total, because `z` subtracts a **uniform** column mean while each feature has
its own loading on the global component (median 0.992, range -0.240 to 0.998). Measured directly - set-level
`z` computed on the raw matrix versus the fully residualised one:

| set | r(z_raw, z_residual) | r(auroc_raw, auroc_residual) |
|---|---|---|
| `glycolysis/cyto_nohub` | 0.846 | 0.845 |
| `dntp_endpoints/all` | 0.841 | 0.767 |
| `tca_endpoints/all_nohub` | 0.800 | 0.771 |
| `amino_acid_endpoints/all` | 0.360 | 0.350 |
| `oxphos/mito` | **0.259** | **0.423** |

So the existing numbers are neither invalid nor safe: the glycolysis and nucleotide statistics are largely
unaffected by full residualisation, while `oxphos/mito` and the amino-acid sets change substantially. Section
10 re-runs the designs on the residual matrix to see which conclusions move.

### 8.4 What it does not mean: the curated sets still carry set-specific signal

The wrong test is "remove PC1 and see what survives". With PC1 at 95%, that removes the dataset; and requiring
a real signal to be orthogonal to the dominant direction is backwards, since biology aligned with that
direction would be discarded along with it. The right question is **specificity**: is the set's non-global part
larger than a size-matched random set's?

| set | n | r(set, global) | residual SD | percentile among 2,000 random sets |
|---|---|---|---|---|
| `dntp_endpoints/all` | 12 | 0.990 | 0.1428 | **99.7%** |
| `ntp_endpoints/all` | 19 | 0.988 | 0.1520 | **100.0%** |
| `nmp_endpoints/all` | 37 | 0.989 | 0.1498 | **100.0%** |
| `tca_endpoints/all_nohub` | 20 | 0.983 | 0.1845 | **100.0%** |
| `amino_acid_endpoints/all` | 88 | 0.998 | 0.0696 | **100.0%** |
| `oxphos/mito` | 13 | 0.992 | 0.1290 | 99.2% |
| `glycolysis/cyto_nohub` | 27 | 0.994 | 0.1053 | 100.0% |
| `etc_endpoints/all_nohub` | 6 | 0.985 | 0.1705 | 94.8% |

Every curated set carries more set-specific movement than arbitrary metabolite collections of the same size.
The mechanism is that a random draw averages over unrelated metabolites and collapses toward the global mean,
while a functionally coherent set has correlated deviations that survive averaging. The test is correctly
size-adjusted: `amino_acid_endpoints` has the smallest absolute residual (0.0696) and the highest percentile,
because 88-member random sets average out further.

So `dntp_endpoints` is not the global component wearing a label. **What this does not establish** is that the
*phenotype associations* ride on the set-specific part rather than the global part - a different question,
and the one the producer-average null of section 6 would settle.

### 8.5 Two incidental findings from the same investigation

**The medium does not explain the amino-acid null.** The hypothesis was: amino acids are in the medium, so
import is the cheapest route and inhibiting biosynthesis cannot move the score. The mechanics check out -
all 60 medium exchanges have `gene_reaction_rule = (none)`, GPR-less reactions are scored at the sample's
**median** penalty (`run_metchange_exe.py`, `results[rxn_id] = median_arr.copy()`), which with `threshold=2`
is ~0.5, and 19 of 20 amino acids are in `Human2-GEM/rpmi_medium.json`. But it is not the binding constraint:
**alanine, the one amino acid with no exchange at all, behaves like the other 19** (HEPG2 basal 1.388e4 vs
supplied median 1.356e4; SD 354 vs 343). Closing an amino acid's own exchange also leaves it producible,
because the model interconverts from the other 18. The null is explained by 8.1 instead: the scores are
dominated by a global term, so metabolite-specific perturbations move them only slightly.

The 8 essential amino acids are pairwise **r = 1.00000** (min 0.99997), effective dimensions 1.00 - but this
is a symptom of 8.1, not a separate defect, and it does not warrant dropping them from the set: redundancy
among members of a set reported as null cannot manufacture a false positive.

**Folate must come out of the NCI-60 comparison.** SD of log2FC is 0.050 in A549 and 0.057 in HEPG2 against a
median of ~290 across all features - the **0th percentile in both cell lines**, i.e. a dead feature, and
consistently so, which makes it structural rather than a cell-line quirk. The section 2 design is therefore
**25 metabolites, not 26**.

### 8.6 What this means for the FNN, which is the actual purpose

The project's thesis is that MetChange supplies a metabolic signature carrying information transcriptomics
does not. Section 8.1 does not refute that, but it raises the bar and supplies the decision criterion this
project has been missing:

1. **The raw scores are a poor input representation.** Feeding an FNN 5,707 features whose pairwise
   correlation is 0.97 spends the network's capacity on one scalar. The global component should be removed
   (per-perturbation centring across metabolites) and supplied as a single separate feature - it may well be
   biologically real, but it should not be smeared across 5,707 inputs.
2. **The comparison that decides the investment** is whether the *residual*, after that removal, carries
   drug-target-relevant information beyond what expression provides. Expression enters with 14-18 effective
   dimensions; MetChange enters with ~1.1 plus a residual. The residual is where the entire case lives.
3. **This is now testable and was not before.** 8.4 shows the residual is structured along curated pathway
   lines, which is the necessary condition. The sufficient condition is that it predicts.

Nothing here says MetChange is invalid. It says the method compresses hard, that the compression is a
property of the implementation rather than of the data, that the surviving signal is real but small, and
that the project's existing `n_sig` controls do not address it. Every set-level claim in reports 09 and 10
should be read against its random-set column (`rand=`) rather than its parametric p-value, which the
association runs have been computing all along.

### 8.7 Open, in priority order

1. Build the producer-average null (section 6) and apply it to the existing MMP and PRISM results. It is
   now the gating analysis, not an optional comparator.
2. Re-read the reports 09 and 10 headline claims against `rand=` and record which survive.
3. Test whether the global-removed residual predicts, before committing to the FNN.
4. Diagnose *why* the implementation collapses this hard - the rank-based per-sample penalty is the obvious
   suspect, since it makes every sample's penalty distribution identical by construction - and check whether
   a different penalty transform preserves more dimensions.

---

## 9. The producer-average null: does the LP add anything? (2026-11-25)

Reproduce with `bash scripts/mc_phenotype_assoc/reproduce_report11.sh 9 10 11`.

### 9.1 The comparator

For each metabolite, the mean expression log2FC of the genes catalysing reactions that **can produce it**
(forward flux and a positive coefficient, or reverse flux and a negative one; currency metabolites excluded,
reactions with no GPR contribute no genes). Same gene-level input MetChange starts from, aggregated to the
same metabolites, with **no LP, no stoichiometric balancing, no network propagation** - the GEM is used only
to decide which genes to average. Built by `build_derived_labels.py --mode producer_average`.

HEPG2: 4,403 of 5,707 metabolites scored; 1,304 have no measured producing enzyme. Producer sets are small
(median 2 genes, 90th percentile 12, max 203).

This is the metabolite-centric comparator section 6 called for. It is the right one, and E-Flux / GIMME /
iMAT are not, because those score fluxes while the quantity under test is a per-metabolite production
capacity. It is also a *strong* comparator: if MetChange cannot beat it, the LP adds nothing beyond averaging
the relevant enzymes.

### 9.2 The LP is what collapses the dimensionality

| representation | PC1 | participation ratio | median r, random pair |
|---|---|---|---|
| MetChange (LP) | 95.5% | **1.10** | 0.975 |
| producer-average (same genes, same metabolites, **no LP**) | 39.8% | **5.18** | 0.367 |
| global residual | 31.6% | 6.21 | -0.055 |
| expression (genes) | 19.6% | 13.96 | 0.066 |

Aggregating expression to metabolites by producer set does induce correlation - metabolites share producing
enzymes, so r = 0.37 is expected - but it leaves 5.18 effective dimensions. **The LP takes that to 1.10.**
Section 8's collapse is attributable to the LP step specifically, not to metabolite-level aggregation.

## 10. Four arms, one design: the result

Identical phenotype (`attene_ramos_mmp.tsv`), identical sets, identical flags, all four arms recomputed in
the same run rather than compared against numbers quoted from report 09 - the sets file gained four
subsystems and four curated sets after that report was written, so a cross-report comparison could have
reflected the pipeline instead of the arm.

`compare_arms.py` prints only keys present in **every** arm, so the table cannot become "each arm's best
row". n = 144 cases / 475 controls throughout. Random-set p floor 0.004975 (201 draws).

### 10.1 MMP_down, the `z` statistic

| set | variant | MetChange | producer-avg | residual | expression |
|---|---|---|---|---|---|
| tca | cyto_nohub | 0.393 **rand 0.005** | 0.460 rand 0.43 | 0.375 **rand 0.040** | 0.405 rand 0.12 |
| tca | mito_nohub | 0.362 **rand 0.005** | 0.422 rand 0.10 | 0.382 **rand 0.030** | 0.405 rand 0.12 |
| tca | all_nohub | 0.383 **rand 0.005** | 0.432 rand 0.26 | 0.378 **rand 0.025** | 0.405 rand 0.12 |
| pyruvate | all_nohub | 0.389 **rand 0.010** | 0.461 rand 0.41 | 0.379 rand 0.065 | 0.477 rand 0.72 |
| pyruvate | mito_nohub | 0.365 **rand 0.015** | 0.490 rand 0.92 | 0.355 **rand 0.030** | 0.477 rand 0.72 |
| tca | cyto | 0.391 rand 0.020 | 0.457 rand 0.49 | 0.352 **rand 0.010** | 0.405 rand 0.12 |
| tca | all | 0.395 rand 0.020 | 0.431 rand 0.22 | 0.370 **rand 0.010** | 0.405 rand 0.12 |

The rank-based `auroc` statistic gives the same ordering (tca mito_nohub: MetChange 0.360 rand 0.005,
producer-avg 0.427 rand 0.16, residual 0.370 rand 0.015, expression 0.412 rand 0.15), so this is not an
artefact of the heavy tails of limma log2FC.

**Two conclusions, both of which the design was set up to be able to refute.**

1. **MetChange beats the producer-average, decisively.** On TCA and pyruvate the producer-average is null
   (rand 0.10-0.92) where MetChange sits at the random-set floor. Same genes, same metabolites, no LP - no
   signal. **The LP is doing real work**, and this is the first result in the project that isolates the
   integration itself rather than comparing MetChange to gene-level expression.
2. **The signal is not the global component.** The residual arm keeps it: 0.352-0.382 at rand 0.010-0.040.
   Removing the dominant direction costs some power - 0.005 to 0.025-0.040 - so the global component carries
   *part* of the signal, but not the substance of it. Section 8 does not invalidate the report 09 section 4
   result.

### 10.2 The producer-average is not inert, which is what makes 10.1 meaningful

A null comparator proves nothing if it cannot detect anything. It can: in the **MMP_up** contrast (n = 76) it
is the *strongest* arm, with `nucleotide/all` at AUROC 0.259 and `purine/mito` at 0.296, both at the
random-set floor, where MetChange is null there (0.630-0.649, rand 0.085-0.13). So the two representations
succeed on different contrasts, and 10.1 is a genuine comparison rather than a live arm beaten against a dead one.

### 10.3 The sign convention, which inverts between the arms

`run_metchange_exe.py` minimises the penalty sum (`ModelSense = GRB.MINIMIZE`) and its own comment at line
266 states that **"a metchange score of zero means the metabolite can be easily produced"**. The MetChange
label is therefore a **cost**: higher means *harder* to produce. Expression log2FC runs the other way.

So the apparent MMP_up contradiction in 10.2 - producer-average 0.259 against MetChange 0.630 - is
**agreement**: enzyme expression down and production cost up both say nucleotide capacity falls.

### 10.3b The MMP_down direction CONFLICTS with both expression arms, and is unexplained

No sign flip is applied anywhere in the pipeline: the statistics are computed on the raw log2FC matrix. So
AUROC 0.362 means MMP-down drugs have **more negative** TCA log2FC, i.e. a **gain** of consistency with the
transcriptome - the TCA metabolites became *easier* to produce from the drug-treated transcriptome than from
the control's. Put the three arms side by side on the same set, same contrast, same drugs:

| arm | AUROC | median case | median control | rand p |
|---|---|---|---|---|
| MetChange, `tca mito_nohub` (cost) | 0.362 | **-0.013** | +0.186 | **0.005** |
| producer-average, expression of the TCA metabolites' producing enzymes | 0.422 | -0.154 | -0.053 | 0.10 |
| expression, `tca gene` set | 0.405 | -0.059 | +0.077 | 0.12 |

**Both expression arms say TCA genes go DOWN in MMP-down drugs.** Lower expression means higher rank-based
penalties, which should make TCA metabolites *more* expensive and give a POSITIVE log2FC, AUROC > 0.5.
MetChange says the opposite, and it is the only one of the three that is significant.

This is a real inconsistency, not merely a counterintuitive direction, and it is not resolved. What it is
*not*: a global-component artefact (section 10.1's residual arm retains the signal) or a responsiveness
artefact (report 09 section 4's `n_sig` matching retains it). Two candidate mechanisms, each with a specific
diagnostic:

1. **The LP takes a minimum over routes.** The cheapest route to a TCA metabolite need not use TCA-cycle
   enzymes. If a drug suppresses TCA genes while an alternative route's genes improve in rank, the minimum
   can fall even as the canonical route gets worse, and the binding route can differ between samples. The
   producer-average cannot show this, because it averages over ALL producing enzymes rather than taking the
   cheapest subset - which is also why MetChange beating it in 10.1 and disagreeing with it here are the same
   fact seen twice. *Diagnostic:* extract the reactions carrying flux in each metabolite's optimal solution
   and check which genes actually bind. Needs Gurobi, so cluster only.
2. **Rank-relativity.** Penalties are `1 - rank/n` on ABSOLUTE expression within each sample, while the
   expression arm reports log2FC. A gene can fall in log2FC yet RISE in within-sample rank if other genes
   fall further. The two quantities can therefore move in opposite directions by construction.
   *Diagnostic:* for MMP-down drugs against their controls, compute the change in within-sample rank of the
   TCA genes and compare its sign to their log2FC. Needs the absolute expression matrix
   (`data/perturbation/lincs/cell_matrices/<CELL>/expr.csv`): present locally for A549 (10 GB), absent for
   HEPG2.

**Until one of these is settled, the TCA result supports the statistical claim but not a mechanistic one.**
The statistical claim - that MetChange carries information neither expression arm has - is if anything
strengthened by the disagreement, since a proxy for TCA gene expression could not disagree with TCA gene
expression. But what it is tracking instead is currently unknown, and section 11's framing should be read
with that limit in mind.

### 10.4 Caveats that limit how far 10.1 can be pushed

- **Multiple testing.** 82 set x variant keys are scored per statistic. At the random-set floor of 0.005 that
  is ~0.4 expected false positives, so the three `tca` rows at the floor are credible; the rows at
  rand 0.04-0.07 are not individually safe.
- **One cell line, one phenotype.** This is HEPG2 against Attene-Ramos MMP. A549 has both label arms locally
  and the launcher takes `cell=A549`; that run has not been done.
- **The sets are not independent.** Section 8.4 showed effectively ~1 dimension across the whole matrix, and
  `tca`, `pyruvate` and `glycolysis` overlap in membership, so the seven rows above are not seven
  independent confirmations.
- **`--mode value` only.** `qvals.csv` belongs to the original limma fit and does not describe a derived
  matrix; `build_derived_labels.py` symlinks it and writes a `PROVENANCE.txt` saying so, so a later
  `--mode sig` run cannot silently use the wrong q-values.

## 11. Where this leaves the FNN question

Section 8 said the case for MetChange lives entirely in its residual, and that whether the residual predicts
was untested. Section 10 tests it on the project's strongest existing result and the answer is **yes, it
survives** - weakened, but significant, and still beating both expression comparators.

Combined with 9.2, the picture is coherent: the LP compresses the representation hard (1.10 effective
dimensions against the producer-average's 5.18), and it simultaneously adds signal that the producer-average
does not have. Compression and informativeness are not the same axis. That is the first evidence in this
project that the stoichiometric integration contributes something expression alone does not, which is the
claim the whole FNN plan rests on.

What it does **not** settle: whether that signal is useful for drug-target prediction specifically. The MMP
design tests whether the signature separates drugs by a measured mitochondrial phenotype, not whether it
identifies targets. The remaining step is a grouped cross-validated target-prediction task - drug-level
folds, Jaccard-grouped to stop a target's other ligands leaking across the split - comparing the four
representations. That is the test that should gate the FNN investment, and it is now the only one left.

Still true from section 8: feed the FNN the residual plus the global component as one separate feature, not
5,707 features correlated at 0.97.

---

## 12. A549 as the cross-cell CONTROL, and it passes (2026-11-26)

`cell=A549 bash scripts/mc_phenotype_assoc/reproduce_report11.sh 9 10 11`

### 12.0 Read this run for what it is: a control, not a replication

**The Attene-Ramos phenotype is measured in HepG2** (PubChem AID 720637) whatever cell line supplies the
labels. `attene_ramos_mmp.tsv` does not change cell line between these runs - only the MetChange labels do.
So the two runs are:

| run | phenotype | labels | design |
|---|---|---|---|
| section 10 | HepG2 MMP | **HepG2** | **matched** |
| this section | HepG2 MMP | **A549** | **cross-cell control** |

Report 08 section 8.3 states the reading: "Running them against A549 labels is the **cross-cell-line
control** [...] not a replication: if the association is as strong with A549 labels against a HepG2-measured
phenotype, it is a compound property rather than HepG2 metabolism."

**A first version of this section got that backwards**, calling the A549 null a non-replication that cost the
section 10 result its generality. It is the opposite: a null here is the control passing. The error is
recorded because the project had already documented this exact trap, in those words, and it was walked into
anyway.

### 12.1 The arms are not co-extensive in A549, and that had to be fixed first

A549 has **44,396** perturbations in the MetChange arm against **7,636** in the expression arm, sharing
**7,633**. Unrestricted, the comparison would have scored MetChange on nearly six times as many
perturbations - and therefore a different drug set - as the expression arms. HEPG2 does not have this
problem: both arms carry the identical 15,032 perturbations in identical order, so section 10 needed no
restriction and stands as reported.

All four A549 arms are cut to the 7,633 shared perturbations, giving a matched **n = 78 cases / 428 controls**
in every arm. `reproduce_report11.sh` now compares the two headers itself and restricts when they differ.

**A bug this exposed.** `build_derived_labels.py` was symlinking `qvals.csv` and `fit_info_nsig.csv` from the
unrestricted source, so a 7,633-column matrix sat beside 44,397-column q-values. The expression arm died
loudly (`matrix q has a different header to value`); the MetChange arm did **not**, because its symlinked
`fit_info_nsig.csv` satisfied the `n_sig` lookup before `qvals` was read - so it silently took `n_sig` from a
different perturbation set. One cause, a loud failure in one arm and a silent inconsistency in another. Aux
files now carry over only when the perturbation set is identical, and their absence is explained in each
derived directory's `PROVENANCE.txt`.

### 12.2 The control passes: A549 labels carry no information about HepG2 MMP

MetChange arm, `z`, MMP_down vs inactive:

| set | variant | HepG2 labels (n=144/475) | A549 labels (n=78/428) |
|---|---|---|---|
| tca | mito_nohub | 0.362 **rand 0.005** | 0.504 rand 0.93 |
| tca | cyto_nohub | 0.393 **rand 0.005** | 0.483 rand 0.82 |
| tca | all_nohub | 0.383 **rand 0.005** | 0.498 rand 0.96 |
| tca_endpoints | cyto_nohub | 0.357 **rand 0.005** | 0.515 rand 0.91 |
| tca_endpoints | all_nohub | 0.358 **rand 0.005** | 0.512 rand 0.86 |
| pyruvate | mito_nohub | 0.365 **rand 0.015** | 0.491 rand 0.94 |
| pyruvate | all_nohub | 0.389 **rand 0.010** | 0.465 rand 0.58 |

All 18 TCA and pyruvate rows fall between 0.451 and 0.536, random-set p 0.39 to 0.995. The strongest row on
any set is `glycolysis mito` at 0.566, rand 0.095.

This is a null and not an attenuation - the point estimates are centred on chance, and with n = 78/428 the
AUROC standard error is about 0.036, so the HepG2 value of 0.362 sits ~3.6 SE away. **Therefore the section
10 association is not a compound property.** A drug's effect on A549 metabolism says nothing about its effect
on HepG2 mitochondria, which is what should happen if the signal is HepG2 metabolism rather than chemistry.

### 12.3 It reproduces report 09 section 3

Report 09 section 3 ran this control across all eight cell lines with labels and found the MMP-down
association HepG2-only: A549 at 0.477 (own drugs), 0.493 (common set), 0.463 (responsiveness-matched). This
run returns 0.504 through a different pipeline - restricted perturbation set, freshly built derived arms,
random-set nulls - and agrees. A passed reproducibility check on an established control.

### 12.4 What it does and does not establish

**Does:** the section 10 result survives the cross-cell control, so "MetChange beats the producer-average on
HepG2 MMP" is not explained by compound chemistry. Combined with section 10.1, the producer-average null and
the global-residual test, that is three distinct challenges the HepG2 TCA association has now survived.

**Does not:** generality across phenotypes. Every one of these tests uses the same HepG2-measured MMP
phenotype, so they probe different confounds of one association rather than several associations. A second
*phenotype* measured in HepG2 would test that; a second *label* cell line cannot, because the phenotype is
fixed in HepG2.

**And it does not open the 10.3b diagnostics.** The sign conflict lives in the matched HepG2 design, so
A549's absolute expression matrix - the one that happens to be local - is the wrong one. Resolving 10.3b
needs `cell_matrices/HEPG2/expr.csv`, by transfer or by running the diagnostic on the cluster.

### 12.5 Status of the claims in this report

| claim | status |
|---|---|
| MetChange is ~1-dimensional, and the LP is what collapses it (8, 9.2) | **holds**, both cell lines (PR 1.10 / 1.11) |
| MetChange beats the producer-average; the signal is not the global component (10.1) | **holds**, and now also survives the cross-cell control |
| the association is HepG2 metabolism, not compound chemistry (12.2) | **supported** |
| the MMP-down direction has a transcriptional reading (10.3b) | **open, blocked** on HepG2 absolute expression |
| any of this generalises to another phenotype | **untested** - all tests share one HepG2-measured phenotype |
