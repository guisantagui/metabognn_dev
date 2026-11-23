# Report 09 — The MetChange / membrane-potential association is specific to the cell line the phenotype was measured in

*2026-11-23. Cross-cell-line run of the label-level association of reports 07-08, over the 8 LINCS core cell
lines with MetChange labels (plus MCF7, expression only). Inputs: `results/mc_phenotype_assoc/<cell>/`,
all scored against the single union phenotype table `phenotypes_ALL/` so the groups are identical everywhere.
Analysis scripts for this report are the ad-hoc ones listed in section 6; the association itself is
`scripts/mc_phenotype_assoc/` unchanged.*

**Summary.** For **MMP-down**, the MetChange association is **present in HepG2 and absent in the other seven
cell lines** on the identical set of drugs. Since the compounds are the same in all eight comparisons, this
**rules out the compound-class / chemotype explanation that report 08 left open** - a chemical property cannot
act in one cell line and not in seven. HepG2 is where Attene-Ramos measured the phenotype, and it has the
*fewest* drugs of any line, so it is not a power artefact either.

**Two qualifications, added after the first version of this report and corrected below.** (1) **MMP-up is NOT
HepG2-specific**: HCC515 survives responsiveness matching as well (0.624, p = 0.011) and A549 is borderline, so
the specificity claim holds for MMP-down only (section 3.1). (2) The expression arm is **not** "cell-line
general" either - that reading came from pooled numbers; under matching it mostly collapses, in HepG2 to
*exactly* 0.500, surviving only in PC3 and HCC515 (section 4). The sharper and better-supported statement is
within one cell line: **in HepG2, MetChange survives responsiveness matching (0.365, p = 3.4e-04) while
expression goes to chance (0.500, p = 0.99).**

## 1. The prediction, stated before the run

Report 08, section 7 item 2 and the working discussion that preceded this run: the phenotypes are measured in
**HepG2** whatever cell line's labels are used, so every other line is a cross-cell-line **control**, not a
replication. The two outcomes were given in advance:

- **uniform across cell lines** -> the signal is a property of the compounds (chemotype / general bioactivity),
  which is the alternative report 08, section 2 could not exclude;
- **HepG2-specific** -> the signal is a property of HepG2 metabolism.

The result is the second. This is a confirmed prediction rather than a pattern found by searching, which
matters for how much weight it carries.

## 2. Comparability: same phenotype table, same drugs

All eight runs used the union table (`run_info.json`: `phenotype = .../phenotypes_ALL/...`), so group
assignment is identical. But coverage is not: drugs per cell range from **717 (HepG2) to 1,364 (PC3)**, because
each cell line has labels for a different subset. Comparing raw per-cell numbers would therefore confound
coverage with biology.

**Restricted to the drugs present in all 8 cell lines with a single consistent group: 641 drugs, 125 MMP-down
and 427 inactive.** Every number below is on exactly those drugs, so the only thing that varies between
columns is the cell line.

## 3. MetChange, MMP-down: HepG2 only

`tca mito_nohub`, statistic `z`, MMP-down vs inactive:

| cell | own drugs | **common set (125/427)** | responsiveness-matched | responsive share case / control |
|---|---|---|---|---|
| **HEPG2** | 0.362, p = 5e-07 | **0.350, p = 3.3e-07** | **0.365, p = 3.4e-04** (118 pairs) | 28% / 7% |
| A549 | 0.477, p = 0.34 | 0.493, p = 0.80 | 0.463, p = 0.32 | 62% / 37% |
| A375 | 0.507, p = 0.76 | 0.510, p = 0.73 | 0.565, p = 0.082 | 42% / 21% |
| HA1E | 0.470, p = 0.20 | 0.434, p = 0.025 | 0.454, p = 0.21 | 64% / 37% |
| HCC515 | 0.480, p = 0.47 | 0.481, p = 0.51 | 0.472, p = 0.44 | 42% / 23% |
| HT29 | 0.487, p = 0.59 | 0.464, p = 0.22 | 0.473, p = 0.47 | 30% / 16% |
| PC3 | 0.454, p = 0.036 | **0.478, p = 0.46** | 0.453, p = 0.21 | 46% / 25% |
| VCAP | 0.498, p = 0.95 | 0.485, p = 0.61 | 0.480, p = 0.59 | 22% / 10% |

Four things worth drawing out.

1. **HepG2 is the only cell line with a real effect**, and the restriction makes it *stronger* (0.362 -> 0.350),
   i.e. the drugs it shares with the other lines are the ones carrying the signal.
2. **It is not power.** HepG2 has the fewest drugs of any line and the same 125/427 as everyone else here.
3. **PC3's apparent signal was coverage.** Nominally significant on its own drugs (p = 0.036), gone on the
   common set (p = 0.46) - exactly what the restriction exists to catch.
4. **It is not responsiveness.** The responsive share differs a lot between cell lines and HepG2 has the
   largest *ratio* (28% vs 7%), so this had to be checked: matching case and control drugs on
   `log1p(n_sig)` within each cell leaves HepG2 at 0.365 (p = 3.4e-04) and every other line null. HA1E, whose
   pooled value is nominally significant (p = 0.025, uncorrected, and ~0.2 over 8 cell lines), goes to
   p = 0.21 matched.

### 3.1 MMP-up is not HepG2-specific

The first version of this report ran the cross-cell comparison for **MMP-down only** and generalised from it.
That was wrong. `energy_currency mito` (z), MMP-up vs inactive, common set 74/427:

| cell | pooled | matched |
|---|---|---|
| **HEPG2** | 0.665 (5.6e-06) | **0.662 (7.4e-04)** |
| **HCC515** | 0.637 (1.6e-04) | **0.624 (0.011)** |
| A549 | 0.589 (0.014) | 0.589 (0.062) |
| VCAP | 0.553 (0.15) | 0.509 (0.86) |
| HA1E | 0.541 (0.26) | 0.543 (0.37) |
| A375 | 0.524 (0.51) | 0.541 (0.40) |
| HT29 | 0.507 (0.84) | 0.502 (0.97) |
| PC3 | 0.457 (0.24) | 0.412 (0.071) |

**HCC515 survives matching** and A549 is borderline, so the hyperpolarisation association is shared by at
least two or three cell lines. `tca mito_nohub` under MMP-up is null everywhere including HepG2 (0.433,
p = 0.064), confirming report 08's finding that the two phenotype directions live in different metabolite sets.

So **cell-line specificity is a property of the MMP-down / TCA association, not of the dataset as a whole.**
Whether that is meaningful or just the weaker contrast being noisier (74 cases against 125) is not resolved
here.

### 3.2 Why HA1E looks nominally significant

HA1E is the one line with a nominal MMP-down signal (0.434, p = 0.025) and it **does not survive matching**
(0.454, p = 0.21). Two reasons it was always weak: HA1E has the **largest responsive share of any line**
(64% of cases against 37% of controls), so it is the most exposed to the confounder of report 08 section 3;
and p = 0.025 uncorrected over 8 cell lines is about 0.2 after Bonferroni. It is responsiveness plus noise,
not a second positive.

## 4. The expression arm, corrected: not "cell-general", mostly responsiveness

The first version of this report read the expression arm from **pooled** numbers and concluded it was
cell-line general. Under responsiveness matching that does not hold. `tca gene` (z), MMP-down, common set:

| cell | pooled | **matched** |
|---|---|---|
| HEPG2 | 0.418 (0.005) | **0.500 (p = 0.99)** |
| HT29 | 0.377 (2.7e-05) | 0.440 (0.11) |
| **PC3** | 0.400 (6.7e-04) | **0.404 (0.0087)** |
| **HCC515** | 0.430 (0.017) | **0.408 (0.012)** |
| A375 | 0.426 (0.012) | 0.481 (0.61) |
| HA1E | 0.447 (0.073) | 0.505 (0.88) |
| VCAP | 0.446 (0.067) | 0.473 (0.46) |
| A549 | 0.472 (0.34) | 0.514 (0.70) |

`pyruvate gene` (z), MMP-up: HEPG2 0.623 -> matched **0.651 (0.0018)**; A549 0.572 -> 0.553 (0.26);
HA1E 0.581 -> 0.568 (0.15). HepG2 only.

So the expression MMP-down association is **almost entirely responsiveness**: in HepG2 it goes to *exactly*
0.500, i.e. the matched design removes it completely, and it survives only in PC3 and HCC515 - two cell lines
where MetChange shows nothing. That is not "general"; it is a confounder distributed unevenly across cell
lines.

**The claim that survives, and it is a better one.** Within HepG2, on the same drugs and the same matched
design:

| HepG2, MMP-down, matched | AUROC | p |
|---|---|---|
| **MetChange `tca mito_nohub`** | **0.365** | **3.4e-04** |
| expression `tca gene` | 0.500 | 0.99 |

MetChange survives the confounder control that removes the expression association entirely. This needs no
cross-cell argument at all, and it is the cleanest MetChange-vs-expression contrast in the project - report 08
section 13 compared magnitudes *before* matching and found them similar, which is why it concluded there was
no advantage.

## 5. The O-glycan negative control behaves like expression, not like TCA

`negctrl_o_glycan all` (z), same common set: HEPG2 0.590 (p = 0.0021), HA1E 0.582 (p = 0.0055),
VCAP 0.578 (p = 0.0076), HT29 0.564 (p = 0.029), A549 0.552 (p = 0.076), the rest null.

So the set that failed as a negative control in report 08, 3.1 is **not** cell-line specific either. That
supports the reading of report 08, 3.3 - its signal sits in low-degree terminal glycan nodes and behaves like
glycosyltransferase expression - and it shows the cross-cell design **discriminates between the two
findings** rather than flagging everything: TCA is HepG2-only, O-glycan is not.

## 6. What this does and does not establish

**Excluded, and this is the main gain:** the compound-class / chemotype explanation. Identical compounds in
all eight comparisons cannot give a signal in one cell line and nothing in seven. Report 08, section 2 listed
this as the live alternative and as the reason the Hallinger null and the ToxCast-subset finding were
worrying; it is now ruled out for the cross-cell contrast.

**Still not established:**

1. **WITHDRAWN as stated.** The first version said cell specificity "rules out chemistry, not transcription",
   because a HepG2-specific transcriptional programme projecting onto TCA would look identical. That caveat
   was close to vacuous: **MetChange is by construction an LP over the transcriptional signature and the
   stoichiometric matrix**, so "transcription projected onto TCA" is a description of the method, not an
   alternative to the finding.

   **The real limitation is the asymmetry of the design: MMP was measured only in HepG2.** Two readings remain
   and nothing here separates them:
   - **the phenotype itself is cell-specific** - these drugs collapse the potential in HepG2 but perhaps not
     in A549 - and MetChange tracks it correctly in each line, in which case the seven nulls are *correct
     negatives* and the result is strong;
   - **MetChange only works well in HepG2** and the phenotype is in fact shared, in which case the nulls are
     failures.

   We do not know what MMP these compounds would have produced in A549 or HT29. The first reading is
   biologically plausible - mitochondrial toxicity is cell-line dependent, and HepG2 is the standard line for
   it precisely because hepatocytes are mitochondria-rich and drug-metabolising - but that is an argument, not
   evidence.
2. **One phenotype, one measured cell line.** The design has a single positive cell, so "HepG2 metabolism"
   cannot be separated from anything else peculiar to HepG2 (the only liver line here, and the line whose
   labels have been examined most). Confirming the label pipeline was identical across cells is worth doing.
3. **HepG2 was the discovery cell**, so its own significance is not independent evidence. The evidence is the
   *pattern* - one line positive, seven null, on identical drugs.
4. **Effect size is unchanged**, |AUROC - 0.5| ~ 0.15. Specific, not strong.

**The decisive next experiment** is a second (phenotype, cell line) pair: a functional readout measured in a
*different* line, predicting *that* line specifically. A second HepG2 assay (`APR_HepG2_MitoMembPot`) cannot do
it - same line, so it tests reproducibility, not specificity.

**Schuhknecht is not a substitute, although it is A549 metabolomics already on disk.** It was proposed for this
and then withdrawn. The Attene design works because the phenotype is an *independent modality* - a ratiometric
dye - while the predictor is inferred production capacity. A Schuhknecht-derived phenotype would be *measured
metabolite levels* scored over the same metabolites the MetChange set summarises, so the group test would be
the already-null per-drug correlation of report 05, 8.2 (drug-specific rho = -0.013) coarsened into two bins -
the same information at lower resolution. Report 05, 9.3's argument applies at full strength: capacity labels
are not validated by pools.

**What preserves the design: the LINCS genetic perturbations** (STATUS 6.2b item 1), promoted above everything
else for this question because `trt_sh` / `trt_xpr` / `trt_oe` signatures exist **per cell line**. That supplies
(phenotype, cell) pairs across many lines with no new data: the ground truth is exact and cell-independent
(which gene was knocked down) while the expected metabolic consequence is cell-specific, so the symmetric test
is available - does cell X's MetChange read the X knockdown better than cell Y's MetChange reads the same
knockdown? It is also not a metabolite-level readout, so the objection to Schuhknecht does not apply. The known
caveat stands: MetChange reads the perturbed transcript itself, so the 0-reaction version is partly circular
and the 1-2-reactions-away version is the real test.

**Higher-powered but a different phenotype: PRISM / DepMap viability**, thousands of compounds across all nine
core lines. Not metabolic, so it cannot validate the content of the claim, but it tests the *mechanism* of
cell-specificity at large n - does cell X's MetChange predict cell X's sensitivity better than cell Y's does?

**And an honest gap:** for a mitochondrial phenotype in a non-HepG2 LINCS core line we are not aware of one at
screen scale. Everything touched so far - Attene-Ramos, Hallinger, the Apredica `APR_HepG2_MitoMembPot` panel,
the CCTE_Simmons Seahorse assays - is HepG2, which is itself informative: HepG2 is *the* standard line for
mitochondrial toxicity, so a matched second line may simply not exist. Worth a targeted literature search
rather than an assumption.

## 7. Analysis provenance

**`scripts/mc_phenotype_assoc/cross_cell_compare.py`** (added with this report) performs the common-drug
restriction and both designs, for any run / set / statistic / contrast, reading only `per_drug_scores.tsv`:

```
python3 -I scripts/mc_phenotype_assoc/cross_cell_compare.py \
    --run attene_allmodes --set tca --variant mito_nohub --case MMP_down \
    --out results/mc_phenotype_assoc/cross_cell/tca_mito_nohub_MMP_down.tsv
```

It reproduces every number in section 3 exactly. The first version of this report was produced by throwaway
scripts in a session temporary directory, which would have left the central result unreproducible once the
session was cleaned up; that is why this exists.

The underlying runs are reproducible with `use_union=true bash scripts/mc_phenotype_assoc/submit_assoc.sh`
(report 08, 8.4-8.5). The expression and negative-control arms in sections 4-5 use the same script with
`--run attene_expr` / `--run attene_negctrl` and the corresponding set names.
