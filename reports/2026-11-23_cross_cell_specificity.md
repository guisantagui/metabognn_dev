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

**Added 2026-11-25: `scripts/mc_phenotype_assoc/reproduce_report09.sh`** runs every number in sections 3, 3.1,
4 and 5, one block per section (`bash ... reproduce_report09.sh` for all, or `... 4` for one). Verified: the
section-4 head-to-head returns HEPG2 MetChange 0.350 pooled / **0.365 matched** and expression 0.418 pooled /
**0.500 matched**, exactly as reported. For sections 9 and 9.6 use `run_naia_assoc.sh`, which always existed.

The gap it closes: the HepG2 MMP association itself was always scripted - `launch_assoc_HEPG2_slurm.sh` runs
both the MetChange and expression arms with their stratified and responsiveness-matched checks - but the
**cross-cell comparison and the head-to-head were not**, because both apply the common-drug restriction that
the per-cell launchers do not. Section 4 holds the project's most load-bearing number for the central thesis
and was reproducible only from an ad-hoc command line.

One thing the launcher surfaces that the report did not record: the **responsive shares differ enormously
between the two label types**. For MMP-down in HEPG2 the MetChange arm is 28% of cases against 7% of controls,
while the expression arm is 62% against 50%. Expression labels call far more drugs responsive, with a much
smaller case-control gap, which is relevant context for why the matched designs behave so differently - and
worth bearing in mind before reading the 0.500 as a clean null.



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

## 8. Three new candidate datasets, and what each can and cannot do

*Added 2026-11-23 after `data/validation/{naia,nyffeler,xu}/` arrived. Section 6 asked for a second
(phenotype, cell line) pair; these are the first candidates, and one of them can supply it. Overlap figures
below are InChIKey-first-block matches against `lincs_unique_perts.tsv`, which is a **floor** - the proper
mapping needs `map_smiles.R` with salt stripping and stereo handling, which typically recovers more.*

| | naia | nyffeler | xu |
|---|---|---|---|
| cell line | **SH-SY5Y** (differentiated neuroblastoma) | **U2OS** | AC16 cardiomyocytes |
| in LINCS `cellinfo_beta`? | **yes** (`SHSY5Y`) | **yes** (`U2OS`) | **no** |
| phenotype | cellular ATP + cytotoxicity, z-scores | cell painting, 10 `Mito_*` feature categories | **mitochondrial membrane potential** |
| design | single concentration | **concentration-response (BMD)** | AC50 + active/inactive call |
| compounds | 1,199 | 462 (APCRA, 24 h, test chemicals) | 247 |
| **matching LINCS** | **683** (+29 by name) | 87 with a mito-selectivity value | 53 |
| source | Naia et al., *BMC Biol* 2021, PMC7989211 | EPA Nyffeler, U2OS HTPP | Xu et al., `10.1021/EHP.6c00052` table S10 |

### 8.1 naia - DISQUALIFIED for the symmetric test: 10 compounds overlap

**Resolved 2026-11-23 by the gating check of 8.4, with `data/LINCS/lincs_unique_perts_SHSY5Y.tsv`.** LINCS has
only **224 perturbations in SHSY5Y** against 5,494 in HEPG2, and the Prestwick FDA-approved library barely
touches that panel:

| naia compounds | n |
|---|---|
| screened | 1,199 |
| resolved to a structure | 1,193 |
| **matching LINCS SHSY5Y** | **10** |
| matching LINCS HEPG2 (reference) | 340 |

Of those 10: **1 ATP-down, 0 ATP-up, 4 ATP-flat.** The design of the rest of this section is impossible - and
the single ATP-down compound is **D-cycloserine** (z ATP -2.11), a peptidoglycan-synthesis inhibitor, i.e. one
of the implausible hits flagged below as probable single-concentration noise.

**The gating check earned its keep.** The 683 whole-LINCS structure matches said nothing about whether those
compounds were ever profiled in the cell line that matters, and the differentiation worry below turned out to
be moot.

**What naia is still good for: it replaces xu as the falsification test, with 6x the power.** 340 of its
compounds are in the HEPG2 panel, so the converse question can be asked - *do HepG2 labels predict an ATP
phenotype measured in a different cell line?* That is a better version of 8.2: a metabolic readout rather than
membrane potential, an independent cytotoxicity covariate, and 340 compounds instead of 53. **Read the
direction carefully**, because it is not the same design as section 3:

- section 3 fixes the phenotype (HepG2 MMP) and varies the **label** cell line - only HepG2 labels worked;
- this test fixes the labels (HepG2) and varies the **phenotype** cell line.

If HepG2 labels do predict SH-SY5Y ATP, the two results together are an odd asymmetry - other lines' labels
cannot predict HepG2's phenotype, yet HepG2's labels predict other lines' phenotypes - which would point at
**reading 2 of section 6**: HepG2 labels are simply the better-quality set, and MetChange works where the LINCS
data is deepest. That is the alternative this report could not exclude, so the test is worth running for the
same reason xu is: it can damage the conclusion.

### 8.1a The original assessment, retained for the record

SH-SY5Y is in LINCS and **683 of 1,199 compounds match**, so this is the one dataset that could give a second
(phenotype, cell line) pair. Two further advantages:

- **The phenotype is ATP itself** - a metabolite in the `energy_currency` set, now producer-only in the
  network (report 08, 9). The correspondence between readout and predictor is more direct than for membrane
  potential.
- **It ships an independent cytotoxicity z-score per compound**, which is a *better* confounder control than
  the `n_sig` responsiveness proxy used throughout reports 07-09. And it is genuinely independent:
  **Spearman(z ATP, z cytotox) = -0.083** (p = 0.004) over all 1,199 compounds, so ATP loss here is *not*
  simply cell death. That is unusual and valuable.

Weaknesses, which are real:

- **Single concentration**, so no potency and more noise than a BMD design. The noise is visible: three
  peptidoglycan-synthesis inhibitors sit in the ATP-down non-cytotoxic set, and beta-lactams have no mammalian
  target, so some of that group is measurement noise.
- **The case groups shrink fast.** ATP-down (z <= -2) is 116 compounds, 55 of them LINCS-matched; requiring
  non-cytotoxic (|z cytotox| < 1) leaves 48, of which roughly 27 are LINCS-matched. Against a flat-ATP control
  group of about 280 (160 matched) that is much less power than Attene's 125/427 - an effect the size of the
  HepG2 one (AUROC 0.35) would land near p = 0.02 rather than 3e-07. **A null would therefore be weak evidence
  of absence.**
- It was designed as an **enhancer** screen, so the authors' own finding is ATP-*up* (60 compounds, 42
  non-cytotoxic, about 40 LINCS-matched). That direction is as usable as the down direction and should be
  tested as its own contrast, not pooled.

**And a mismatch that may matter more than the power problem.** Naia screened **differentiated** SH-SY5Y
(the abstract says "differentiated human neuroblastoma cells"), while LINCS `SHSY5Y` signatures are
undifferentiated. This is not a cosmetic difference for a mitochondrial ATP phenotype: differentiation shifts
SH-SY5Y from predominantly glycolytic to oxidative metabolism, which is precisely the axis the phenotype
measures. So "cell-matched" here means matched on cell *line*, not on metabolic state - and a null could then
mean the labels describe a glycolytic cell while the phenotype was measured in an oxidative one. Worth
confirming the differentiation protocol in the paper's methods before treating this as the symmetric test;
it is the one caveat that could disqualify naia even with good LINCS depth.

### 8.2 xu - cannot do the cell-matched test, but is the sharpest test *against* this report

AC16 is **not in LINCS**, so no AC16 labels can be computed and the symmetric design is impossible. Its value
is the opposite: it is **membrane potential, the same phenotype as Attene-Ramos**, measured by a different
group in a different cell type (cardiomyocytes, mitochondria-rich).

That makes it a **falsification test of section 3**. If HepG2 MetChange labels predict *AC16* membrane
potential about as well as they predict HepG2 membrane potential, then the signal travels with the compounds
after all, the chemotype explanation revives, and the cell-specificity claim of this report is in serious
trouble. Worth running **because it can only damage the current conclusion** - and note the asymmetry when
reading the outcome: with 53 matched compounds (61 active / 123 inactive / 63 inconclusive before matching) a
**positive result is informative and a null is not**.

### 8.3 nyffeler - best design, not enough overlap

The strongest experiment of the three and the weakest for us. It is **concentration-response**: a BMD per
feature category, so potency rather than a hit call. And it has an internally controlled mitochondrial
contrast - **mito-selectivity**, `min(Mito_* BMD) - min(other BMD)`, which asks whether the mitochondrial
features respond at a *lower* concentration than everything else and therefore nets out general bioactivity,
the confounder that has dominated reports 07-09.

But: only **87 of 457** compounds with a mito-selectivity value match LINCS, and of those only **19 are
mito-selective** against 68 not, so a binary design is dead. There is also effectively **no inactive group** -
457 of 462 test chemicals have at least one responding `Mito_*` category, because this is a curated bioactive
library. The usable form is a *continuous* rank correlation between mito-selectivity and a MetChange set score
over those 87 drugs, which is worth an afternoon but will not settle anything on its own.

### 8.4 The gating check, before any of this is worth doing

`cellinfo_beta.txt` lists **every** cell line LINCS ever touched, including lines with a handful of
signatures. It says nothing about depth, and we have no `siginfo` locally. **If SHSY5Y has only a few dozen
compounds profiled, naia is dead regardless of its 683 structure matches**, and the same for U2OS and
nyffeler. On the HPC, where the full LINCS data lives:

```bash
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}
            $h["pert_type"]=="trt_cp" && ($h["cell_iname"]=="SHSY5Y" || $h["cell_iname"]=="U2OS") \
            {n[$h["cell_iname"]]++; seen[$h["cell_iname"]" "$h["pert_id"]]=1}
            END{for(c in n) print c, n[c] " signatures";
                for(k in seen){split(k,a," "); d[a[1]]++} for(c in d) print c, d[c] " distinct pert_ids"}' \
    /scratch/users/gsantamaria/test_large_files/metabognn_dev/data/LINCS/siginfo_beta.txt
```

The number that matters is **distinct `pert_id`s overlapping the validation set**, not total signatures. For
naia that needs to be in the hundreds for the design in 8.1 to have the power described; for U2OS, 87 is
already the ceiling.

## 9. The falsification test, run: POSITIVE, but through a different pathway

*Naia SH-SY5Y ATP phenotype scored against **HEPG2** MetChange labels, the converse design of section 8.1.
Built with `scripts/mc_phenotype_assoc/build_naia_phenotype.py` (new), run with the unchanged
`mc_phenotype_assoc.py`, controls with `scripts/mc_phenotype_assoc/covariate_matched_check.py` (new).
Phenotype `results/mc_phenotype_assoc/phenotypes_ALL/naia_atp.tsv`, results
`results/mc_phenotype_assoc/HEPG2/naia_atp/`.*

Groups, from the ATP z-score with a symmetric dead band (case |z| >= 2, control |z| < 1), mapped to the HEPG2
panel: **25 ATP_down, 27 ATP_up, 226 inactive** (340 of 1,199 naia compounds matched by InChIKey, 18 more by
name, 90 excluded as intermediate).

### 9.1 Result

`ATP_down vs inactive`, statistic `z`, with the size-matched random-set control that reports 07-09 treat as the
real test:

| set | variant | AUROC | MW p | perm p | **random-set p** | BH q |
|---|---|---|---|---|---|---|
| **glycolysis** | cyto | **0.265** | 1.1e-04 | 2e-04 | **0.005** | **0.005** |
| glycolysis_endpoints | all | 0.280 | 3.1e-04 | 6e-04 | 0.005 | 0.006 |
| ppp | all | 0.303 | 0.0012 | 0.001 | 0.030 | 0.012 |
| pyruvate | cyto | 0.306 | 0.0015 | 0.0014 | 0.025 | 0.012 |
| tca | cyto | 0.337 | 0.0075 | 0.008 | **0.080** | 0.025 |
| tca | mito_nohub | 0.352 | 0.015 | 0.014 | **0.194** | 0.032 |
| oxphos | mito | 0.428 | 0.24 | 0.24 | 0.61 | 0.28 |
| **energy_currency** | all | **0.482** | 0.77 | 0.77 | 0.89 | 0.81 |
| redox_couples | all | **0.648** | 0.015 | 0.016 | 0.100 | 0.032 |

**`ATP_up vs inactive` is null everywhere** - every BH q = 0.94, AUROC 0.44 to 0.56. So only ATP *loss* is
predicted, not the enhancement the paper was looking for.

Three features of that table matter more than the headline:

1. **It is glycolysis and the cytosolic sets, not TCA or OXPHOS.** TCA is nominally significant by
   Mann-Whitney but **fails the random-set control** (p = 0.08 to 0.19), which is the control that has
   governed every claim in reports 07-09. `oxphos mito` is flatly null.
2. **`energy_currency` is null (0.482, p = 0.77)** even though the phenotype *is* ATP and ATP is a member of
   that set. Whatever is being detected, it is not a direct read-out of the ATP node.
3. **`redox_couples` moves in the opposite direction** (0.648). A single "this compound is bioactive" factor
   would push every set the same way; these do not move together.

### 9.2 It is not cytotoxicity, and it is not a few compounds

naia ships a measured per-compound cytotoxicity z-score, so viability can be controlled **directly** rather
than through the `n_sig` proxy used in reports 07-09. The imbalance was severe - the ATP-down group averaged
+6.62 z-units more cytotoxic than the controls - and matching removes it:

| `glycolysis cyto`, ATP_down | AUROC | p | n |
|---|---|---|---|
| pooled | 0.265 | 1.1e-04 | 25/226 |
| **matched on measured cytotoxicity** | **0.173** | **0.0025** | 15 pairs |
| matched on `n_sig` responsiveness | 0.296 | 0.014 | 25 pairs |
| matched on both | 0.222 | 0.010 | 15 pairs |

**The association survives cytotoxicity matching.** The first write-up said it "gets stronger, so cell death
was suppressing it" - that over-read the numbers and is **withdrawn**: the matched estimate is a different
comparison (15 cases against 15 *matched* controls, not 25 against 226), so its being lower than the pooled
value is not evidence of anything. The claim that holds is **survival**, not strengthening.

What matching does and what it costs. The groups are badly imbalanced on cytotoxicity - cases median
**+1.91** (IQR -0.31 to 7.10, max 31.0), controls median **-0.71** (IQR -1.16 to -0.24, max 2.7) - so the
pooled test cannot distinguish metabolism from cell death. Matching pairs each case with an unused control of
similar cytotoxicity and compares only within pairs; the balance line `z_cytotox +6.62->-0.01` is the
diagnostic that it worked. But the controls only span up to 2.7, so **10 of the 25 cases are too cytotoxic to
have any counterpart and are discarded**. The matched test therefore answers a narrower question: *among
compounds whose cytotoxicity falls inside the range the controls span*, the association holds (0.173,
p = 0.0025). It says nothing about the 10 most cytotoxic cases.

**A stratified check, which keeps all 25 cases, is the better support** and was run after the fact:

| cases, split by cytotoxicity | n | median glycolysis z | vs the 226 controls |
|---|---|---|---|
| z_cytotox < 1 | 10 | -0.17 | AUROC 0.236, p = 0.0048 |
| z_cytotox >= 1 | 15 | -0.10 | AUROC 0.284, p = 0.0052 |
| (controls) | 226 | +0.06 | |

**Both strata separate from the controls, to the same degree**, and within the cases glycolysis z barely
tracks cytotoxicity at all (Spearman +0.156, p = 0.46). So the association is not carried by the dying
compounds, and it is not carried only by the surviving ones either - it is present across the cytotoxicity
range, which is a stronger statement than the 15-pair matched test alone. `energy_currency` stays null under every design (0.40 to 0.48). Robustness:
leave-one-out never moves p above 3.1e-04; dropping the five most extreme cases still gives AUROC 0.321
(p = 0.008, n = 20); and the MoA composition is heterogeneous - 9 of 25 unannotated, no class contributing
more than 2.

The case group is also chemically plausible in a way the full screen was not: it contains **hexachlorophene**
and **triclosan** (both classical uncouplers), **proscillaridin, digoxin and digoxigenin** (cardiac
glycosides, which raise ATP demand through Na+/K+-ATPase), **chloroxine** and **parthenolide**. The
peptidoglycan-synthesis inhibitors that made the unmapped portion look noisy (8.1) are not in the mapped set.

### 9.3 What this does to section 3

**The strong reading of this report is refuted.** HepG2 labels are *not* specific to HepG2-measured
phenotypes: they predict an ATP phenotype measured in a different cell line, at q = 0.005, surviving a direct
viability control. And the asymmetry predicted in 8.1 is now observed - other lines' labels cannot predict
HepG2 MMP (section 3), yet HepG2 labels do predict SH-SY5Y ATP.

**But the two associations are not the same signal**, which is what stops this from collapsing into the
chemotype explanation:

| | section 3 (MMP, HepG2) | section 9 (ATP, SH-SY5Y) |
|---|---|---|
| strongest set | `tca mito_nohub` 0.350 | `glycolysis cyto` 0.265 |
| TCA | **yes**, q = 4e-06 | **no** (random-set p = 0.08-0.19) |
| glycolysis | not the hit | **yes**, q = 0.005 |
| `energy_currency` | 0.525, null | 0.482, null |
| compartment | mitochondrial | cytosolic |

Different phenotype, different pathway, different compartment. A generic bioactivity or chemotype axis would
have produced the *same* sets in both.

**And the depth explanation - reading 2 of section 6 - is now dead.** `data/LINCS/lincs_unique_perts_*.tsv`:

| cell | compounds profiled | section 3 result |
|---|---|---|
| **A375** | **13,314** | 0.510, p = 0.73 - the most null of all eight |
| HEPG2 | 5,494 | 0.350, p = 3.3e-07 |
| SHSY5Y | 224 | (too thin to test, 8.1) |

A375 has **2.4x more profiled compounds than HepG2 and shows nothing**. HepG2 also had the *lowest* responsive
share of the eight lines (28% of cases vs 7% of controls, section 3), so its labels are not "more
significant" either. "MetChange works where the LINCS data is deepest" does not survive contact with those
numbers.

**Revised statement.** Section 3 establishes that among label cell lines, only HepG2's predict a
HepG2-measured membrane-potential phenotype, and that this is not explained by chemistry, power, coverage,
responsiveness or panel depth. It does **not** establish that HepG2 labels are inert elsewhere - section 9
shows they are not. The coherent reading of both is that MetChange tracks compound-induced metabolic shifts
with the *pathway matching the phenotype* - mitochondrial sets for membrane potential, glycolytic sets for
ATP content - rather than that one cell line's labels are magic. **"The association is HepG2-specific" is no
longer sayable without naming which association.**

### 9.4 Caveats, including one against myself

- **n = 25 cases and 15 matched pairs.** Small, and the matched estimate in particular rests on 15 pairs.
- **The pathway assignment is post hoc.** 160 set x variant x statistic combinations were scored per contrast.
  BH q is reported and glycolysis survives at q = 0.005, but *highlighting glycolysis* was data-driven; the
  comparison table in 9.3 would be far stronger if the glycolysis prediction had been registered first.
- **The glycolysis reading is biologically odd and unexplained.** Naia used *differentiated* SH-SY5Y, which
  shifts toward oxidative metabolism, so mitochondrial sets were the a priori expectation for an ATP
  phenotype. Getting glycolysis instead may mean the labels (HepG2, a glycolytic hepatoma line) are reporting
  HepG2's own metabolic wiring rather than the phenotype cell's - which would be a label-side artefact that
  happens to correlate with ATP loss. This is the most likely remaining confounder and it is not controlled.
- **Single concentration**, so no potency axis and no dose-response confirmation.
- The phenotype and label cell lines differ, so nothing here is a cell-matched test; section 8.1 established
  that no cell-matched test is possible with these data.

### 9.5 The two phenotypes are largely different compounds

Worth recording, because it bears on whether sections 3 and 9 are one signal or two. The Attene-Ramos MMP
table and the naia ATP table share **238 drugs**. Cross-tabulating their group calls:

| naia \ Attene | MMP_down | MMP_up | cytotox | inactive | incon |
|---|---|---|---|---|---|
| **ATP_down** | 9 | 3 | 2 | 9 | 0 |
| **ATP_up** | 7 | 0 | 0 | 14 | 2 |
| **inactive** | 32 | 13 | 2 | 120 | 25 |

The two readouts are **correlated but far from the same thing**: on the case/control cells alone,
[[9, 9], [32, 120]], an ATP-down drug is **3.75x** more likely to be MMP-down than an inactive drug is
(Fisher p = 0.016) - yet only 9 of the 23 ATP-down drugs with an MMP call are MMP-down, and 9 are
MMP-*inactive*. Of the 48 shared MMP-down drugs only 9 are ATP-down.

So section 3 and section 9 are not the same experiment re-run. They use different compound libraries
(Tox21 vs Prestwick FDA-approved), different assays (ratiometric dye vs ATP luminescence), different cell
lines for the phenotype, and substantially non-overlapping active sets. That they implicate different
metabolite sets is therefore unsurprising - and it is the reason a single shared factor, chemical or
stress-related, is a poor explanation for both at once.

### 9.6 Head to head with the expression baseline: MetChange survives the viability control, expression does not

*Run after 9.5, on the same phenotype table and the same 25/226 drugs.
`results/mc_phenotype_assoc/HEPG2/naia_atp_expr/`, HEPG2 expression labels with the GPR gene sets.*

**Expression also predicts SH-SY5Y ATP loss - and through a different pathway than MetChange.** Its best sets
are mitochondrial: `oxphos gene` AUROC 0.314 to 0.348 depending on statistic (random-set p = 0.015 for
`sign_mean`), then `ppp` 0.355, `tca` 0.347, while `glycolysis gene` is null (0.445). MetChange's ordering was
the reverse: glycolysis first, `oxphos mito` null. So the two label types are **not** reading the same thing,
and expression's ordering is the biologically expected one for uncoupler-driven ATP loss.

Then the control that matters. Matching on naia's **measured per-compound cytotoxicity**, `ATP_down`:

| labels | set | pooled | **cytotox-matched** | matched on both |
|---|---|---|---|---|
| **MetChange** | glycolysis cyto | 0.265 (p = 1.1e-04) | **0.173 (p = 0.0025)** | **0.222 (p = 0.010)** |
| **MetChange** | ppp all | 0.303 (p = 0.0012) | **0.156 (p = 0.0014)** | **0.231 (p = 0.013)** |
| **MetChange** | pyruvate cyto | 0.306 (p = 0.0015) | **0.187 (p = 0.0037)** | **0.240 (p = 0.016)** |
| **MetChange** | tca cyto | 0.337 (p = 0.0075) | 0.280 (p = 0.042) | 0.311 (p = 0.082) |
| MetChange | energy_currency all | 0.482 (p = 0.77) | 0.400 (p = 0.36) | 0.444 (p = 0.62) |
| expression | **oxphos gene** | 0.348 (p = 0.013) | **0.329 (p = 0.115)** | 0.329 (p = 0.115) |
| expression | ppp gene | 0.355 (p = 0.017) | 0.302 (p = 0.068) | 0.316 (p = 0.089) |
| expression | tca gene | 0.397 (p = 0.091) | 0.311 (p = 0.082) | 0.378 (p = 0.26) |
| expression | glycolysis gene | 0.445 (p = 0.37) | 0.422 (p = 0.48) | 0.409 (p = 0.41) |

**Expression's association with ATP loss is viability-driven: it does not survive matching on measured
cytotoxicity (oxphos p = 0.013 -> 0.115). MetChange's does survive (p = 1.1e-04 -> 0.0025).** See 9.2 for why
"survives" and not "strengthens", and for the stratified check that keeps all 25 cases.

This is the **same pattern as section 4, independently replicated**, and the independence is the point:

| | section 4 | section 9.6 |
|---|---|---|
| phenotype | membrane potential | **cellular ATP** |
| phenotype cell line | HepG2 | **SH-SY5Y** |
| compound library | Tox21 | **Prestwick FDA-approved** |
| confounder control | `n_sig` responsiveness **proxy** | **measured cytotoxicity** |
| MetChange, matched | 0.365 (p = 3.4e-04) | 0.173 (p = 0.0025) |
| expression, matched | **0.500 (p = 0.99)** | **0.329 (p = 0.115)** |

Four things differ between the two rows - phenotype, phenotype cell line, compound library, and the nature of
the confounder control - and the same asymmetry appears in both. That is a much better basis for the claim
that **the stoichiometric layer carries metabolic information the transcriptome does not** than report 08
section 13, which compared unmatched magnitudes in a single cell line and concluded there was no advantage.

> **ADDED 2026-11-24, and it qualifies the above. BOTH ROWS USE HEPG2 LABELS.** There are no SH-SY5Y labels and
> none are possible: LINCS has 224 perturbations in SHSY5Y and only 10 Naia compounds are among them (8.1).
> `results/labels/` contains no SHSY5Y directory. Both `naia_atp` and `naia_atp_expr` have
> `labels = results/labels/HEPG2/{metchange/rank,expr}`; SH-SY5Y is the cell the **phenotype** was measured in,
> nothing more.
>
> So the two instances are independent in four respects and **dependent in the one that may matter most**. If
> anything about HepG2's label generation systematically favours MetChange over expression - its reference
> level, its limma fit, or the fact that it has the lowest responsive share of the eight cell lines (section 3)
> - both instances inherit it. Describing them as simply "independent" overstated the case.
>
> The missing axis is a MetChange-versus-expression comparison on **non-HepG2 labels**. The PRISM runs cover
> five such cell lines (A375, A549, HT29, MCF7, PC3) but were executed with the MetChange arm only - there is
> no `prism_viability_expr` - so they do not yet speak to this claim at all. Adding that arm is the cheapest
> available test of the project's central thesis (report 10, 10.5).

**What it does not license.** The effect sizes are small (|AUROC - 0.5| = 0.15 to 0.33), the matched estimates
rest on **15 pairs**, and the pathway was chosen post hoc in both arms. And two failures remain unexplained
and must not be dropped from the account: the Schuhknecht per-drug metabolite-level correlation was null
(report 05, 8.2, rho = -0.013) and the Hallinger OCR association was null (report 07). A claim that MetChange
tracks real metabolism has to accommodate those, not just the positives - the most likely reconciliation
being that set-level group contrasts are a far weaker demand than per-drug profile agreement (section 8.1's
reasoning about why Schuhknecht could not substitute), but that is an argument, not a result.

### 9.7 How to run section 9

One command from the repository root reproduces every number in 9.1 to 9.6:

```bash
bash scripts/mc_phenotype_assoc/run_naia_assoc.sh          # default label cell line HEPG2
bash scripts/mc_phenotype_assoc/run_naia_assoc.sh A375     # or any cell with labels + a panel file
```

It chains five steps, each runnable on its own:

| step | script | what it produces |
|---|---|---|
| 1 | `build_naia_phenotype.py` | `phenotypes_ALL/naia_atp.tsv` - groups from the ATP z-score, drugs mapped to that cell's LINCS panel |
| 2 | `compare_phenotype_tables.py` (new) | the 9.5 cross-tab and its Fisher test |
| 3 | `mc_phenotype_assoc.py` | 9.1, MetChange labels with the metabolite sets |
| 4 | `mc_phenotype_assoc.py` | 9.6, expression labels with the GPR gene sets |
| 5 | `covariate_matched_check.py` (new) | the matched controls of 9.2 and 9.6; `--robustness` adds the leave-one-out, drop-5, MoA composition and case listing |

`covariate_matched_check.py` generalises `matched_check.py` to any numeric phenotype column, which is what
allows matching on naia's **measured** cytotoxicity rather than the `n_sig` responsiveness proxy - the
difference that makes 9.6 a stronger test than report 08 section 13.

**Note on provenance.** The first write-up of 9.2 and 9.5 was produced by inline one-off commands rather than
scripts, repeating the mistake recorded in section 7. `--robustness` and `compare_phenotype_tables.py` exist
so that nothing in section 9 depends on a command that is not in the repository; the numbers above were
re-derived from these scripts, and the Fisher p-value in 9.5 is new because the one-off had only estimated the
ratio.
