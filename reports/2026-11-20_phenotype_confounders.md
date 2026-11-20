# Report 08 — Is the MetChange / membrane-potential association chemistry or metabolism?

*2026-11-20. No new model runs and no new label runs. Confounder analysis of the association found in
report 07 (`2026-11-19_validation_datasets.md`, sections 11–17), plus the compound-subset finding that
explains why Hallinger is null. Code: `scripts/mc_phenotype_assoc/{confounder_checks,matched_check}.py`.
Statistics and conventions: report 07, section 14.2 (three levels of aggregation; the test is per set).*

**Summary.** The association survives responsiveness matching and coarse mechanism-family stratification. It
does **not** show a potency or dose-exposure gradient, it is **absent in the ToxCast portion of the compound
library** — which is why Hallinger is null — and it cannot be shown to persist within narrow mechanism classes,
where power is gone. **Decisively for the interpretation (3.1): of five pathways declared irrelevant in
advance, four are null but one (O-glycan metabolism) shows an equal and opposite association that passes every
specificity test the TCA result passes** — including beating its own size-matched random null. So what the data
support is a multi-pathway difference between phenotype-active and inactive compounds, **not** a TCA-specific
metabolic response. Net: a real group difference in the labels, with its pathway interpretation unsupported and
a compound-class explanation still live.

> **Correction carried from the working session.** On the first pass I ran the within-MoA test with strata
> that double-counted drugs with several MoA labels, no coarse grouping, an unpaired test on matched pairs and
> no panel-level correction. On that version every result attenuated to non-significance and I reported that
> the pattern "fits a compound-class effect better than a graded metabolic response". **That reading was an
> artefact of the broken test and is withdrawn** (section 3).

## 1. Why this had to be done

Two facts moved the chemotype question from optional to decisive:

1. **Matching on `n_sig` equalises how *much* a drug responds, not *what kind* of compound it is**
   (report 07, 16.1). A chemotype or mechanism class shared by the phenotype-active drugs could project onto
   TCA metabolites through the GEM with no real metabolic effect.
2. **The signal is carried by one part of the library** (section 2).

## 2. The compound-subset finding, and why Hallinger is null

All 115 Hallinger-mapped compounds are inside the 849 Attene-Ramos-mapped set, so the Attene contrast can be
restricted to them. On `tca mito_nohub` (statistic `z`), MMP-down vs inactive:

| population | n | AUROC | p | power for 0.364 |
|---|---|---|---|---|
| all mapped compounds | 147/480 | **0.364** | 5.9e-07 | 99% |
| only compounds also in Hallinger (ToxCast) | 38/55 | **0.447** | 0.39 | 60% |
| only compounds **not** in Hallinger | 109/425 | **0.352** | 1.8e-06 | 99% |

**The association is carried entirely by the non-ToxCast compounds.** In the ToxCast subset it is 0.447 here
and 0.483 in Hallinger itself — two independent estimates agreeing that the effect is not there for those
chemicals.

**So Hallinger's null is not primarily an n problem.** Power at 29/86 to detect AUROC 0.371 is 54%, so a null
p-value alone would be uninformative — but the point estimates are the information, and for the sets that
carry the Attene signal they sit at 0.47–0.55, not 0.37:

| set / statistic | Hallinger active [CI] | Attene |
|---|---|---|
| `tca mito_nohub` z | 0.483 [0.356, 0.597] | 0.364 |
| `tca_endpoints cyto_nohub` z | 0.468 [0.338, 0.582] | 0.359 |
| `pyruvate mito_nohub` z | 0.484 [0.357, 0.603] | 0.363 |
| `tca cyto` sign_dev | 0.550 (opposite sign) | 0.415 |

Hallinger **is** the ToxCast Phase I/II library, so this is one consistent observation, not two independent
failures. Two readings remain open: the effect is compound-class-specific for a biological reason (drugs with
defined protein targets give coherent transcriptional responses that the GEM can integrate; pesticides and
industrial chemicals act promiscuously or non-transcriptionally), or it is a chemotype artefact.

## 3. Mechanism-class stratification

`confounder_checks.py`, van Elteren stratified Wilcoxon (mid-ranks, tie-corrected variance, stratum weight
1/(n+1)), with a cluster bootstrap over strata for the interval. MoA and SMILES from
`compoundinfo_beta.txt`, joined on `cmap_name`. **One MoA per drug** (`--moa-pick common`: the most frequent
of its labels), because 92 of the 480 annotated drugs carry several and the first version counted them in
every matching stratum. Two granularities:

- **fine** — the free-text labels themselves, 161 classes in this population, 27 with both groups present;
- **coarse** — 9 mechanism families built from the vocabulary actually present (GPCR/neurotransmitter 166
  drugs, nuclear receptor ligand 81, other 77, metabolic enzyme 68, ion channel 28, DNA/cytoskeleton/cytotoxic
  27, kinase inhibitor 21, transporter/pump 9, protease 3).

p-values are Holm- and BH-corrected across the 6-spec panel.

| spec | reference | fine MoA (27 strata) | **coarse MoA (9 strata)** | Holm |
|---|---|---|---|---|
| `tca mito_nohub` z, MMP-down | 0.364 | 0.419, p = 0.13 | **0.380** [0.302, 0.442], p = 0.0014 | **0.0054** |
| `tca_endpoints cyto_nohub` z | 0.359 | 0.440, p = 0.26 | **0.388** [0.300, 0.457], p = 0.0028 | **0.0083** |
| `tca cyto` sign_dev | 0.415 | 0.447, p = 0.32 | 0.439 [0.351, 0.495], p = 0.10 | 0.20 |
| `energy_currency mito` z, MMP-up | 0.660 | 0.608, p = 0.073 | **0.657** [0.607, 0.752], p = 1.7e-04 | **0.0010** |
| `adenylates all` z, MMP-up | 0.643 | 0.534, p = 0.57 | **0.653** [0.602, 0.776], p = 2.6e-04 | **0.0013** |
| `oxphos mito` z (power diagnostic, **not** a negative control - see 3.1) | 0.476 | 0.504, p = 0.94 | 0.469 [0.401, 0.527], p = 0.41 | 0.41 |

**Within coarse mechanism families the effects survive essentially unattenuated** (0.364 -> 0.380;
0.660 -> 0.657) and survive Holm correction. `oxphos mito` is null in every stratification (but see 3.1 for
what that does and does not establish).

**Which to believe.** The coarse test controls mechanism *family* with real power; the fine test controls
mechanism better but has almost none — 27 strata of a few drugs each, and the population is dominated by
receptor-targeting drugs (of the fine strata, essentially one has a metabolic MoA). So the fine test is both
underpowered and a different biological population, and its attenuation is ambiguous between "confounder
removed" and "population changed". Note also that for `tca mito_nohub` the fine stratification's bootstrap CI
still excludes 0.5 ([0.332, 0.494]) while its van Elteren p is 0.13; **the van Elteren p is the one to quote**,
since the cluster bootstrap of a weighted AUROC understates uncertainty with small strata.

Conclusion: **mechanism family does not explain the association, and narrow mechanism class cannot be tested at
this n.**

### 3.1 Declared negative controls: four behave, one does not

**A correction first.** Through sections 3-5 above and in report 07 I repeatedly called `oxphos mito` a
"negative control". **It is not one, and the label was wrong in a way that mattered.** It is the 13
mitochondrial metabolites of the Oxidative phosphorylation subsystem (ADP, ATP, FAD, FADH2, H+, H2O, NAD+,
NADH, O2, Pi, PPi - all hubs - plus ubiquinol and ubiquinone). Three objections:

1. it was never pre-specified; I designated it post hoc *because* it came out null, which is circular;
2. biologically it is the **most on-target** set for a membrane-potential phenotype, so it is a failed
   *positive* control being used as a negative one - and its being null is a point against the finding, not a
   reassurance;
3. 11 of its 13 members are hubs, so it may be insensitive for structural reasons (report 07, 11.2).

What it legitimately shows is that the matched and stratified designs **have not simply lost power**: if
matching pushed every set to 0.5 one could not distinguish a removed confounder from a dead test. It is a
**power diagnostic** and is now labelled as such throughout.

**Proper negative controls, declared before looking.** `build_metabolite_sets.py --emit negative-controls`
writes five Human-GEM subsystems selected on biology alone, by criteria fixed in advance: not mitochondrial;
no role in respiration; no role in supplying redox cofactors to the ETC (so folate, nicotinate/NAD,
riboflavin/FAD, pantothenate/CoA, ubiquinone synthesis, lipoic acid, porphyrin/heme and ROS detoxification are
excluded); not central carbon catabolism; not plausibly linked to the drug classes in this library (so the
steroid pathways are excluded, the library being rich in glucocorticoid and estrogen receptor ligands); and
large enough to test. What survives those filters is glycan and glycosaminoglycan metabolism. They are written
to their own file so they do not enlarge the BH family of the hypotheses of interest.

MMP-down vs inactive, variant `all`, same pipeline (`attene_negctrl/`):

| set | statistic | AUROC [95% CI] | p | rand |
|---|---|---|---|---|
| `negctrl_keratan_degradation` (77) | z | 0.476 [0.421, 0.522] | 0.38 | 0.64 |
| `negctrl_keratan_biosynthesis` (70) | z | 0.481 [0.422, 0.531] | 0.49 | 0.75 |
| `negctrl_blood_group` (59) | z | 0.489 [0.435, 0.539] | 0.69 | 0.85 |
| `negctrl_cs_hs_biosynthesis` (58) | z | 0.527 [0.469, 0.576] | 0.33 | 0.57 |
| **`negctrl_o_glycan` (30)** | **z** | **0.596 [0.543, 0.647]** | **4.3e-04** | **0.055** |
| `negctrl_o_glycan all_nohub` (26) | z | 0.592 [0.540, 0.646] | 6.9e-04 | 0.040 |
| `negctrl_o_glycan all_nohub` | auroc | **0.631 [0.575, 0.682]** | **1.4e-06** | **0.0099** |
| *for comparison* `tca mito_nohub` | z | 0.364 [0.316, 0.409] | 5.9e-07 | 0.0099 |

**Four of the five are null, which is what a negative control should do** - and it matters that they are,
because it rules out the worst possibility, that any set of metabolites separates these drugs.

**But `negctrl_o_glycan` is not null.** O-glycan metabolism separates MMP-active from inactive drugs at
|AUROC - 0.5| = 0.13, the *same magnitude* as TCA (0.136) in the opposite direction, it **beats its
size-matched random null** (rand = 0.0099, identical to TCA's), and it partly survives coarse-MoA
stratification (0.613, p = 0.0024, Holm 0.0097 on the `auroc` statistic; 0.570, p = 0.063 on `z`) and
similarity matching (p = 3e-04). On every criterion this report and report 07 used to argue that the TCA
signal is specific, O-glycan metabolism performs as well.

**What this does to the interpretation.** The honest reading is no longer "MetChange shows a TCA-capacity
shift". It is that **MMP-active and inactive drugs differ across multiple metabolic pathways at once** - TCA
down, O-glycan up, 34 of 40 primary sets leaning below 0.5 - and that beating a size-matched random null is
**not sufficient** to establish that a particular pathway is mechanistically implicated. The specificity
argument of report 07, 17.2 ("MetChange's difference is concentrated where expression's is diffuse") survives
only in the weak sense that MetChange's difference is concentrated in *a few* pathways rather than all of
them; it does not support TCA as the pathway.

It is tempting to argue that O-glycan metabolism is not really irrelevant - its metabolites are nucleotide
sugars (UDP-GalNAc, UDP-Gal, CMP-NeuAc), which depend on UTP/CTP and on glucose flux, so an energy-compromised
cell might plausibly show altered capacity there. That may well be true, but **the set was declared irrelevant
before it was run, and explaining it away afterwards is exactly what pre-declaration exists to prevent.** It
is recorded as a failure of the negative-control test.

### 3.2 Where inside each set the signal sits — and why TCA and O-glycan are not the same phenomenon

The obvious rescue for 3.1 is that O-glycan metabolism is **not** really irrelevant: its nucleotide sugars
(UDP-GlcNAc, UDP-GalNAc, UDP-galactose, CMP-NeuAc) are made by the hexosamine pathway from
fructose-6-phosphate, glutamine, UTP and **acetyl-CoA**, and sialic acid synthesis needs PEP and acetyl-CoA, so
the set is metabolically coupled to the TCA cycle. The direction even fits: mitochondrial dysfunction drives a
glycolytic shift with carbon diverted from oxidation into biosynthesis, which would lower TCA capacity and
raise nucleotide-sugar capacity.

**That hypothesis makes a prediction, and the prediction fails.** Splitting the 30-member set by where it sits
in the network (`metabolite_sets_oglycan_split.tsv`, `attene_oglycan_split/`):

| subset of `o_glycan` | members | AUROC (z / auroc) | p | rand |
|---|---|---|---|---|
| nucleotide sugars — CMP, CMP-NeuAc, UDP-GlcNAc, UDP-GalNAc, UDP-Gal, serine (**the TCA-coupled part**) | 6 | 0.488 / 0.475 | 0.66 | 0.83 |
| the same without serine | 5 | 0.498 / 0.483 | 0.95 | 0.99 |
| **glycan cores** — core 2-8, T / Tn / sialyl-T / disialyl-T antigens (**metabolic dead ends**, degree 2-12) | 21 | **0.601 / 0.649** | **4.4e-08** | **0.005** |

The TCA linkage runs entirely through the nucleotide sugars and **those are null**; the signal is in the
terminal glycan structures. So O-glycan's association is not via its coupling to central carbon metabolism.

What it is instead is more worrying than the hypothesis it replaces: **for a low-degree dead-end metabolite,
MetChange "production capacity" is close to a readout of the few local enzymes that make it** - here
glycosyltransferases (GALNT, C1GALT1, ST3GAL, B3GNT families). On that reading the O-glycan association is a
transcriptional-program difference between phenotype-active and inactive drugs, consistent with 34 of 40
primary sets leaning one way.

**The same split applied to TCA goes the other way** (`metabolite_sets_tca_split.tsv`,
`attene_tca_split/`):

| subset of `tca mito_nohub` | members | AUROC (z) | p | rand |
|---|---|---|---|---|
| **`tca_core_cycle`** - citrate, cis-aconitate, isocitrate, AKG, succinyl-CoA, succinate, fumarate, malate, OAA (degree 8-39) | 9 | **0.367** | 1.1e-06 | 0.005 |
| `tca_peripheral` - thiamin-PP, ubiquinol / ubiquinone, acetoacetate, glyoxalate, glycolate, hydroxypyruvate … (degree 2-17) | 10 | 0.393 | 8e-05 | 0.045 |

**The TCA signal is carried by the cycle's own intermediates**, high-degree metabolites whose capacity depends
on many reactions, and it is if anything stronger there than in the set's periphery. That is the opposite
architecture to O-glycan.

**So the two results are different in kind** in where inside the set the signal sits, which partially restores
the TCA interpretation that 3.1 had removed (but see 3.3: per member the two effects are the same size, so the
difference is interpretability, not strength):

- TCA: high-degree central metabolites, signal in the cycle proper — consistent with a throughput reading;
- O-glycan: low-degree terminal structures, signal only in the dead ends, connected precursors null —
  consistent with local glycosyltransferase expression.

What 3.1 establishes still stands: beating a size-matched random null does not by itself implicate a pathway,
and a pre-declared irrelevant set did show a comparable association. What 3.2 adds is that the two
associations have different internal structure, and only one of them looks like metabolism. **A useful general
diagnostic follows: for any set that comes out significant, check whether the signal sits in its high-degree
or its low-degree members.** Low-degree members make MetChange a proxy for a handful of local genes; high-degree
members are the ones for which a capacity interpretation means anything.

### 3.3 What the O-glycan set actually contains, and a correction to 3.2

3.2 called the O-glycan association the weaker of the two. **Tested per member, that is wrong: they are the
same size.** Testing each metabolite separately removes the effect of within-set correlation:

| | median \|AUROC - 0.5\| per member | members with p < 0.05 |
|---|---|---|
| TCA core cycle | **0.100** | 8/9 |
| O-glycan cores | **0.097** | **21/21** |

TCA per member: malate 0.349 (p = 2.9e-08), citrate 0.383, fumarate 0.389, AKG 0.395, isocitrate 0.400,
OAA 0.409, succinate 0.432, cis-aconitate 0.443, **succinyl-CoA 0.473 (p = 0.32, null)**. O-glycan per member:
every one of the 21 between 0.567 and 0.614.

**So the difference is interpretability and internal structure, not effect size or robustness.**

**Collinearity, and what it does to the random-set control.** Mean pairwise Spearman between members, across
drugs:

| set | members | mean pairwise rho | median |
|---|---|---|---|
| O-glycan cores | 21 | **0.799** | 0.871 |
| TCA core cycle | 9 | 0.511 | 0.645 |

The glycan cores are near-collinear - sequential nodes on one linear route, so effectively **one degree of
freedom measured 21 times**, which is redundancy rather than independent corroboration. TCA's members are
correlated but far less, and they behave heterogeneously (succinyl-CoA null while malate is the strongest),
which is what genuine differential behaviour along a cycle looks like.

**This is a caveat on a control used throughout reports 07 and 08.** A size-matched random set draws members
that are weakly correlated with each other, so its *average* is noisy; a near-collinear pathway set averages
near-copies, so its average is clean. **Beating a size-matched random null is therefore partly a function of
within-set correlation, not only of biological specificity** - the control implicitly assumes comparable
internal correlation and that assumption does not hold here. Wherever `rand_p` is quoted in reports 07-08 it
should be read with that in mind, and **set-level statistics should be accompanied by per-member ones**, which
are immune to it. A better null would match random sets on effective degrees of freedom rather than member
count.

**What the 21 compounds are.** Mucin-type O-glycosylation structures: the Thomsen-Friedenreich antigen family
(**Tn-antigen, T-antigen, sialyl-T, sialyl-Tn, disialyl-T**), the **core 2-8** classes,
`glcnac-alpha-1,4-core 1/2`, and the KEGG glycan identifiers `F1alpha`, `G00031`, `G00032`. Note that most are
**lumped glycan classes rather than defined molecules** - "core 4" is a family, not a compound one could
measure - and **four are not compounds at all**: `[protein]`, `[protein]-L-serine` (cytosolic and Golgi) and
`Protein-Linked Serine Or Threonine Residue (O-Glycosylation Site)`, which are GEM bookkeeping entities for
the acceptor protein. They carry the same AUROC as the real glycans (0.567-0.571), which is the clearest
indication that this set reports enzyme expression rather than metabolite capacity.

**Does this reach the model? No, and only partly for the reason I first assumed.** 291 of the 5,707 MetChange
features (5.1%) are placeholders or lumped pools (`PC-LD pool`, `PI pool`, `SM pool`, `ceramide pool`,
`fatty acid pool`, `acyl-CoA pool`, `[protein]-L-serine`, ...). They are in the label matrix because MetChange
is computed for every metabolite of the model and limma was run over the whole matrix. **They are removed
during network construction** - together with the currency metabolites - because they are hubs connecting
everything and carry no real biology, so they cannot be selected as model sources or used in a query
(consistent with report 07, section 5: 3,747 of 4,086 features are exact network nodes, the 339 dropped
including the currency metabolites).

Checked against `data/networks/combined/nodes.tsv` for this particular set: the three true `[protein]`
entities are indeed **absent** from the network. But 18 of the 21 remain metabolite nodes, including
`Protein-Linked Serine Or Threonine Residue (O-Glycosylation Site)` and all the lumped glycan classes. So the
filter removes the currency metabolites and the `[protein]` placeholders; it does not remove lumped
pseudo-species. They are low-degree (2-12) so they will not act as hubs, but they can in principle be selected
as sources when they are significant and large. **Worth one cheap check: how often pool / lumped names appear
among the 100 selected sources across perturbations.**

**Status of the O-glycan finding.** A real group difference, the same magnitude as the TCA one and more
internally consistent, but with roughly one effective degree of freedom and carried by nodes that are largely
lumped classes and modelling entities. The most parsimonious reading is altered **glycosyltransferase
expression** (GALNT, C1GALT1, ST3GAL, B3GNT families) reported through terminal glycan nodes - a
transcriptional observation, reinforced by the fact that the nucleotide-sugar precursors are null, so substrate
supply is unchanged (3.2). It is not usable as a negative control, it is not dismissible as weak, and it is not
evidence for the metabolic claim. The test that would settle it is the expression arm on an O-glycan GPR gene
set: if the genes separate the groups at least as well, the finding is transcriptional and MetChange adds
nothing there.

## 4. Chemical-similarity matching: ran, uninformative

ECFP4 (Morgan radius 2, 2048 bits) on `canonical_smiles`, Tanimoto, greedy 1:1 without replacement, cases
ordered hardest-first, tested **paired** (Wilcoxon signed-rank on within-pair differences; the first version
used an unpaired test and so discarded the pairing).

| spec | fraction of pairs with case lower | p | median Tanimoto |
|---|---|---|---|
| `tca mito_nohub` z | 0.673 | 6.8e-07 | **0.268** |
| `tca_endpoints cyto_nohub` z | 0.660 | 2.0e-06 | 0.268 |
| `energy_currency mito` z (MMP-up) | 0.307 (i.e. 69% higher) | 0.0059 | 0.300 |

The effects hold, but **median Tanimoto 0.27–0.30 means the matched pairs are structurally unrelated**, so
this controls for almost nothing. It is reported for completeness and should not be cited as a chemotype
control. A useful version needs a similarity floor (e.g. only pairs above 0.6), which this library cannot
supply at n = 147.

## 5. Potency and dose exposure: no gradient, and the design cannot provide one

Within actives only (controls have no potency by construction), Spearman of the per-drug score against
log10 of the measured `Ratio Potency (uM)` (the most potent record per compound).

| spec | rho | p | Holm |
|---|---|---|---|
| `tca cyto` sign_dev | +0.205 | 0.013 | 0.076 |
| **`oxphos mito` z (power diagnostic)** | **+0.191** | **0.020** | 0.102 |
| `tca_endpoints cyto_nohub` z | +0.080 | 0.34 | 1.0 |
| `energy_currency mito` z (MMP-up) | +0.116 | 0.32 | 1.0 |
| `tca mito_nohub` z | +0.007 | 0.94 | 1.0 |
| `adenylates all` z (MMP-up) | -0.030 | 0.80 | 1.0 |

**Nothing survives correction**, and decisively: the one nominal gradient (`tca cyto sign_dev`) is matched by
a gradient of the same size in the **negative-control set**, which has no group effect at all. So potency
correlates with something non-specific, and **there is no dose-response evidence supporting the association**.

A graded relationship between measured potency and the label would have been hard for a chemotype to produce,
and it is the one test here that could have spoken *for* the result rather than only guarding against an
artefact. It is absent.

### 5.1 What potency is here, and why the range is the problem

`Ratio Potency (uM)` is the **AC50** of the ratiometric MMP readout: the concentration at half of that
compound's own fitted maximal response, from a Hill fit to the qHTS titration (every compound tested across
~14-15 concentrations up to ~92 uM, which is why potency exists for ~10,000 compounds). `Ratio Efficacy (%)`
is the fitted maximum, median **-92%** - most actives collapse the potential almost completely, so potency
encodes *where on the dose axis* rather than *how much* biology. A *decade* below means a factor of 10 on a
log scale (one order of magnitude).

| | values | ratio | decades |
|---|---|---|---|
| AC50 across actives, full range | 0.0096 - 72 uM | 7,500x | **3.9** |
| AC50, interquartile | 5.05 - 28.5 uM | 5.6x | **0.75** |
| max LINCS dose, interquartile | 10 - 10 uM | 1x | **0.00** |

So the middle half of the actives sits inside a single 10-fold window, and **the LINCS profiling dose does not
vary at all** (39% of the relevant perturbations are at 10 uM). Add the censoring - only actives have an AC50,
because a compound whose curve does not fit below the top concentration is *called inactive* - and the weakest
half of the potency distribution is structurally missing from the correlation.

### 5.2 Dose exposure: 48% of the cases were profiled below their own AC50

The phenotype call and the transcriptional profile are different experiments at different concentrations.
Comparing each case drug's AC50 with the highest dose LINCS actually profiled it at:

| MMP-down drugs (n = 147) | n | % |
|---|---|---|
| AC50 <= max LINCS dose (profiled at or above its potency) | 77 | 52% |
| **AC50 > max LINCS dose (profiled below its potency)** | **70** | **48%** |

Median log10(dose / AC50) = +0.04: the profiling dose sits essentially *at* the AC50, with half the compounds
below it. **So about half the case group was transcriptionally profiled at a concentration where it does not
measurably collapse the membrane potential.** That is a real property of the design and worth knowing
independently of any test.

`dose_exposure_check.py` splits the cases on it. Prediction if the association is a dose-dependent metabolic
response: the exposed half should show a stronger effect (i.e. a *smaller* AUROC) and the sub-threshold half
should resemble inactives.

| set / statistic | all cases | exposed (77) | sub-threshold (70) | exposed vs sub | rho(log10 dose/AC50) |
|---|---|---|---|---|---|
| `tca mito_nohub` z | 0.364 | **0.361** | **0.367** | p = 0.68 | +0.002, p = 0.98 |
| `tca_endpoints cyto_nohub` z | 0.359 | 0.336 | 0.384 | p = 0.61 | -0.066, p = 0.43 |
| `tca cyto` sign_dev | 0.415 | 0.385 | 0.449 | p = 0.20 | **-0.193, p = 0.019** |
| `oxphos mito` z **(power diagnostic)** | 0.476 | 0.444 | 0.510 | p = 0.21 | **-0.217, p = 0.0082** |

**The prediction fails for the primary set**: exposed and sub-threshold are identical (0.361 vs 0.367) with no
continuous relationship whatever. Where a relationship does appear (`tca cyto sign_dev`), the **negative
control shows a stronger one**, so it is non-specific - the same pattern as the potency correlations in the
table above.

### 5.3 Why this test is weaker than it looked when proposed

Two limitations, both of which should have been anticipated before running it:

1. **It is not independent of the potency test.** Because the maximum LINCS dose has an interquartile spread of
   0.00 decades while AC50 spans 0.86, log10(dose / AC50) is driven almost entirely by AC50. The test is
   therefore close to the potency correlation dichotomised at 10 uM - which is exactly why it reproduces its
   results, `oxphos mito` included.
2. **The MMP AC50 is not a threshold for a metabolic effect.** It comes from a 1-5 h ratiometric dye assay,
   while the LINCS signature is 6-24 h transcription, and a drug may perturb metabolism well below the
   concentration at which the potential fully collapses. "Profiled below its AC50" therefore does not mean "no
   metabolic effect was there to detect", so the failed prediction does not carry the weight it would in a
   matched-exposure design.

Minor: exposure is assigned per drug from its highest profiled dose. 71% of these drugs have exactly one
perturbation so the assignment is exact for most; for the 28% profiled at 4-7 doses the averaged score mixes
exposed and unexposed perturbations and dilutes the contrast.

**Conclusion for this line of attack: exhausted, and not by sample size.** LINCS profiled essentially
everything at one concentration, so the data cannot separate potency from compound identity, and the MMP AC50
is the wrong threshold for the quantity of interest. No further dose or potency analysis of this dataset will
settle the question; a design where exposure is matched to the phenotype concentration would be needed.

## 6. Where the result stands

**Survives:** responsiveness matching (report 07, 16; unchanged by the paired test added here —
`tca mito_nohub` 0.636 of pairs, p = 0.0023), the size-matched random-set null, the matched-universe-fraction
null-shape comparison (report 07, 17), and coarse mechanism-family stratification (section 3).

**Seriously qualified by the declared negative controls (section 3.1):** one of five
biologically-irrelevant pathways shows an association of the same magnitude, in the opposite direction, which
also beats its size-matched null and also survives coarse-MoA stratification. **Partially restored by 3.2:**
that set's signal lives entirely in its metabolic dead ends (its TCA-coupled nucleotide sugars are null),
whereas the TCA signal lives in the cycle's own high-degree intermediates - so the two are different in kind,
and only one is structured like a metabolic effect.

**Does not survive / absent:** any potency gradient and any dose-exposure gradient (section 5), with the one
nominal appearance of each matched or exceeded by `oxphos mito`; the ToxCast compound subset
(section 2); narrow-MoA stratification, which is untestable at this n rather than negative (section 3).

**Neither supports nor refutes:** the dose-exposure split, once its two limitations are taken into account
(5.3). It is reported because the 48 / 52 exposure imbalance is a fact about the design that anyone using this
dataset should know, not because the test discriminates.

**Defensible statement:**

> In HepG2, the MetChange labels of drug-like LINCS compounds differ between compounds measured as
> mitochondrial-membrane-potential active and inactive. The difference is multi-pathway: TCA capacity lower
> (AUROC 0.38 within coarse mechanism families, Holm p = 0.005) and O-glycan capacity higher (0.61,
> Holm p = 0.010) - the latter a pathway declared irrelevant in advance (3.1). The transcriptome read over the
> same GEM subsystems does not reproduce the TCA difference (report 07, 16). The effect is confined to the
> drug-like part of the library, shows no dose-response relationship with measured potency or with whether the
> compound was profiled above its own AC50, and has not been tested within narrow mechanism classes. Note that
> this dataset cannot provide a dose-response test: LINCS profiled essentially every compound at 10 uM (5.3).

**Not defensible:** that MetChange predicts mitochondrial toxicity; that it models per-drug metabolic change
(the per-drug test in report 05, 8.2 is negative); that the finding generalises to environmental chemicals
(section 2 says it does not); or that the association is **uniquely** about the TCA
cycle, since a pre-declared irrelevant pathway shows an equal and opposite association passing the same tests
(3.1) - though 3.2 shows the two have different internal structure and only the TCA one is carried by
high-degree central metabolites.

## 7. Next

1. **LINCS genetic perturbations of GEM enzymes** (`trt_sh` / `trt_xpr` / `trt_oe`) — now clearly the decisive
   experiment, for a reason this report sharpens: **a knockdown has no chemotype and no MoA**, so every
   confounder in sections 2–5 is unavailable by construction. Readouts as in STATUS 6.2b item 1: does
   knockdown of an enzyme read as less production of its products, with the 1–2-reaction-away version
   reported separately from the (partly circular) 0-reaction one.
2. **A second cell line** for the confounder structure, once correct A549 labels exist.
3. **`APR_HepG2_MitoMembPot` 24h/72h** with `MitoMass` as a specificity control (orthogonal MMP readout, same
   cell line, 1,108 chemicals, in `RawMC5_ToxCast_by_aeid.csv`, mappable via `casn`). Note it is also ToxCast,
   so by section 2 it should come out **null** — which makes it a test of section 2's reading rather than of
   the association.
4. A similarity-matched analysis with a Tanimoto floor, if a library with closer analogue pairs can be
   assembled (section 4).
5. Not a priority: more phenotype datasets of the same kind, and no further potency or dose analysis of this
   one (5.3). Sections 2–5 show the limit here is confounding, compound-class scope and single-concentration
   profiling, not sample size.
6. If a dose-response test is ever wanted, it needs the exposure matched to the phenotype concentration -
   e.g. LINCS perturbations of the same compound at several doses spanning its AC50, which exists for the 28%
   of these drugs profiled at 4-7 doses and could be analysed per perturbation rather than per drug (the
   per-column statistics are already in `set_member_values.npz`).

## 8. Running it: the launcher chain after the 2026-11-20 refactor

The cell line is now a parameter of the outermost launcher, and nothing below it names a cell line.

```
sbatch launch_assoc_HEPG2_slurm.sh   (or launch_assoc_A549_slurm.sh)
   └── export cell=<CELL>; exec launch_assoc_slurm_generic.sh   # --array=0-1
          ├── task 0 → launch_assoc_metchange.sh      # MetChange arm
          └── task 1 → launch_assoc_expr.sh           # expression baseline
                 └── both source mc_phenotype_assoc_common.sh
                        └── mc_phenotype_assoc.py / stratified_check.py / matched_check.py
```

To add a cell line, copy a slurm wrapper and change one line; the per-cell wrappers differ only in
`export cell=`. Run an arm directly with `cell=A549 bash launch_assoc_metchange.sh`, or one modality with
`cell=A549 modes="--mode sign" bash launch_assoc_metchange.sh`.

### 8.1 Where things live, and why

`mc_phenotype_assoc_common.sh` resolves `scratch_proj_dir` to
`/scratch/users/gsantamaria/test_large_files/metabognn_dev` when it exists and to the repo otherwise, so
nothing large lands on the home partition:

| what | path | scope |
|---|---|---|
| labels | `<scratch_proj_dir>/results/labels/<cell>/{metchange/rank,expr}` | per cell |
| outputs | `<scratch_proj_dir>/results/mc_phenotype_assoc/<cell>/` | per cell |
| **phenotype tables** | `<scratch_proj_dir>/results/mc_phenotype_assoc/<cell>/phenotypes/` | **per cell** |
| **metabolite / gene sets** | `<scratch_proj_dir>/results/mc_phenotype_assoc/metabolite_sets/` | **per project** |

**Phenotype tables are per cell line** because they are keyed on LINCS perturbation names and built from that
cell line's structure mapping — the mappable compound set differs, so sharing them would silently drop
compounds. That failure mode produces a plausible wrong answer rather than an error, so `require_inputs`
hard-fails on a missing table instead of falling back, and every run prints a provenance block naming the
cell, the label dirs, the phenotype dir and the resampling counts.

**Sets are per project** because they are built from Human-GEM alone. To make that correct rather than merely
convenient they are now built **without `--features`**, i.e. not restricted to any one label matrix;
`mc_phenotype_assoc.py` drops members absent from the matrix it is given and **reports what it found** (e.g.
`oxphos mito  13 of 16 members found (81%)`), and says when a set falls below `--min-set-size` for that cell.
Previously `--features` baked one cell line's feature list into the sets, which is what made them
cell-coupled. Verified equivalent: on HEPG2 the new layout scores the same 40 sets and reproduces the previous
per-set numbers exactly (`glycolysis cyto_nohub` z 0.379 / auroc 0.406).

### 8.2 One-off preparation per cell line

```
cell=<CELL> bash scripts/mc_phenotype_assoc/prepare_assoc_inputs.sh
```

writes the project-level sets (if absent; `rebuild_sets=true` forces) and the per-cell phenotype tables,
including the by-mechanism Hallinger split. It needs the structure mapping for that cell line, and says so
and stops if it is missing:

```
Rscript scripts/parse_lincs_data/get_lincs_pert_info.R          # -> lincs_unique_perts_<CELL>.tsv
cell=<CELL> bash scripts/benchmarking/map_attene_lincs_perts.sh      # -> mapping_<CELL>/
cell=<CELL> bash scripts/benchmarking/map_hallinger_lincs_perts.sh   # -> mapping_<CELL>/
```

Both mapping scripts now take `cell` and write to `mapping_<cell>/`, falling back to the original HEPG2
locations when those are what exist. `build_phenotype_tables.py` takes `--cell` and a required `--out-dir`.

**Read the A549 run for what it is.** The Attene-Ramos and Hallinger phenotypes were measured in **HepG2**
whatever cell line is mapped. Running them against A549 labels is the **cross-cell-line control** of
section 7 item 2, not a replication: if the association is as strong with A549 labels against a
HepG2-measured phenotype, it is a compound property rather than HepG2 metabolism. The mapping scripts say so
in their headers.

### 8.3 Superseded locations

`data/validation/metabolite_sets/` and `data/validation/phenotypes/` were the locations used for every result
in reports 07 and 08 and are referenced as such there. They are left in place so those references resolve;
new runs read the paths in 8.1.
