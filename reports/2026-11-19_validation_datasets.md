# Report 07 — Validation datasets: compilation, structure mapping, and what each can actually test

*2026-11-19. No model runs. Compilation and triage of external validation data, the rewrite of
`scripts/benchmarking/map_smiles.R`, and the chemical-identification work for Hallinger (Seahorse) and
Attene-Ramos (mitochondrial membrane potential). Design for using Attene-Ramos against the HepG2 MetChange
labels. Continues report 06 section 7 (Seahorse query design) and STATUS 6.2b.*

## 1. The datasets and what each can test

| dataset | readout | cell line | compounds | mapped to LINCS | verdict |
|---|---|---|---|---|---|
| Schuhknecht 2025 | intracellular metabolomics (pools) | A549 | 1,520 | 750 | **done, null** (report 05, 8.2) |
| Sarvin 2026 | metabolomics, 6 doses x 6 times | HCT116 | 76 (30 enzymes) | ~22-30 | raw mzXML only; enzyme ground truth, small |
| **Attene-Ramos 2014** | **mito membrane potential, qHTS** | **HepG2** | **~10,000** | **904** | **best powered; HepG2 labels exist** |
| Hallinger 2020 | Seahorse OCR (flux) | HepG2 | 1,042 tier 1 / 270 tier 2 | 42 of 258 named (17%) | only tier 2 identifiable (section 4) |
| Lackner 2022 | isotope labelling | breast cancer line | 1 drug class | — | case study, not a benchmark |
| TraceBase / Metab. Workbench | time-resolved metabolomics | various | scattered | — | harmonisation cost; later |
| LINCS shRNA / CRISPR / OE | genetic perturbation of enzymes | many | — | — | **the one test of the labels** (STATUS 6.2b item 1) |

Flux readouts (Hallinger, Lackner) test the production / consumption claim that pools cannot (report 05, 9.3);
pool and potential readouts are characterisation. Power, however, runs the other way: Attene-Ramos is the only
one with enough mapped compounds for a per-drug statistic.

## 2. `map_smiles.R` rewritten

The original (`intersect` + `match` on rcdk canonical SMILES) was sound but lost matches four ways. Fixed:

1. **Salts / mixtures reduced to the largest fragment** before canonicalising, so a hydrochloride matches the free
   base. The biggest source of loss: 2,040 of Attene-Ramos's 10,337 SMILES carry a `.` (LINCS: 143 of 34,419).
2. **All pairs kept** (merge on the canonical structure instead of `match`, which returned only the first hit on
   each side and silently dropped legitimate one-to-many links).
3. **NA can no longer match NA** (R's `intersect` keeps `NA` when both sides have it).
4. **Unique SMILES canonicalised once**, which also makes it much faster on tables that repeat rows per well.

**Stereochemistry is no longer decided silently.** `--stereo both` (default) matches exactly first, then retries
only the unmatched compounds without stereochemistry; every pair is labelled `exact` / `stereo_insensitive` and
carries `n_pert_forms` / `n_val_forms` (how many stereo forms share that flat structure per side) and
`ambiguous_stereo`. So the cost of merging enantiomers is measurable per analysis instead of assumed.
`--stereo exact|flat` force one behaviour. A `mapping_summary.txt` records parsed / matched / by type.

## 3. Attene-Ramos mapped to LINCS HepG2

`map_smiles.R` on `lincs_unique_perts_HEPG2.tsv` vs `AID_720637_datatable.csv`: **1,509 pairs, 904 LINCS
compounds, 1,379 Attene-Ramos CIDs**, of which 718 pairs exact and 791 stereo-insensitive (225 flagged
ambiguous).

**The 580 "failed to parse" perturbations are not a parsing problem:** 564 of them have **no SMILES and no
InChIKey** in LINCS at all (76% are unnamed `BRD-…` codes; in the full LINCS file 5,785 of 34,419 are like this).
They are unmappable to any external dataset by structure, permanently. So the real denominator is **4,930
mappable HepG2 compounds**, and 904 is 18% of those, not 16%. The remaining ~16 are genuinely unparseable.
(The script's summary should distinguish "no structure recorded" from "unparseable"; small fix, not done.)

Two decisions before using the mapping: (a) 1,379 CIDs collapse onto 904 compounds (salts merging, repeated
CIDs), so the Attene-Ramos measurements need an aggregation rule per LINCS compound (mean, or most potent);
(b) run first with `exact` + unambiguous pairs, then check whether adding the 225 ambiguous ones changes
anything.

## 4. Hallinger: only the 270 tier-2 actives can be identified

**The tiers** (Dryad readme): tier 1 = single-concentration screen of **1,042** chemicals (`sc_seahorse`, 1,051
spids); tier 2 = multi-concentration follow-up of the **249 actives** (`mc_seahorse`, 273 spids); tier 3 = 149
putative ETC inhibitors (`EFA`).

**Only tier 2 is named.** All three annotation files in the deposit (`RawMC3`, `RawMC5`,
`mc5_mc6_ncct_mito_nov2019`) cover the same 270 of 1,051 spids. The rest are EPA plate IDs only
(`EPAPLT0033A02`…); `cpid` is just the source plate (two plates), so the deposit cannot identify them.

**Everything public was checked and none of it has the plate maps:**

| source | result |
|---|---|
| `analytical_qc_invitrodb_v4_3` (2.9 MB) | has `spid` + `dsstox_substance_id` for 33,690 samples, 97 EPAPLT plates — **but not 0033 / 0034** |
| `assay_annotations_invitrodb_v4_3` | no chemical / sample columns at all |
| `pubchem_invitrodb_v4_3` (per-aeid) | aeids 2442 / 2444 / 2446 exist, but their Results sheets hold the same 272 tier-2 spids |
| **`sc1_sc2_invitrodb_v4_3` (from the 7.5 GB summary zip)** | **0 of 1,051 spids; no `CCTE_Simmons_MITO` rows; aeids 2442/2444/2446/2450 absent from sc2; no `EPAPLT0033/0034` sample anywhere** |

Only four `CCTE_Simmons_MITO` endpoints are registered in invitrodb (2442, 2444, 2446 multi-concentration,
2450 viability), all HepG2. **The tier-1 screen was never registered**, so no public invitrodb release contains
it. The paper's supplementary tables give chemical names with activity calls but no plate or sample IDs.

**Remaining options:** (a) ask the corresponding author (Simmons, EPA/CCTE) for the plate maps of `EPAPLT0033`
and `EPAPLT0034` — now the only real source; (b) the `sample` table of the 17 GB MySQL dump, which registers
samples independently of assay data, but given tier 1 was never registered the odds are low and it is not worth
17 GB on spec.

**What the 270 give**: 266 resolved to structures via the DSSTox dump (6 failures are mixtures: abamectin,
milbemectin, a surfactant), 42 match LINCS by InChIKey first block (17%), 29 of those have annotated targets,
11 are also in A549 LINCS. The unmatched are pesticides and industrial chemicals (allethrin, buprofezin,
cyazofamid, propargite, metalaxyl): ToxCast is an environmental-chemical library, LINCS a drug library, and they
barely intersect. Extrapolating 17% to 1,042 would give ~170 LINCS compounds even if the identification problem
were solved. **Hallinger is therefore tier 2 of the validation list, not tier 1.**

**Tools acquired on the way:** `DSSTox_CCD_dump_12092025` (CSV dump,
`DTXSID / CASRN / INCHIKEY / SMILES / QSAR_READY_SMILES`; a dated release, so identifier → structure is now
reproducible offline, which the Dashboard batch search is not) and
`data/validation/hallinger/hallinger_chems_dsstox.csv` (the 272 named chemicals with structures). Clowder files
download programmatically as `https://clowder.edap-cluster.com/api/files/<file-id>/blob` — note **HEAD returns
404 on that endpoint, only GET works**, so link checkers wrongly report it broken.

## 5. Currency metabolites: in the labels, not in the network

Checked because an energy-metabolism query needs them. **They are MetChange features** — ATP (7 compartments),
ADP (7), NAD+ (5), NADH (4), NADPH (5), NADP+ (4), O2 (7), H+ (8), H2O (7), CoA (6) — but they are among the 339
features that are **not network nodes** and are dropped before `select_sources` (report 05, 5.5). Ubiquinol has
one node and no MetChange feature; ubiquinone neither. Consequences:

- **label-level tests can use the energy readout** (NADH / NAD+, ATP / ADP, O2), which is exactly what a
  membrane-potential or respiration phenotype predicts;
- **model queries cannot**: a query can only be built from network nodes, so the glycolysis / TCA proxies
  (lactate 8 features, pyruvate 8, succinate 10, citrate 7, fumarate 3, acetyl-CoA 7, glucose 19) carry it;
- therefore **the enrichment test is richer than the ranking test** for these phenotypes, which reinforces doing
  the label-level analysis first (report 06, section 7).

## 6. Design: Attene-Ramos against the HepG2 MetChange labels

**Phenotype.** A drop in mitochondrial membrane potential means a compromised ETC (direct inhibition or
uncoupling). The expected metabolic consequence is a shift from oxidative to glycolytic metabolism: lactate and
pyruvate production up, TCA intermediates down, glucose consumption up, NADH / NAD+ shifted — and all of those
are MetChange features (section 5).

**Test, one value per drug** (so ~900 drugs, properly powered, unlike the Schuhknecht correlation): summarise each
drug's MetChange signature over the glycolytic and the TCA / ETC metabolite sets, then compare MMP-active with
MMP-inactive drugs (Mann-Whitney / AUROC).

**Three controls decide whether it means anything:**
1. **Separate mitochondrial toxicity from plain cytotoxicity.** The qHTS MMP assay reads viability in the same
   wells, and most compounds that collapse MMP are simply killing cells. The informative contrast is
   **MMP-down without viability loss** vs inactive. Without it the test mostly predicts cytotoxicity. *Check the
   viability columns of `AID_720637_datatable.csv` first.*
2. **Expression baseline** read the same way, as in the Schuhknecht analysis: does the network reading add
   anything over the genes?
3. **Do not fix the expected sign in advance.** The compensation hypothesis (report 05, 9.2) applies: a drug
   blocking the ETC post-translationally may upregulate OXPHOS genes, so MetChange could read "more TCA
   production" while flux falls. Report the direction found.

**Prerequisite, and it is cheap:** the HepG2 MetChange labels exist on the cluster
(`<scratch>/results/labels/HEPG2/metchange/rank/`, ~1.5 GB per matrix), not on the laptop — so this analysis runs
on the HPC, like the Schuhknecht association. The usable n is not 904 but the intersection of the mapped
compounds with **active HepG2 perturbations**, which is the first number to compute there.

## 7. Next

1. Attene-Ramos viability columns → is the cytotoxicity split possible? (laptop, now)
2. Mapped compounds x active HepG2 perturbations → the real n. (HPC)
3. The enrichment analysis of section 6, reusing `scripts/schuhknecht_mc_assoc/` machinery. (HPC)
4. Hallinger: email Simmons for the two plate maps; use the 270 meanwhile.
5. Unchanged from STATUS 6.2b: the LINCS genetic-perturbation test remains the one test of the labels.

## 8. Hallinger: the plate map is not needed (added 2026-11-19, later)

**Section 4 above is superseded in its conclusion, not in its facts.** The plate maps of `EPAPLT0033` /
`0034` / `0035` are indeed not public — but identifying the wells was never the only route to the tier-1
chemicals. **Supplemental Table 2 of the paper lists all 1,042 tier-1 chemicals by name with their activity
calls**, and a second supplemental table characterises the 243 actives. Only the *per-well raw OCR values*
need the plate map; a per-compound active / inactive contrast does not.

**Where the supplement is.** The journal supplement is deposited as EPA figshare
`doi:10.23645/epacomptox.12123735` ("Supplemental files for Respirometric Screening…"), which holds no files
itself: it is a single link-only entry pointing at an EPA Clowder folder in dataset
`61147fefe4b0856fdc65639b` (the same dataset the DSSTox dump came from). The files, by Clowder file id,
fetched with the `api/files/<id>/blob` GET trick of section 4:

| file | id | content |
|---|---|---|
| `SuppTable2_SC.summ.table.pdf` | `66969672e4b0a7c65d1a38a1` | 19 pp: **all 1,042 tier-1 chemicals**, test concentration, active / inactive |
| `SuppTable2_Study_summary.pdf` | `668451dae4b0a7c65d175744` | 26 pp: the **243 actives**, four per-aeid calls, proposed mechanism, electron-flow-assay activity, inhibited ETC complex |
| `SuppTable1_QC_summary.pdf` | `66969672e4b0a7c65d1a3899` | plate QC metrics per tier |

Copies in `data/validation/hallinger/supplemental/`.

**Parsed** by `scripts/benchmarking/parse_hallinger_supp.py` (needs `pypdf`; re-runnable, every count below is
printed by it). Both tables parse cleanly, with no unmatched data lines:

- `hallinger_tier1_suppT2.csv` — 1,042 chemicals, **243 active / 799 inactive**, exactly the paper's counts.
- `hallinger_tier2_mechanism.csv` — the 243 actives. Proposed mechanism: **ETCi 149, ATP synthase 29,
  uncoupler 15, redox 10, interference 1, "inactive" at tier 2 39**. Inhibited complex where the electron
  flow assay resolved it: **I 59, III 15, II 5, mixed I/III 4, IV 1**, `nd` 65. The per-aeid columns are
  2442 (decreased basal), 2444 (decreased maximal), 2443 (increased basal), 2447 (increased inhibited) —
  note 2443 / 2447 appear here but section 4 found only 2442 / 2444 / 2446 / 2450 registered in invitrodb.
- `hallinger_tier1_all1042_annotated.csv` — the two joined, plus DSSTox structures and the LINCS HepG2 match.

**Yield.** 960 of 1,042 names resolve to a structure through the dated DSSTox CCD dump (`PREFERRED_NAME` plus
the `IDENTIFIER` synonym column). Matching the first InChIKey block against `lincs_unique_perts_HEPG2.tsv`
gives **98 LINCS HepG2 compounds: 23 active, 75 inactive** — against 42 actives-only from the tier-2 route.
That is a floor on three counts: 82 names are still unresolved (`Disulfiram`, `Niclosamide`, `MEHP`,
`Pirinixic acid`, `4-tert-Butylphenol` … — all real drugs or well-known chemicals, so synonym or
PubChem resolution should recover most), the match ignores salts and stereochemistry, and `map_smiles.R`
has not been run on the `SMILES` column yet.

**What changes.** The tier-2 route gave actives only, so there was no contrast to test. There is now one:
23 vs 75 at worst, probably ~30 vs ~100 after proper mapping — enough for the label-level enrichment test
(one value per drug, Mann-Whitney), which is the test that was planned for Attene-Ramos anyway and can now
be run on both datasets with the same code. Drug *ranking* is still far too small here.

**And it is mechanistically labelled, which Attene-Ramos is not.** Membrane potential is one phenotype;
Hallinger separates ETC inhibition from uncoupling from ATP-synthase inhibition, and names the complex.
Those predict *different* MetChange patterns (an uncoupler raises respiration while collapsing potential),
so the handful of LINCS-matched actives — rotenone, troglitazone, triclosan, pentachlorophenol,
phenolphthalein, 4-hexylresorcinol, butylparaben, nilutamide, diethylstilbestrol, indomethacin — are the
positive controls for reading the Attene-Ramos result, even at n ~ 10. Rotenone (complex I, EFA-active) is
the textbook case.

**Still worth asking Simmons for the plate maps** (section 7, item 4), but it is no longer blocking: it would
only add per-well raw OCR and concentration-response for the ~780 unnamed tier-1 wells, and tier 2 already
carries the quantitative curves for every active.

## 9. The real Hallinger n, and two bugs in `map_smiles.R` (added 2026-11-19, later)

### 9.1 Name resolution: 1,040 of 1,042

Section 8 left 82 names unresolved. The cause was in my own reader, not in the data: `IDENTIFIER` in the
DSSTox CCD dump is a **pipe-separated synonym list** (CASRNs and names), and the first pass hashed the whole
string as one name instead of splitting it. Splitting on `|` takes resolution from 960 to **1,040 of 1,042**
(80 via a synonym, e.g. Disulfiram = DSSTox `Tetraethylthiuram disulfide`, Niclosamide, MEHP, Pirinixic acid).
The two left are `6:2 Fluorotelomer alcohol` and `N,N,N',N'-Tetrakis(2-hydroxypropyl)ethylenediamine`, both
tier-1 inactive and neither a drug. PDF text also needed a mojibake repair (`â€²` → `'`).

### 9.2 Two bugs in `map_smiles.R`, both found by checking a pair by hand

The rewrite of 2026-11-19 (section 2) was mapping the right structures to the **wrong compounds**. The first
output paired LINCS `1-NITRONAPHTHALENE` with Hallinger `Diazinon` and printed a DDT SMILES on the row.

**Bug 1 — subset row numbers used as table row numbers.** In `pairs_on`,
`prt[!is.na(get(key)), .(pert_row = .I, k = get(key))]` returns `.I` numbered **1..n over the filtered
subset**, not the rows of `prt` (verified: `data.table(x=c(NA,"a","b",NA,"c"))[!is.na(x), .I]` gives 1,2,3,
while `[, .I[!is.na(x)]]` gives 2,3,5). Those positions were then used to index the full table
(`prt$pertid_dataset[p$pert_row]`), so every id and name was read from the wrong row, shifted by however many
unparseable SMILES preceded it — 580 on the LINCS side. The matched structures were right; the identities
attached to them were not. Fixed by taking the row numbers before filtering.

**Bug 2 — explicit hydrogens defeat the canonicalisation.** DSSTox writes stereocentres as
`[H][C@@]12COC3=...` and LINCS as `COc1cc2OC[C@H]3...`. With explicit `[H]` kept as graph nodes, CDK
canonicalises the two to **different strings**, in both the exact and the flat pass, so true matches were
silently lost — including rotenone and 17beta-estradiol, which are byte-identical by InChIKey across the two
sources. Fixed with `remove.hydrogens()` before writing the canonical SMILES. Checked that it does not merge
what should stay apart: both notations of D- and L-alanine still differ under `--stereo exact` and still merge
under `flat`.

Bug 1 affected the identity of every pair; bug 2 cost matches. **Both corrupt the Attene-Ramos mapping of
section 3** (1,509 pairs / 904 compounds), which was produced by the same script before either fix: the pair
*count* was roughly right, since matching keys off canonical SMILES was never the broken part, but the
compound identities were wrong. It has been re-run; the previous output is kept as
`data/validation/attene_ramos/mapped_drugs_BUGGY_2026-11-19.csv.bak`. **Any analysis built on the earlier
Attene-Ramos mapping must be redone.**

### 9.3 Hallinger mapped to LINCS HepG2: 115 compounds

`scripts/benchmarking/map_hallinger_lincs_perts.sh` (new; `--stereo both`, name column `chnm`), on the
`SMILES` column of `hallinger_tier1_all1042_annotated.csv`. Output in
`data/validation/hallinger/mapping/`:

| | |
|---|---|
| Hallinger tier-1 chemicals mapped | **115** of 1,041 parsed (102 exact, 13 stereo-insensitive) |
| LINCS HepG2 perturbations | 120 |
| **tier-1 active / inactive** | **29 / 86** (28 / 83 dropping the 9 ambiguous-stereo pairs) |

Progression, same data, as the bugs came out: 42 (tier-2 route, actives only) → 94 (buggy mapper) →
**115**. The InChIKey-first-block floor of section 8 gave 105, so the structure mapping is now strictly
better than the cruder key match, as it should be.

**The 29 actives carry their mechanism**, which is what makes this set useful at this n:

- **ETC inhibitors, 19** — rotenone (complex I, EFA-confirmed), troglitazone, nilutamide, phenolphthalein,
  butylparaben, 4-hexylresorcinol, diethylstilbestrol, corticosterone (all complex I), plus bisphenol A,
  fenofibrate, kepone, DMBA, tamoxifen, clomifene, progesterone, cyproterone, androstenedione and both
  estradiols (complex not resolved).
- **Uncouplers, 3** — niclosamide, triclosan, pentachlorophenol.
- **ATP synthase, 4** — indomethacin, MEHP, triiodothyronine, testosterone propionate.
- **Tier-2 inactive, 3** — disulfiram, pirinixic acid, fipronil (tier-1 active, not confirmed at tier 2).

Uncouplers and ETC inhibitors predict **opposite** respiration changes while both collapsing membrane
potential, so these labels are the thing Attene-Ramos cannot give on its own: they are how to read an MMP
result, and the 3 vs 19 split is the first contrast to look at. Rotenone is the textbook positive control.

n = 115 is still far too small for drug ranking; it is enough for the label-level enrichment test (one value
per drug), which is what section 6 proposed for Attene-Ramos and can now run on both sets with the same code.

### 9.4 Attene-Ramos re-mapped: the counts survive, not one identity does

Re-run of `map_attene_lincs_perts.sh` with both fixes:

| | before (buggy) | after |
|---|---|---|
| pairs | 1,509 | **981** (718 exact, 263 stereo-insensitive, 88 ambiguous) |
| LINCS HepG2 compounds | 904 | **890** |
| Attene-Ramos CIDs | 1,379 | 918 |

**Pairs in common between the two outputs: 0.** The headline of section 3 is therefore nearly unchanged as a
*number* — 890 mapped HepG2 compounds, 18% of the 4,930 mappable — while every compound-to-compound link in
the old file was wrong. The pair count fell because the shifted row indices were also inflating it through the
cartesian merge, not because matches were lost.

Nothing downstream had been built on the old file yet, so the only cost is the re-run. The aggregation
decisions of section 3 stand as stated: 918 CIDs collapse onto 890 compounds, so Attene-Ramos measurements
still need a per-compound rule, and the 88 ambiguous-stereo pairs should be tested as a sensitivity check
rather than merged silently.

### 9.5 What to do with this

1. **Attene-Ramos viability columns** (section 7 item 1) — unchanged and still first; it decides whether the
   test measures mitochondrial toxicity or plain cytotoxicity.
2. **Mapped x active HepG2 perturbations** — now computable for both sets once the HepG2 MetChange labels are
   down: Attene-Ramos 890, Hallinger 115.
3. **Run the enrichment test on both**, Hallinger second but read first: its mechanism labels (19 ETC
   inhibitors vs 3 uncouplers vs 4 ATP-synthase) say which direction to expect, and Attene-Ramos supplies the
   power. If MetChange separates ETC inhibitors from inactives in the 115 but not in the 890, the difference
   is the phenotype definition, not the labels.
4. **Re-check the Schuhknecht mapping.** `scripts/schuhknecht_map/` is separate code and the Schuhknecht
   inputs predate the `map_smiles.R` rewrite, so this bug is probably not in that path — but the 750 matched
   drugs behind the null result of report 05, 8.2 were mapped by *something*, and after two identity bugs in
   one day that deserves a direct check rather than an assumption.

## 10. The Schuhknecht mapping checked: sound, but a 20% undercount (added 2026-11-19, later)

Section 9.5 item 4 said the Schuhknecht mapping deserved a direct check rather than an assumption. Done.
The original is `scripts/schuhknecht_map` (rcdk `"Canonical"` + `intersect` / `match` on the canonical SMILES
of LINCS `compoundinfo_beta.txt` and sheet 2 of the supplement).

**Bug 1 is not in that code path.** It matches on *values*, never on row positions, and assembles its output
with `match(overlapping_smiles, ...)`, so every row keeps its own structure. The identity scrambling that hit
Attene-Ramos cannot happen there, and **the null of report 05, 8.2 is not a mapping artefact.**

**What it does lose**, all in the permissive-or-missing direction:

1. **Explicit hydrogens — and here they are not hypothetical: 228 of 1,520 Schuhknecht SMILES (15%) carry
   `[H]`, against 0 of 33,515 on the LINCS side.** Exactly the asymmetry of bug 2 (9.2), and it breaks the
   *flat* comparison as well as the stereo one, so dropping `"Isomeric"` was no protection.
2. **No salt stripping** (175 LINCS and 31 Schuhknecht SMILES carry a `.`).
3. **`match()` keeps the first hit per side**, so a drug present under several BRD ids contributes one.
4. **`intersect` lets `NA` match `NA`**, and `match(NA, ...)` then resolves to each side's first unparseable
   row. **The predicted single fabricated pair exists: `BRD-U01690642` (unparseable SMILES) is paired with
   Schuhknecht's `Gramicidin`** — 1 bad row in 1,040.
5. **Stereochemistry dropped silently** throughout (`"Canonical"` with no `"Isomeric"`), so every match was
   flat and racemate-vs-enantiomer pairs were never flagged.

**Re-mapped** with the fixed `map_smiles.R` (`scripts/benchmarking/map_schuhknecht_lincs_perts.sh`, new;
output `data/validation/schuhknecht/mapping_remap/`):

| | original | re-map |
|---|---|---|
| pairs | 1,040 | **1,537** (906 exact, 631 stereo-insensitive, 350 ambiguous) |
| Schuhknecht drugs | 1,040 | **1,251** |
| LINCS pert_ids | 1,040 | 1,520 |

**998 of the 1,040 original pairs reproduce exactly**, 1,039 of the 1,040 drugs are retained, and the single
drug lost is **Gramicidin** — the NA artefact, i.e. every genuine link survived. The 42 differing pairs are
drugs kept but re-assigned or extended to further `pert_id`s, which is item 3.

**212 drugs gained (+20%)**: 89 by an exact, stereo-preserving match (these were pure losses, mostly the
explicit-`[H]` and salt cases) and 123 only stereo-insensitively, 90 of them flagged ambiguous. The gained
names are visibly the salts and single enantiomers the old script could not reach — `(R) Propranolol
hydrochloride`, `(-) Isoproterenol hydrochloride`, `(S,-) Cycloserine`, `Adamantamine fumarate`.

**So: nothing to redo, something to gain.** The Schuhknecht conclusions stand as published in report 05
(1 spurious drug in 1,040 cannot move a drug-specific ρ of −0.013), and the benchmark is deprioritised anyway
(report 05, 9.3). But if it is ever revisited, the re-map offers 1,251 drugs instead of 1,040, with the
stereo-ambiguous 90 separable. By contrast Attene-Ramos (9.4) **must** be redone, because there bug 1 applied.

**Third bug in `map_smiles.R`, found by this re-run** (fixed): `largest_fragment` guarded only
`parse.smiles`, while an empty fragment (`"CC."`, a stray separator) parses to a non-NULL molecule with no
atoms and `get.atoms()` then throws a Java `IndexOutOfBoundsException` that **aborts the entire run** instead
of dropping the compound. The Hallinger and Attene-Ramos inputs happen not to contain one; Schuhknecht does.
The whole body is now guarded, with warnings treated as a parse failure.

## 11. First results: MetChange vs measured mitochondrial phenotypes (added 2026-11-19, later)

The label-level test of section 6, run on the HepG2 MetChange labels
(`results/labels/HEPG2/metchange/rank`, downloaded 2026-11-19). **This is the first external
association for the MetChange labels that is not null** — with one confound that has to travel with it.

### 11.1 Code

`scripts/mc_phenotype_assoc/`:

- `build_metabolite_sets.py` — 40 sets from Human-GEM subsystems (oxphos, TCA, glycolysis, pyruvate, PPP)
  plus curated endpoint sets (`etc_endpoints`, `energy_currency`, `adenylates`, `redox_couples`,
  `*_endpoints`), each in compartment variants (`all` / `mito` / `cyto`) and with / without hub metabolites.
- `build_phenotype_tables.py` — drug → phenotype tables keyed on the LINCS pertname, which is what the
  label columns are named after. Attene-Ramos groups: `MMP_down` (147), `inactive` (480), `MMP_up` (75),
  `cytotox` (24), `incon` (123, excluded). Hallinger: `active` (29) / `inactive` (86), with mechanism.
- `mc_phenotype_assoc.py` — one streaming pass over `log2FC.csv`; two per-drug statistics, `z` (the set's
  mean log2FC standardised by that perturbation's own mean and SD) and `auroc` (the set's mean background
  ECDF, rank-based); averaged over the drug's perturbations; then Mann-Whitney, bootstrap CI, label
  permutation, **and a size-matched random-set control** (200 random sets per set size, same pipeline), BH
  across sets. Both statistics reduce to a `nanmean` over the set's rows, which is what makes the random-set
  control affordable. Dumps `set_member_values.npz` so GSEA / NES can be added later without re-reading
  1.5 GB.
- `stratified_check.py` — the confound check of 11.4, from `per_drug_scores.tsv` alone.

**No expected sign is imposed** (the compensation hypothesis, report 05, 9.2) and **the model's activity
filter is deliberately not applied**: only 2,075 of 15,032 HepG2 perturbations have a significant MetChange
feature, and requiring one would reduce 147 vs 480 drugs to 45 vs 44. This test reads continuous log2FC, not
a selected feature list, so it does not need the filter — and responsiveness is then handled explicitly as a
confounder (11.4) rather than by throwing away 85% of the data.

### 11.2 Attene-Ramos (membrane potential): two signals

147 MMP-down vs 480 inactive drugs, 2,186 perturbations, 726 drugs scored. AUROC is of the per-drug statistic,
case vs control, so 0.5 is no difference and **below 0.5 means the set is lower in the phenotype-positive
drugs**; q is BH across the 40 sets; `rand` is the size-matched random-set p.

| contrast | set | AUROC [95% CI] | q | rand |
|---|---|---|---|---|
| MMP down | `tca_endpoints cyto_nohub` | **0.363** [0.314, 0.411] | 7e-06 | 0.005 |
| MMP down | `tca mito_nohub` (19 members) | **0.364** [0.315, 0.410] | 7e-06 | 0.010 |
| MMP down | `pyruvate mito_nohub` | 0.372 [0.325, 0.420] | 2e-05 | 0.015 |
| MMP down | `glycolysis cyto_nohub` | 0.428 [0.377, 0.476] | 0.019 | 0.119 |
| MMP down | `redox_couples` (NAD(P)+/NAD(P)H) | **0.605** [0.557, 0.655] | 4e-04 | 0.060 |
| MMP down | `oxphos mito` | 0.491 [0.440, 0.544] | 0.75 | 0.88 |
| MMP down | `etc_endpoints mito` | 0.476 [0.425, 0.526] | 0.40 | 0.59 |
| **MMP up** | `energy_currency mito` | **0.663** [0.586, 0.733] | 2e-04 | 0.030 |
| **MMP up** | `adenylates` | 0.643 [0.565, 0.718] | 6e-04 | 0.075 |
| cytotoxic (24) | best of 40 | 0.417 | 0.75 | 0.045 |

Read as a gradient, not a single hit: **TCA lowest (0.36), pyruvate intermediate, glycolysis weak
(0.43, and it does not beat random sets), PPP and OXPHOS / ETC flat.** That asymmetry is the specificity
contrast of section 9.5 behaving as designed — the mitochondrial sets move and the glycolytic ones barely do.
The effect is strongest in the **non-hub** variants, so it is not merely the currency metabolites, and
`tca mito_nohub` (19 real TCA metabolites, no hubs) carries it.

**Direction.** MMP-collapsing drugs have *lower* MetChange TCA production capacity and *higher*
NAD(P)+/NAD(P)H capacity. The TCA direction is the naive expectation (ETC block → less TCA throughput), i.e.
this is not the compensation pattern. **Hyperpolarising drugs go the other way**, with higher energy-currency
capacity — and the two phenotype directions giving opposite-signed label differences is itself a coherence
check that a single artefact would not produce.

**The cytotoxic group is flat** (24 drugs, nothing at q < 0.74), so the association is not a cytotoxicity
readout. Underpowered, but it is the right sign of evidence.

**Caveat on the on-target sets.** `oxphos` and `etc_endpoints` — the sets most directly about membrane
potential — are the ones that show *nothing*. And the OXPHOS subsystem set is 24 of 26 hub metabolites
(ATP, H+, H2O, Pi, NAD+/NADH, O2) plus ubiquinol / ubiquinone, so `oxphos_nohub` does not exist above the
size floor: for respiration the testable features essentially *are* the currency metabolites. A positive
result there could not have been shown to be pathway-specific anyway.

### 11.3 Hallinger (OCR): null at this n

29 active vs 86 inactive, and 19 ETC inhibitors vs 86 inactive. Nothing survives BH; best raw p = 0.13
(active vs inactive, `glycolysis cyto_nohub`) and 0.24 (ETCi vs inactive, `pyruvate mito`). The
mitochondrial sets do lean consistently one way in the ETCi contrast (AUROC 0.57–0.59 for `pyruvate mito`,
`etc_endpoints mito`, `oxphos all`), which is the direction Attene-Ramos shows, but every CI crosses 0.5.
**At n = 19–29 this is uninformative, not contradictory**: it cannot resolve an effect of the size seen in
Attene-Ramos (|AUROC − 0.5| ≈ 0.1). Hallinger's value stays what section 9.3 said — mechanism labels for
interpretation, not a test of its own.

### 11.4 The confound, and what survives it

**MMP-active compounds are bioactive compounds**: 31% of MMP-down drugs have a responsive MetChange
signature (a perturbation with ≥ 1 significant feature) against 9% of inactives, a 3.4-fold difference. Both
statistics are normalised within each perturbation, so a larger overall response does not mechanically raise
a set's score — but if responsive signatures have a different *shape*, responsiveness alone would reproduce
the association. `stratified_check.py` tests this directly by repeating the contrast inside each
responsiveness stratum, and by measuring the confounder on its own: responsive vs non-responsive drugs
**within the control group**.

| set | pooled | responsive (45/44) | non-responsive (102/436) | confounder: resp vs non, controls only (44/436) |
|---|---|---|---|---|
| `tca_endpoints cyto_nohub` | 0.363 (4.5e-07) | 0.383 (0.059) | **0.383 (2.3e-04)** | **0.355 (0.0015)** |
| `tca mito_nohub` | 0.364 (5.9e-07) | 0.409 (0.14) | **0.370 (4.2e-05)** | 0.389 (0.015) |
| `pyruvate mito_nohub` | 0.372 (2.8e-06) | 0.395 (0.090) | 0.392 (6.5e-04) | 0.355 (0.0015) |
| `glycolysis cyto_nohub` | 0.428 (0.008) | 0.477 (0.71) | 0.416 (0.008) | 0.455 (0.33) |
| `redox_couples` | 0.605 (1.2e-04) | 0.609 (0.078) | 0.588 (0.006) | 0.578 (0.090) |
| `oxphos mito` | 0.491 (0.75) | 0.515 (0.81) | 0.490 (0.76) | 0.468 (0.48) |
| **MMP up:** `energy_currency mito` | 0.663 (5.9e-06) | 0.565 (0.45) | **0.698 (8.4e-07)** | 0.439 (0.18) |
| **MMP up:** `adenylates` | 0.643 (6.5e-05) | 0.555 (0.52) | 0.670 (2.2e-05) | 0.460 (0.38) |

**What this says, honestly:**

1. **The TCA association is real but partly confounded.** It survives in both strata at a consistent
   0.37–0.41, and holds at p = 2e-4 among the 538 non-responsive drugs. But the confounder is of comparable
   size and the same sign: *within the inactive drugs alone*, responsive ones already show the low-TCA
   pattern (AUROC 0.355, p = 0.0015). So part of the pooled 0.363 is "this drug does something in HepG2"
   rather than "this drug depolarises mitochondria". The honest effect size is the stratified one, ~0.38, not
   the pooled one.
2. **The MMP-up / energy-currency association is clean.** Pooled 0.663, and 0.698 among non-responsive
   drugs, while the responsiveness confounder is absent and if anything points the other way (0.439, ns).
   This is the stronger of the two results, which is the opposite of what one would expect from an artefact,
   since `MMP_up` is the smaller group (75).
3. **`oxphos` is null in every stratum**, so the nulls are not a power artefact of the stratification.
4. The responsive stratum is underpowered throughout (45 vs 44), so its p-values are not the evidence; the
   agreement of its *effect sizes* with the non-responsive stratum is.

### 11.5 Next, in order

1. **The expression baseline, which is now the decisive control.** The HepG2 expr labels are on the HPC. The
   same script runs unchanged on them (`--labels <expr dir>`) given gene sets instead of metabolite sets —
   the same subsystems' GPR genes, which `build_metabolite_sets.py` can emit from the `gene_reaction_rule`
   field it already parses (note the GEM rules are Ensembl IDs while the prediction pipeline keeps symbols,
   so the ID space of the expr matrix decides whether a mapping step is needed). **If the transcriptome
   separates these groups as well as MetChange does, the stoichiometric reading has added nothing** — this
   is the same question the kNN baseline asked of the model (report 06, 6.1), and it should be answered
   before the result is presented anywhere.
2. **A responsiveness-matched design** rather than stratification: match each MMP-down drug to inactive
   drugs on `n_sig` (or regress the statistic on `n_sig` and test the residual). Cheap, and it turns 11.4's
   two strata into one estimate.
3. **A second cell line.** A549 MetChange labels exist and Attene-Ramos is HepG2 only, so this is not a
   replication of the phenotype — but the *confounder* (responsiveness shaping set scores) can be measured
   in any cell line, and should be.
4. `APR_HepG2_MitoMembPot_24h` / `_72h` (1,108 chemicals, 197 / 143 hits) with `APR_HepG2_MitoMass`
   alongside it, in `RawMC5_ToxCast_by_aeid.csv` of the Hallinger deposit — an orthogonal MMP readout in the
   same cell line, with a mass control separating "fewer mitochondria" from "depolarised mitochondria".
   Mappable through `casn`.
5. GSEA / NES per drug from `set_member_values.npz`, as the secondary statistic. Not expected to change the
   conclusion: `z` and `auroc` already agree on every row above.

**What may not be claimed yet:** that MetChange "predicts" mitochondrial toxicity. What the data support is
that the HepG2 MetChange labels carry a TCA-related signal that differs between measured MMP-active and
inactive drugs, at |AUROC − 0.5| ≈ 0.12 after the responsiveness confound is stratified out, with a cleaner
and larger signal for the hyperpolarising group — and that the pathway-specific direction (TCA, not
OXPHOS / ETC) is not the one a membrane-potential assay most obviously predicts.

## 12. Three modalities, and the expression baseline started (added 2026-11-19, later)

`mc_phenotype_assoc.py` now takes `--mode value | sign | sig` (repeatable). The refactor was verified against
section 11: **all 80 statistics shared with the pre-refactor run reproduce identical AUROCs.**

| modality | statistics | question |
|---|---|---|
| `value` | `z`, `auroc` | does the pathway move *at all* (magnitude) |
| `sign` | `sign_mean`, `sign_dev` | does it move UP or DOWN (magnitude discarded) |
| `sig` | `sig_dev`, `sig_signed_dev` | are the pathway's features *significant* more often than the drug's others |

`sign_dev` subtracts the perturbation's own mean sign, so a drug that pushes everything one way does not score
on every set; `sign_mean` is the absolute direction.

### 12.1 The sign modality confirms the TCA result, and is the *cleanest* evidence for it

Attene-Ramos, MMP-down vs inactive, `results/mc_phenotype_assoc/HEPG2/attene_allmodes/`:

| statistic | best TCA set | AUROC | q | rand |
|---|---|---|---|---|
| `z` | `tca_endpoints cyto_nohub` | 0.359 | 5e-06 | 0.005 |
| `auroc` | `tca_endpoints cyto_nohub` | 0.363 | 7e-06 | 0.005 |
| `sign_dev` | `tca cyto` | 0.415 | 0.041 | 0.020 |
| `sign_mean` | — | nothing at q < 0.05 | 0.87 | — |

All three informative statistics agree: **the TCA metabolites are relatively DOWN in MMP-collapsing drugs.**
`sign_mean` being null while `sign_dev` is significant says the drugs' overall direction varies, and TCA is
down *relative to each drug's own baseline* — which is the comparison that means something.

**And the sign modality is less confounded than the magnitude one.** Repeating 11.4's check:

| statistic / set | pooled | responsive (45/44) | non-responsive (102/436) | confounder (controls only) |
|---|---|---|---|---|
| `sign_dev` / `tca cyto` | 0.415 (0.0019) | 0.409 (0.14) | **0.432 (0.032)** | **0.453 (0.30) — absent** |
| `sign_dev` / `tca_endpoints cyto_nohub` | 0.425 (0.0062) | 0.410 (0.14) | 0.439 (0.056) | 0.513 (0.78) — absent |
| `z` / `tca_endpoints cyto_nohub` (11.4) | 0.359 | 0.383 | 0.383 (2e-04) | 0.355 (0.0015) — **present** |

So the responsiveness confound that muddies the magnitude result **does not apply to the directional one**,
which makes sense: `sign_dev` is referenced to the drug's own mean sign, and responsiveness affects magnitudes
far more than relative direction. **The directional reading is now the result to lead with**, even though its
effect size is smaller (0.415 vs 0.359).

### 12.2 The `sig` modality is a responsiveness artefact — do not use it as run here

It produced the largest apparent effects of the whole analysis (`tca mito_nohub` AUROC 0.397, q = 4e-10) and
they are **entirely spurious**:

| | pooled | non-responsive | confounder (controls only) |
|---|---|---|---|
| `sig_dev` / `tca mito_nohub` | 0.397 (4e-10) | **0.500 (p = 1)** | **0.091 (p = 2e-71)** |
| `sig_dev` / `oxphos mito` | 0.442 (5e-04) | 0.500 (p = 1) | 0.159 (p = 3e-50) |

With no significant feature, the set fraction and the column fraction are both 0, so `sig_dev` is identically
zero — hence AUROC exactly 0.5 and p = 1 in the non-responsive stratum, all ties. Since MMP-down drugs are
3.4x more likely to be responsive, the pooled test is detecting **"has any significant feature"**, not pathway
enrichment. Re-run on responsive perturbations only (`--min-nsig 1`, `attene_sig_responsive/`, 45 vs 44):
**nothing survives, best p = 0.079 and every `rand` ≥ 0.09.**

Two consequences:
1. `--mode sig` must always be run with `--min-nsig 1`; the script now says so and the launchers do it.
2. **The apparent cytotoxic-group signal is the same artefact.** In section 11 the 24 cytotoxic drugs were
   flat (value modality). In `sig` they look strong (`sig_signed_dev tca_endpoints mito` AUROC 0.660,
   q = 3e-05) — but cytotoxic drugs are 42% responsive against 9% of inactives. **The clean statement from
   11.2 stands only for the value and sign modalities.**

### 12.3 Expression baseline: started, and not yet decisive

`results/labels/HEPG2/expr/` holds `qvals.csv` and `basal.csv` but **no `log2FC.csv` and no
`fit_info_nsig.csv`**, so only the `sig` modality can run — the one modality 12.2 just showed to be
uninformative. The script now copes with the missing `fit_info_nsig.csv` (drug / dose / time parsed from the
column names, `n_sig` counted from the q-value matrix).

Gene sets: `build_metabolite_sets.py --emit genes` writes the subsystems' GPR genes, 5 sets of 23–106 genes.
The GEM states its rules in Ensembl IDs and the expression matrices are indexed by Ensembl IDs, so **no
identifier mapping is involved**. Worth noting the asymmetry: `oxphos` has **106 genes and no hubs**, so the
expression arm can test OXPHOS properly, which the metabolite side cannot (24 of its 26 members are hubs).

Result (`attene_expr_sig_responsive/`, 98 vs 265 responsive drugs): best is `glycolysis gene` AUROC 0.579,
p = 0.021, **but `rand` = 0.254 — it does not beat size-matched random gene sets**, and nothing survives BH.
On the pooled version the confounder is overwhelming as expected (0.087–0.23, p to 6e-60).

**So the baseline is not yet answered.** The suggestive comparison is that MetChange's TCA sets beat their
size-matched random controls (`rand` 0.005–0.02) while no gene set does (`rand` ≥ 0.25) — but that is
MetChange's *value/sign* result against expression's *sig* result, which is not like for like. **What is
needed is `log2FC.csv` for the HepG2 expression labels**, then `launch_assoc_expr_HEPG2.sh` runs value and
sign automatically and the comparison is set-for-set.

### 12.4 Launchers (reproducible on the HPC)

`scripts/mc_phenotype_assoc/`:

- `mc_phenotype_assoc_common.sh` — sourced config. Detects the cluster and reads labels from
  `/scratch/users/gsantamaria/test_large_files/metabognn_dev/results/labels/<cell>/...`, falling back to the
  repo path on the laptop; the structure below the base is identical, so only the base differs. Every
  variable is overridable by export (`cell`, `labels_base`, `n_perm`, ...).
- `launch_assoc_metchange_HEPG2.sh` — Attene-Ramos (3 contrasts) and Hallinger (activity and by mechanism)
  in `value`+`sign`, then `sig` with `--min-nsig 1`, then the responsive-subset sensitivity run, each
  followed by the stratified check over every statistic produced.
- `launch_assoc_expr_HEPG2.sh` — the same against the expression labels with the gene sets; skips value/sign
  automatically if `log2FC.csv` is absent and says so.
- `launch_assoc_HEPG2_slurm.sh` — `sbatch` wrapper, a 2-task array (0 = MetChange, 1 = expression),
  4 cpus / 32 GB / 6 h.

Run with `bash scripts/mc_phenotype_assoc/launch_assoc_metchange_HEPG2.sh` or
`sbatch scripts/mc_phenotype_assoc/launch_assoc_HEPG2_slurm.sh`. One modality only:
`modes="--mode sign" bash .../launch_assoc_metchange_HEPG2.sh`.

Note on cost: six statistics triple the group-test work (3 contrasts x 40 sets x 6 statistics, each with
bootstrap + permutation + 200 random sets). The laptop runs above used `--n-perm 1000 --n-boot 500`; the
launchers default to 5000 / 2000, which is what the cluster should use. The Mann-Whitney and random-set
p-values do not depend on those counts.

### 12.5 Where this leaves the result

**Stands, and is now better supported:** MMP-collapsing drugs have relatively lower MetChange TCA scores.
Three statistics agree, it beats size-matched random sets, and in the *directional* reading it is free of the
responsiveness confound (pooled 0.415, non-responsive 0.432, confounder absent).

**Withdrawn from section 11's reading:** nothing, but two qualifications. The `sig`-based numbers are not
evidence, and the clean cytotoxic control holds only for value and sign.

**Still open:** the expression baseline needs `log2FC.csv`; the effect is modest (|AUROC − 0.5| ≈ 0.09–0.14);
one cell line, one assay; `oxphos` / `etc_endpoints` remain null in every modality, so the signal is TCA-shaped
rather than respiration-shaped.

## 13. The expression baseline, run (added 2026-11-19, later still)

`log2FC.csv` for `results/labels/HEPG2/expr/` arrived, so the baseline of 11.5 / 12.3 ran properly:
`value` + `sign` on the GPR gene sets, same phenotypes, same contrasts, same pipeline
(`attene_expr_valuesign/`, 20,186 genes, 2,186 perturbations, 726 drugs; `hallinger_expr_valuesign/`).
A code gap was fixed on the way: a value/sign-only run never read the q-values, so `n_sig` was absent and the
responsiveness stratification could not run. The q-matrix is now read for `n_sig` whenever there is no
`fit_info_nsig.csv` (`--no-nsig` opts out). The rerun reproduced the first run's numbers exactly.

### 13.1 Like for like

| contrast / statistic | MetChange (metabolite sets) | expression (GPR gene sets) |
|---|---|---|
| MMP-down, `z` | `tca_endpoints cyto_nohub` **0.359**, q = 5e-06, **rand = 0.005** | `oxphos gene` **0.410**, q = 0.0035, rand = 0.154 |
| MMP-down, `auroc` | `tca_endpoints cyto_nohub` 0.363, q = 7e-06, rand = 0.005 | `oxphos gene` 0.414, q = 0.0074, rand = 0.149 |
| MMP-down, `sign_dev` | `tca cyto` **0.415**, q = 0.041, **rand = 0.020** | `oxphos gene` **0.421**, q = 0.018, rand = 0.129 |
| MMP-up, `z` | `energy_currency mito` **0.660**, q = 1.9e-04, rand = 0.025 | `pyruvate gene` **0.615**, q = 0.0056, rand = 0.164 |
| cytotoxic | flat (value / sign) | flat (best q = 0.074) |

Hallinger's expression arm is null as well (best `glycolysis gene` z, AUROC 0.549, p = 0.43), matching the
MetChange side: n = 29 resolves nothing either way.

**The headline, stated plainly: the transcriptome separates these drugs about as well as MetChange does.**
0.410 vs 0.359 for MMP-down, 0.421 vs 0.415 directionally, 0.615 vs 0.660 for MMP-up. **On raw separation the
stoichiometric reading has no advantage over reading the genes** — the same verdict the kNN baseline returned
for the model (report 06, 6.1), now for the labels.

### 13.2 Two criteria on which MetChange does look more specific

**(a) The size-matched random-set control.** MetChange's TCA sets beat random metabolite sets of their own
size (`rand` 0.005–0.020); the gene sets mostly do not (`rand` 0.129–0.164). The natural reading is that in
expression there is a *global* group difference — MMP-active drugs have broadly different transcriptomes, so
an arbitrary 100-gene set separates them too — while in MetChange the TCA sets do something their size-matched
peers do not.

**(b) The responsiveness confound.** Repeating 11.4 / 12.1 on the gene sets:

| statistic / set | pooled | responsive | non-responsive | confounder (controls only) |
|---|---|---|---|---|
| expr `z` / `oxphos gene` | 0.410 (9e-04) | 0.439 (0.076) | 0.395 (0.022) | **0.339 (p = 1.3e-09)** |
| expr `z` / `tca gene` | 0.413 (0.0014) | 0.430 (0.041) | 0.402 (0.032) | **0.363 (p = 2.3e-07)** |
| expr `sign_dev` / `oxphos gene` | 0.421 (0.0037) | 0.459 (0.24) | 0.395 (0.023) | **0.335 (p = 5.3e-10)** |
| MetChange `sign_dev` / `tca cyto` (12.1) | 0.415 (0.0019) | 0.409 (0.14) | 0.432 (0.032) | 0.453 (p = 0.30) |

The expression association is **strongly confounded by responsiveness in both modalities**, where the
MetChange directional one is not. Both survive in the non-responsive stratum (expr `oxphos` 0.395, p = 0.022),
so neither is pure confound.

### 13.3 Why I would hold (a) and (b) loosely

Both criteria compare two different feature spaces, and the differences between those spaces could produce
them without any biology:

1. **Different universes and set sizes.** 20,186 genes with sets of 23–106, against 5,707 metabolites with
   sets of 5–73. A larger random set has a more stable mean, so it is *better* at capturing a global group
   shift — which inflates expression's `rand` p-values for reasons of arithmetic, not specificity.
2. **The confounder tests have very unequal power.** Responsive share of controls: expression 55%
   (530/960 drug-variant rows), MetChange 9% (264/2880). So the confounder comparison rests on 265 vs 215
   controls for expression and only **44 vs 436** for MetChange. At that n, MetChange's 0.453 has a CI of
   roughly [0.36, 0.55] — it excludes a confounder as large as expression's 0.335 only barely. "MetChange's
   directional result is unconfounded" is therefore **weaker than it reads**; the honest version is that its
   confounder is smaller and not resolvable at this n.
3. MetChange's responsiveness *imbalance* between groups is much larger (3.4x: 31% vs 9%) than expression's
   (1.2x: 67% vs 55%), which is the opposite of what criteria (b) alone would suggest.

### 13.4 The biologically interesting part

**The two label types do not point at the same pathway.** Expression's best set for MMP-down is
**`oxphos` genes** (0.410–0.421) — the respiratory-chain genes, exactly on-target for a membrane-potential
phenotype. MetChange's signal is in **TCA**, while its `oxphos` set is null in every modality. And the
asymmetry is structural: the MetChange `oxphos` set is 24 of 26 hub metabolites, so it *cannot* be tested
cleanly, whereas `oxphos` has 106 genes and no hubs.

So the transcriptome sees the ETC genes respond; the stoichiometric labels see TCA capacity shift. That is a
more informative statement than either AUROC, and it is the first thing to follow up: the two are consistent
with one mechanism (ETC transcripts respond, TCA throughput capacity falls) but the analysis here does not
establish the link.

### 13.5 What this does to the result

- **Section 11's association stands** and 12.1's directional version is still the cleanest form of it.
- **Section 11.5's decisive question is now answered, and the answer is not the favourable one:** the
  expression baseline matches MetChange on effect size. MetChange's advantage is confined to two specificity
  criteria that 13.3 shows are partly artefacts of comparing two feature spaces.
- **For a presentation**, the defensible claim is: *the HepG2 MetChange labels carry a TCA-capacity signal
  that distinguishes measured MMP-active from inactive drugs, at a magnitude the transcriptome also achieves
  through its OXPHOS genes.* Anything stronger needs 13.3 resolved — i.e. a design where the two label types
  are compared on matched set sizes and matched responsiveness, not on their native spaces.
- Next, in order: (i) responsiveness-matched design rather than stratification, for both arms; (ii) random
  sets matched on *universe fraction* rather than absolute size, which removes objection 13.3(1);
  (iii) the same test in a second cell line for the confounder structure; (iv) the orthogonal
  `APR_HepG2_MitoMembPot` readout.

## 14. How the sets are built, and a correction on multiple testing (added 2026-11-19, final)

### 14.1 The set count, in full

11 base sets x 3 compartments x 2 hub variants = **66 candidates**, of which **26 fall below
`--min-set-size 5`** after restriction to features present in the label matrix, leaving **40**. The drops are
structural, not incidental, and they are the same fact reported in 11.2 from the other direction:

| reason | examples (size) |
|---|---|
| a currency set has no non-hub member, by construction | `adenylates all_nohub` 0, `energy_currency {all,mito,cyto}_nohub` 0, `redox_couples all_nohub` 0 |
| `oxphos` is almost entirely hubs | `oxphos all_nohub` 2 (ubiquinol + ubiquinone only), `mito_nohub` 2, `cyto_nohub` 0 |
| wrong-compartment variants are empty or tiny | `ppp mito` 0 (PPP is cytosolic), `glycolysis_endpoints mito` 2, `glycolysis mito_nohub` 4 |
| small curated sets split by compartment | `adenylates mito` 4, `redox_couples cyto` 4 |

**Genes:** the same 5 subsystems through `gene_reaction_rule`, but no compartments (genes have none) and no
curated endpoint sets (those are defined by metabolite names), so 5 x 1 x 2 = 10 rows.

### 14.2 Three levels of aggregation; the test is per set

```
level 1   within one set, within one perturbation:  5,707 features -> 1 number
          (e.g. z = (set mean log2FC - column mean) / column SD)
level 2   within one drug: its perturbations (doses x times) -> 1 number per set
          (plain mean; doses of one drug are not independent)
level 3   the test, SEPARATELY PER SET: 147 case drugs vs 480 control drugs
```

Features are aggregated, but **per set, per perturbation**. Sets are never pooled with one another: each is
its own hypothesis, 40 tests per contrast per statistic, BH across them. That is what makes the gradient of
11.2 readable — `tca` 0.36, `pyruvate` 0.37, `glycolysis` 0.43, `oxphos` 0.49 are four independent tests, not
four parts of one score. A single aggregated "metabolic score" could not have shown it.

**The sets do overlap** (TCA and pyruvate share metabolites; the curated endpoint sets are subsets of the
subsystem sets), so the tests within a family are correlated and BH is not the ideal correction. It does not
affect the conclusions at q ~ 1e-6, and the one borderline row (`sign_dev tca cyto`, q = 0.041) is discussed
in 14.3.

### 14.3 Duplicate sets: collapsed, and a correction to what I first said about it

Identical member lists do occur: **every one of the 5 `gene_nohub` sets equals its `gene` set** (no gene of
these subsystems passes the hub cutoff), and on the metabolite side `glycolysis_endpoints all_nohub` equals
`all`. So `gene_sets.tsv` has 10 rows but **5 distinct sets**, and `metabolite_sets.tsv` 40 rows but **39**.

`read_sets` now collapses identical member lists, keeps the plainest variant name (`gene` over `gene_nohub`)
and records the dropped names in a new `aliases` column.

**Correction.** My first reading was that this had inflated the published q-values and that the expression
q-values were therefore wrong. **That was incorrect.** BH computes `q = p x m / rank`; duplicating the
*entire* family multiplies `m` and `rank` by the same factor, so every q-value is unchanged (verified
numerically at 2x and 4x duplication). Because all five gene sets were duplicated exactly once, the
duplication was uniform and **the expression q-values in section 13 are correct as published**. Recomputed
explicitly:

| run | rows | BH family | q-values |
|---|---|---|---|
| expr value/sign (Attene) | 120 -> 60 | 10 -> 5 | **unchanged** |
| expr sig, responsive | 20 -> 10 | 10 -> 5 | **unchanged** |
| MetChange allmodes | 720 -> 702 | 40 -> 39 | 4th significant figure only (1.61e-08 -> 1.57e-08) |

The MetChange side shifts a little because only one of its 40 sets was duplicated, so there the duplication
was *not* uniform; 2% on a q of 1e-08 changes nothing, and `sign_dev tca cyto` (q = 0.041) is untouched.

So the dedupe is kept for three reasons, none of them the one I first gave: honest set counts (39 and 5, not
40 and 10, which overstate independence), half the compute on the expression arm, and output without paired
identical rows. **No result in sections 11-13 changes.**

## 15. Run inventory and two clarifications (added 2026-11-19, final)

### 15.1 Every combination that was run

| phenotype | labels | `value` | `sign` | `sig` | output dir |
|---|---|---|---|---|---|
| Attene-Ramos | MetChange | yes | yes | yes | `attene_allmodes/`, `attene_sig_responsive/` |
| Attene-Ramos | expression | yes | yes | yes | `attene_expr_valuesign/`, `attene_expr_sig{,_responsive}/` |
| Hallinger | MetChange | yes | yes | yes | `hallinger_allmodes/`, `hallinger_etci/` |
| Hallinger | expression | yes | yes | yes | `hallinger_expr_valuesign/`, `hallinger_expr_sig/` |

All under `results/mc_phenotype_assoc/HEPG2/`. `hallinger_expr_sig` was the last gap and is **null**
(best `glycolysis gene` AUROC 0.582, p = 0.25, 22 vs 71 responsive drugs), so **every Hallinger arm is null in
both label types and all three modalities** — one consistent statement about an n of 20-29, not four
independent ones.

### 15.2 Expression `sig`, stated fully

| run | set | AUROC | p | q | rand |
|---|---|---|---|---|---|
| pooled (147/480) | `glycolysis gene` | 0.577 | 0.0032 | **0.016** | 0.144 |
| pooled | `oxphos gene` | 0.555 | 0.034 | 0.086 | 0.239 |
| responsive (98/265) | `glycolysis gene` | 0.579 | 0.021 | 0.106 | 0.254 |
| responsive | `oxphos gene` | 0.566 | 0.052 | 0.130 | 0.323 |

**Clarification on section 12.3.** That section said nothing survives BH, which was true of the *responsive*
run; in the **pooled** run `glycolysis gene` does survive at q = 0.016, and that was not stated explicitly.
It is still not evidence, on two independent grounds: it fails the size-matched random-set control
(rand = 0.144), and the responsiveness confounder for that same set is AUROC 0.230 at p = 1.5e-26 — the
pooled `sig` test measures responsiveness, as 12.2 established for MetChange. Restricting to responsive drugs
takes it to q = 0.106. **The conclusion of 12.3 is unchanged; the statement of it was incomplete.**

### 15.3 The two expression runs are not on identical columns

`attene_expr_sig*` covers 702 drugs / 2,074 perturbations, `attene_expr_valuesign` 726 / 2,186, because the
`sig` runs were given only `--extra-contrast MMP_up:inactive` while the value/sign run also had
`cytotox:inactive`, and column selection keeps only the contrasted groups. Not an error — the MMP-down vs
inactive comparison is identical in both — but the runs are not column-for-column comparable, and a rerun
through `launch_assoc_expr_HEPG2.sh` (which passes the same contrasts to every modality) would make them so.

## 16. The responsiveness-matched design: section 13.5 was wrong (added 2026-11-19, final)

13.5 concluded that the expression baseline "matches MetChange on effect size" and that MetChange's apparent
advantage rested on two specificity criteria that 13.3 argued were partly artefacts of comparing a
5,707-metabolite space with a 20,186-gene one. **The matched design shows that conclusion was wrong.** The
pooled effects were indeed similar in size, but they were similar for different reasons: expression's is a
bioactivity signal, MetChange's is a pathway signal, and controlling responsiveness separates them.

### 16.1 `matched_check.py`

Two ways of removing the confounder while keeping one estimate, instead of 11.4's two strata:

- **matched** — 1:1 greedy nearest-neighbour matching of control drugs to case drugs on `log1p(n_sig)`,
  without replacement, caliper 0.25. Balance achieved: MetChange 140 of 147 cases matched, mean
  `log1p(n_sig)` 2.163 (cases) vs 2.144 (matched controls) vs 0.656 (all controls), residual imbalance
  p = 0.64; expression 139 of 147, 3.761 vs 3.750 vs 2.136, p = 0.97.
- **residual** — least squares of the per-drug statistic on `log1p(n_sig)` over all drugs, then the same test
  on the residuals. Keeps every drug, at the price of assuming the confounder acts linearly.

They fail differently (matching discards data, residualisation assumes a form), so their agreement is the
evidence. Both read `per_drug_scores.tsv`; no access to the label matrix is needed. Each arm is matched on
**its own** `n_sig`, since each arm's confounder is its own.

### 16.2 Result

| label type / set / statistic | pooled | **matched** | residual |
|---|---|---|---|
| MetChange `tca mito_nohub` z | 0.364 (5.9e-07) | **0.371 (2.0e-04)**, q = 0.0039 | 0.388 (4.0e-05) |
| MetChange `tca_endpoints cyto_nohub` z | 0.359 (2.2e-07) | **0.387 (1.0e-03)** | 0.394 (1.1e-04) |
| MetChange `pyruvate mito_nohub` z | 0.363 (4.9e-07) | **0.369 (1.5e-04)**, q = 0.0039 | — |
| MetChange `tca cyto` sign_dev | 0.415 (0.0019) | **0.422 (0.025)** | 0.442 (0.035) |
| MetChange `oxphos mito` z (control) | 0.476 (0.37) | 0.504 (0.91) | 0.489 (0.68) |
| expression `tca gene` z | 0.413 (0.0014) | **0.486 (0.70)**, q = 0.87 | 0.458 (0.13) |
| expression `oxphos gene` z | 0.410 (9.5e-04) | **0.450 (0.15)**, q = 0.74 | 0.448 (0.056) |
| expression `oxphos gene` sign_dev | 0.421 (0.0037) | **0.466 (0.33)** | 0.459 (0.14) |

**The expression association is essentially all responsiveness; the MetChange TCA association is not.**
Matching barely moves MetChange (0.364 -> 0.371, p = 2e-04, q = 0.0039 over the 39 sets) while expression goes
to chance (0.413 -> 0.486, q = 0.87). Matching and residualisation agree in both statistics, and the
`oxphos mito` negative control stays null in all three designs, so the matched test has not simply lost power.

### 16.3 What this does to the universe-fraction objection of 13.3

Three things, in increasing importance:

1. **`rand_p` does not depend on set size** in either space: Spearman(size, rand_p) = +0.10 (p = 0.52) over
   the 40 metabolite sets, -0.36 (p = 0.55) over the 5 gene sets. The arithmetic mechanism I proposed in
   13.3(1) - bigger sets have more stable means, so they capture a global shift better - is **not visible in
   the data**.
2. **At matched absolute size the contrast holds:** MetChange `tca all_nohub` (42 members) AUROC 0.381,
   rand = 0.0099, against expression `tca gene` (45 members) 0.413, rand = 0.209. Same size, different
   universe, different verdict.
3. **The conclusion no longer rests on random sets at all.** The matched design uses no random sets and makes
   no cross-space comparison: each arm is tested only against its own controls on its own responsiveness
   scale. That sidesteps the objection rather than answering it, and it is now the primary evidence.

**Not run:** the direct null-shape comparison - forcing random sets at matched *universe fractions* (about 149
genes to 42 metabolites) and comparing the spread of the two random nulls. It needs a `--random-set-sizes`
option and a re-stream of both matrices; given (3) it would be confirmatory, not decisive.

### 16.4 The claim, as it now stands

**The HepG2 MetChange labels carry a TCA-capacity signal that separates measured MMP-active from inactive
drugs (AUROC 0.37, p = 2e-04, q = 0.004 after responsiveness matching), and the transcriptome read the same
way does not: its apparently similar association is a bioactivity effect that matching removes.**

Still true and still limiting: the effect is modest (|AUROC - 0.5| ~ 0.13); it is one cell line and one
phenotype assay; Hallinger cannot corroborate at n = 20-29; and MetChange's `oxphos` / `etc_endpoints` sets are
null in every modality and design, so the signal is TCA-shaped, not respiration-shaped - while expression's
*pooled* signal sat in the OXPHOS genes (13.4), which is where the mechanistic follow-up should look even
though that association did not survive matching.

Next: the matched design should be applied to the `MMP_up` / energy-currency result as well (not yet done),
and to a second cell line for the confounder structure.

### 16.5 MMP-up, matched: here the baseline survives too (so 16.4 is narrowed)

The matched design applied to the hyperpolarisation contrast (75 cases), which 16.4 listed as not yet done:

| label type / set / statistic | pooled | matched | residual |
|---|---|---|---|
| MetChange `energy_currency mito` z | 0.660 (8.0e-06) | **0.651 (0.0014)** | 0.671 (1.9e-06) |
| MetChange `adenylates all` z | 0.643 (6.5e-05) | **0.635 (0.0044)** | 0.644 (6.2e-05) |
| MetChange `oxphos mito` z | 0.588 (0.015) | 0.594 (0.046) | 0.601 (0.0051) |
| MetChange `energy_currency mito` sign_dev | 0.532 (0.38) | 0.489 (0.82) | 0.529 (0.43) |
| expression `pyruvate gene` z | 0.615 (0.0013) | **0.637 (0.0040)** | 0.630 (2.8e-04) |
| expression `pyruvate gene` sign_dev | 0.623 (6.4e-04) | **0.640 (0.0034)** | 0.640 (9.6e-05) |
| expression `tca gene` z | 0.576 (0.035) | 0.591 (0.057) | 0.603 (0.0042) |

**Both label types survive matching here**, so the MMP-up association is real for the transcriptome as well as
for MetChange. Two further asymmetries worth keeping:

- MetChange's MMP-up signal is **magnitude only** — `sign_dev` is null (0.489-0.532) while `z` is 0.65 — so
  the energy-currency metabolites shift in size but not consistently in direction relative to each drug's own
  mean sign. Its MMP-down signal appears in *both* magnitude and direction. Expression's MMP-up signal appears
  in both.
- MetChange `oxphos mito` is nominally significant here (0.594, p = 0.046) where it is null in every
  MMP-down design. It is the hub-dominated set (24/26), so this is the one place the currency metabolites
  carry something, consistent with the `energy_currency` result rather than independent of it.

**Correction to 16.4.** Its claim — "the transcriptome read the same way does not [separate the groups]" — holds
for **MMP-down only**. Stated properly:

> For **MMP-down** (147 drugs), the MetChange TCA signal survives responsiveness matching (AUROC 0.371,
> p = 2e-04, q = 0.004) and the expression association does not (0.486, p = 0.70). For **MMP-up** (75 drugs),
> both survive: MetChange through the energy-currency metabolites (0.651, p = 0.0014, magnitude only) and
> expression through the pyruvate-pathway genes (0.637-0.640, p = 0.003-0.004, magnitude and direction).

So the stoichiometric reading is distinguishable from the transcriptome on one of the two phenotypes, not on
both — which is a narrower and more defensible statement than 16.4's, and the one to use.

## 17. Null shape at matched universe fraction (added 2026-11-19, final)

16.3 left this unrun and called it confirmatory. Run now, via `--random-set-sizes` and
`--dump-random-null` (new): random sets are drawn at sizes chosen so the two spaces are compared at the same
*fraction of their universe* — 13 of 5,707 metabolites (0.228%) against 45 of 20,186 genes (0.223%), and
42 (0.736%) against 149 (0.738%). Outputs `nullshape_mc/` and `nullshape_expr/`.

### 17.1 The nulls differ in shape, and not because of set size

`|AUROC - 0.5|` over 200 random sets per size, MMP-down vs inactive:

| statistic | MetChange size / frac% | median / p95 | expression size / frac% | median / p95 |
|---|---|---|---|---|
| `z` | 13 / 0.228 | 0.029 / 0.094 | 45 / 0.223 | **0.047 / 0.121** |
| `z` | 42 / 0.736 | 0.033 / 0.090 | 149 / 0.738 | **0.042 / 0.112** |
| `sign_dev` | 13 / 0.228 | 0.027 / 0.058 | 45 / 0.223 | **0.038 / 0.097** |
| `sign_dev` | 42 / 0.736 | 0.029 / 0.072 | 149 / 0.738 | **0.036 / 0.093** |

**The expression null is 1.3-1.6x wider at matched fraction** (median) and ~1.3x wider at the 95th
percentile. Within each space the null is flat in size — the full size sweep (metabolites, 5 to 149 members,
30 sizes) has median |effect| between 0.028 and 0.041 with no trend, confirming the Spearman test of 16.3(1).

So 13.3(1) was **half right**: the mechanism I proposed (bigger sets, more stable means) is not it, but the
two spaces genuinely do have differently-shaped nulls, so a raw `rand_p` comparison across them is not
apples to apples.

### 17.2 Why the direction of that difference is the interesting part

A wider null means **random gene sets separate MMP-active from inactive drugs better than random metabolite
sets do**: the transcriptomic group difference is *diffuse*, spread across the genome, which is exactly the
global effect hypothesised in 13.2(a). It is a property of the data, not an artefact of the comparison.

Each real set against its own null:

| | observed \|AUROC - 0.5\| | own null p95 | verdict |
|---|---|---|---|
| MetChange `tca mito_nohub` (19) | **0.136** | 0.090 | clears it |
| MetChange `tca all_nohub` (42) | 0.119 | 0.090 | clears it |
| expression `tca gene` (45) | 0.087 | 0.121 | inside the null |
| expression `oxphos gene` (106) | 0.090 | 0.112 | inside the null |

**Both label types carry a group difference; they differ in character.** Expression's is diffuse and, per
section 16, largely responsiveness-driven. MetChange's is concentrated in TCA and survives responsiveness
matching. That is a more informative conclusion than either the effect sizes (13.1, which looked equal) or
the `rand_p` values alone (13.2, which looked decisive for the wrong reason).

### 17.3 What may and may not be claimed from sections 11-17

**Supported:**
- The HepG2 MetChange labels separate measured MMP-active from inactive drugs through a **TCA-specific**
  signal that survives responsiveness matching (AUROC 0.371, p = 2e-04, q = 0.004) and exceeds its own
  size-matched null (0.136 vs p95 0.090).
- Read through the same GEM subsystems, **the transcriptome does not** reproduce that for MMP-down: its
  association is diffuse (inside its own null) and disappears under matching (0.486, p = 0.70).
- For **MMP-up** both label types survive matching (16.5), so the distinction is phenotype-specific.
- The stoichiometric projection **concentrates** phenotype-relevant signal into an interpretable pathway
  coordinate where gene space has it diffuse and confounded with general bioactivity.

**Not supported, and should not be written:**
- That "stoichiometry / FBA adds information that expression does not encode." **MetChange is computed from
  the expression data** (and those signatures are themselves mostly CycleGAN-inferred from 978 landmarks), so
  it cannot contain information absent from them. The defensible claim is about *representation* -
  concentration and interpretability - not information content.
- That MetChange "properly models drug-induced metabolic change". This is a group-level shift of
  |AUROC - 0.5| ~ 0.14, not per-drug accuracy; and the one per-drug test available points the other way
  (report 05, 8.2: drug-specific rho = -0.013 against measured metabolomics, below the expression-producers
  baseline). A group difference and no per-drug agreement are compatible.
- That the result is about "energy metabolism". For MMP-down, `oxphos` and `etc_endpoints` are null in every
  modality and every design; the signal is TCA. The energy-currency sets carry the MMP-*up* association
  instead.

**Untested confounder, and the most likely alternative explanation.** MMP-active compounds may be enriched for
a mechanism class or chemotype whose transcriptional response happens to project onto TCA metabolites through
the GEM. Matching on `n_sig` equalises *how much* a drug responds, not *what kind* of compound it is. Nothing
run so far addresses this. A mechanism-class (MoA) or chemical-similarity matched analysis would, and it is
the next thing I would do before presenting the result.

### 17.4 Correction: the "no added information" objection in 17.3 was wrong

17.3 said the claim that stoichiometry adds information expression does not encode "should not be written,
because MetChange is computed from the expression data". **That objection was a category error and is
withdrawn.**

It conflated information in the *data* with information in the *model*. MetChange is not a function of the
expression signature alone; it is

```
MetChange = LP( expression signature , GPR rules , stoichiometric matrix , reaction directionality )
```

solved per metabolite with COBRA + Gurobi (`scripts/metchange/run_metchange_exe{,_lp1}.py`). The GEM is a
**prior that is not derivable from the LINCS data**. A gene-level summary of a pathway cannot represent that
in A -> B -> C a fall in B interrupts flux to C; a constraint-based model over the stoichiometric matrix can.
So "expression plus stoichiometry explains more than expression alone" is a coherent claim, and the
information it refers to lives in the metabolic network, not in the measurements.

**The comparison in sections 13-17 is in fact well posed for exactly this claim**, better than 17.3 implied:
the expression arm already uses part of the GEM — its sets are the GPR genes of each subsystem's reactions —
so both arms share the gene-to-pathway mapping and the subsystem definition, and what distinguishes them is
stoichiometry, topology, directionality and the LP. The contrast isolates the structural contribution rather
than the presence of a GEM.

**So the defensible form of the claim is:**

> On MMP-down in HepG2, combining the expression signature with the GEM's stoichiometry, GPR rules and
> network topology (MetChange) separates measured phenotype-active from inactive drugs — AUROC 0.371 after
> responsiveness matching (p = 2e-04, q = 0.004), clearing its own size-matched null (0.136 vs p95 0.090) —
> where aggregating the *same* expression values over the *same* GEM subsystems does not: that association is
> diffuse (inside its own null) and vanishes under matching (0.486, p = 0.70). The difference between the two
> is the stoichiometric structure, not additional measurement.

**The caveats that do survive**, and they are narrower than 17.3's:

1. **MMP-down only.** On MMP-up both label types survive matching (16.5), so the structural advantage is
   phenotype-specific, not general.
2. **The comparator is one gene-level aggregation** — the mean (or rank mean) of a subsystem's GPR genes.
   That is the natural matched baseline, but it is not the strongest conceivable expression-based method: a
   topology-aware gene-level approach (pathway impact analysis / SPIA-style, which also uses network
   structure) or a learned model could narrow the gap. The claim is therefore "beats subsystem-level
   aggregation of the same genes", not "beats any method that starts from expression".
3. **The MoA / chemotype confounder** of 17.3 still applies to MetChange's own signal, and is untested.
4. One cell line, one phenotype assay.

Claim 1 of the same discussion ("MetChange properly models real drug-induced metabolic change") is
**unaffected by this correction** and the reservations in 17.3 stand: a group-level shift of
|AUROC - 0.5| ~ 0.14 is not per-drug accuracy, the per-drug test available points the other way (report 05,
8.2), and for MMP-down the signal is TCA rather than the energy / ETC sets.
