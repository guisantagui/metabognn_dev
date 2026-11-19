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
