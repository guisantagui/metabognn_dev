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
