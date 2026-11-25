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
