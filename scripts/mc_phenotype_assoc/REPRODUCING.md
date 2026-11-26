# Reproducing the phenotype-association results

Which script produces which reported number. **Scope: `scripts/mc_phenotype_assoc/` only** — the reports
before 2026-11-19 concern model training, network construction and signature selection, which live in
`scripts/predict_targets/`, `scripts/build_network/` and `scripts/parse_lincs_data/` and are not covered here.

Written 2026-11-25 because the mapping existed only in conversation, which is how the gaps below went
unnoticed. Keep it current when a launcher is added.

## The two stages

Everything splits into a **scoring** stage that needs the GB-scale label matrices (so: the cluster) and an
**analysis** stage that reads only the ~2 MB `per_drug_scores.tsv` files (so: anywhere).

```
cluster                                            local
-------                                            -----
prepare_assoc_inputs.sh   sets, phenotype tables
submit_assoc.sh           Attene/Hallinger, 9 cells  ->  reproduce_report09.sh
launch_prism_slurm.sh     PRISM, both arms, 5 cells  ->  reproduce_report10.sh
```

## Report 08 — `2026-11-20_phenotype_confounders.md`

| section | produced by | status |
|---|---|---|
| the per-cell associations, stratified and responsiveness-matched checks | `launch_assoc_HEPG2_slurm.sh`, `launch_assoc_A549_slurm.sh`, or `submit_assoc.sh` for all cells — these call `run_assoc`, `run_strat`, `run_matched` | scripted |
| §3 declared negative controls | `launch_assoc_metchange.sh` (the `attene_negctrl` run) | scripted |
| §3.2 network-position splits, §5 potency/dose, the van Elteren MoA strata | **`confounder_checks.py` and `dose_exposure_check.py` — referenced by NO launcher** | **GAP** |

## Report 09 — `2026-11-23_cross_cell_specificity.md`

| section | produced by | status |
|---|---|---|
| §3, §3.1 cross-cell MetChange | `reproduce_report09.sh 3` / `3.1` | scripted |
| **§4 MetChange vs expression head-to-head** | **`reproduce_report09.sh 4`** | scripted |
| §5 O-glycan negative control | `reproduce_report09.sh 5` | scripted |
| §6.3 cross-cell signature conservation | **`signature_conservation.py` — no launcher**; also needs local label matrices, not just scores | **GAP** |
| §9, §9.6 Naia ATP, incl. its expression arm | `run_naia_assoc.sh` | scripted |

§4 is the most load-bearing number in the project for the central thesis (MetChange 0.365 against expression
0.500 under the same matched design) and had no launcher until 2026-11-25.

## Report 10 — `2026-11-24_prism_gcp_cellpainting.md`

| section | produced by | status |
|---|---|---|
| §5 xu AC16 MMP | **`build_xu_phenotype.py` + an `mc_phenotype_assoc.py` call — no launcher** | **GAP** (the test was underpowered, so low priority) |
| §8–§13 the PRISM differential, all variants | `reproduce_report10.sh` | scripted |
| §12 the conditional-on-expression headline | `reproduce_report10.sh 12` | scripted |
| §14 producer-reachability filter | `reproduce_report10.sh 14`, via `build_target_relevance.py` | scripted |
| §15 which cell lines are worth computing | `reproduce_report10.sh 15`, via `prism_cell_coverage.py` | scripted |

## Report 11 — `2026-11-25_nci60_exchange_flux_feasibility.md`

| section | produced by | status |
|---|---|---|
| §8 label-matrix dimensionality, both arms, HEPG2 + A549 | `feature_dimensionality.py` | scripted |
| §8.4 set specificity against size-matched random sets | `feature_dimensionality.py` (same run) | scripted |
| §1–§6 the NCI-60 feasibility checks themselves | ad-hoc, run inline on 2026-11-25 | **GAP** |
| §9 the producer-average null | `reproduce_report11.sh 9`, via `build_derived_labels.py` | scripted |
| §10 four-arm comparison, residual vs producer-average | `reproduce_report11.sh 10` | scripted |
| §10.1 the four-way table itself | `reproduce_report11.sh 11`, via `compare_arms.py` | scripted |

`feature_dimensionality.py` needs only local label matrices and runs in ~3 min for both cells. Its functions
have working defaults and are importable, so the analysis can be stepped through rather than only run:

```python
from feature_dimensionality import *
X, ids, cols = load_matrix(label_path('HEPG2', 'metchange'))
spectrum(X)                                   # PC1, participation ratio, median pairwise r
set_specificity(X, ids, read_sets()[('dntp_endpoints', 'all')])
```

**Read `_check_no_double_centring` before writing any variant of this analysis.** Centring the standardised
matrix across features removes the shared component being measured, and the first pass at §8 did exactly
that and reported 9.55 effective dimensions instead of 1.10.

## Known gaps, in priority order

1. **`confounder_checks.py` / `dose_exposure_check.py`** — report 08's van Elteren MoA strata, potency and
   dose analyses. These are substantive controls and the only record of how they were invoked is report 08's
   prose.
2. **`signature_conservation.py`** — report 09, §6.3. Needs the label matrices locally, so it cannot join
   `reproduce_report09.sh` as it stands.
3. **`build_xu_phenotype.py`** — report 10, §5. Low priority; the result was underpowered by construction.

## The structural caveat

A launcher pins a flag combination; it does not make it legible. `prism_differential.py` is ~900 lines with
33 flags and the §12 result depends on six of them. **For a load-bearing claim the better form is a short
script with the specification hardcoded and no arguments** — forty lines printing the numbers, reviewable in
five minutes — with the flexible script kept for exploration. Not yet done for §12 or report 09 §4.

## Conventions that must not drift

### Gene symbol -> Ensembl

**One mapper, one file, always with the GEM overrides.**

```python
import sys; sys.path.insert(0, 'scripts/parse_lincs_data')
from gene_id_mapping import symbol_to_ensembl_map, gem_symbol_overrides
mapping = symbol_to_ensembl_map(
    pd.read_csv('data/utilities/hgnc_complete_set.txt', sep='\t', dtype=str),
    gem_symbol_overrides(gem),            # NOT {} - see below
)
```

`symbol_to_ensembl_map` takes HGNC approved symbols first, then previous and alias symbols **only** where they
are unambiguous and are not themselves an approved symbol. `gem_symbol_overrides` then takes precedence, and
exists because HGNC and Human-GEM assign **different Ensembl ids** to four genes: **SOD2, LRTOMT, PPIAL4C,
SLC37A4**. Without the overrides those four fail to match the GEM.

Every production consumer does this correctly: `metchange/create_penalties.py`,
`metchange/prepare_input_matrices.py`, `parse_lincs_data/expr_symb_to_ensembl.py`,
`metpath/build_metpath_pathways.py`, `metpath/score_metpath.py`. An ad-hoc check written on 2026-11-25 passed
`{}` instead and understated GEM gene coverage of the NCI-60 expression data by three genes (98.41% against
98.53%). The error was small because only four genes diverge - but one of them is SOD2, so the same slip
inside a mitochondrial analysis would not be.

**If you are writing a one-off that maps symbols, import the function. Do not reimplement it, and do not pass
empty overrides.**

### Aggregating duplicate genes

Several expression sources have multiple rows per symbol (CellMiner RNA-seq: 24,911 rows, 23,808 distinct
symbols, 1,103 duplicated). MetChange is rank-based, so a monotone transform of the values is harmless - but
**summing duplicates is only valid in linear space**, which is why `prepare_input_matrices.py` converts from
LINCS's log scale before summing. For a log-scale source either convert first, or aggregate with a
rank-preserving function such as max.

### MetChange score direction

`run_metchange_exe.py` minimises the penalty sum (`ModelSense = GRB.MINIMIZE`); its comment at line 266 reads
**"a metchange score of zero means the metabolite can be easily produced"**. The label is therefore a **cost**:

* MetChange log2FC **up**  = production **harder**
* expression  log2FC **up**  = genes **up**, so capacity **easier**

**The two arms are sign-inverted.** Any head-to-head reading of a MetChange AUROC against an expression
AUROC has to flip one of them, and an apparent disagreement between the arms may be agreement. Report 11
§10.3 found exactly that: producer-average 0.259 against MetChange 0.630 on the same contrast are the same
biological statement. A set-level AUROC below 0.5 on the MetChange arm means the case drugs find that set
*easier* to produce, which is not what "the pathway is impaired" would predict.
