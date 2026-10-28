# GNN-based metabolic repurposing tool development

This repository holds the code for developing the prototype of a GNN-based tool for identify chemical perturbations that shift the transcriptome towards a state that is consistent with altered metabolic perturbations queried by the user. 

For the initial prototype, we will work on a single cell line: A549, given the that there is an extensive metabolomic dataset on this same cell line, in response to multiple perturbations that appear in the transcriptomic dataset we use, LINCS L1000 CycleGAN data. ([Schuhknecht et al. 2025](https://doi.org/10.1038/s41587-024-02524-5)).

The repository is split in several modules.

## 1. Parse LINCS data (./scripts/parse_lincs_data)
Extract perturbations of a given cell type from [Ma'ayan lab GCTX files](https://maayanlab.cloud/sigcom-lincs/#/Download) (Level 3).
Convert gene HGNC symbols to ENSEMBL gene IDs, as Human2-GEM model uses them in the GPRs. Generates mdat.csv and expr.csv, which are aligned by mdat.csv id column.

## 2. Perform MetChange (./scripts/metchange)
Compute ranks penalty on the ENSEMBL-named gene expression matrix, splits it in multiple matrices specified by the number of batches determined by the user, and performs metchange.

## 3. Compute metabolic labels (./scripts/compute_labels)
Perform DE on the metchange results using metChange results, accounting for batch effects. Generates four files:
- basal.csv
- log2FC.csv
- pvals.csv
- qvals.csv
This script is also used on the ENSEMBL expression data to compute basal level of the genes for each perturbation.