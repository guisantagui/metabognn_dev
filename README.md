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

### Papers with data for validating the tool.
- Untargeted stable isotope-resolved metabolomics to assess the effect of PI3Kbeta inhibiton on metabolic pathway activities in a PTEN null breast cancer cell line. (10.3389/fmolb.2022.1004602). X13CMS Isotope labelling.
- Respirometric Screening and Characterization of Mitochondrial Toxicants Within the ToxCast Phase I and II Chemical Libraries. (https://doi.org/10.1093/toxsci/kfaa059). Seahorse.
- Profiling of the Tox21 chemical collection for mitochondrial function to identify compounds that acutely decrease mitochondrial membrane potential. (https://doi.org/10.1289/ehp.1408642). Mitochondrial membrane potential change in response to 10K drugs.
- A framework for large-scale metabolome drug profiling links coenzyme A metabolism to the toxicity of anti-cancer drug dichloroacetate. (https://doi.org/10.1038/s42003-018-0111-x). Time-resolved metabolomics.
- Schuhknecht
- Sarvin
- https://metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyType=&StudyID=ST000168
- Metabolic drug survey highlights cancer cell dependencies and vulnerabilities (https://doi.org/10.1038/s41467-021-27329-x). Seahorse of basal state.
- A functional analysis of 180 cancer cell lines reveals conserved intrinsic metabolic programs (https://pubmed.ncbi.nlm.nih.gov/36321552/). Only basal state--> no drugs arm.
- https://metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyType=&StudyID=ST004447
- https://metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyType=&StudyID=ST004448
- https://metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyType=&StudyID=ST004449
- https://metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyType=&StudyID=ST004450
- https://metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyType=&StudyID=ST004457
- https://metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyType=&StudyID=ST004458
- https://metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyType=&StudyID=ST004460
- https://metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyType=&StudyID=ST003944
- https://metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyType=&StudyID=ST003945
