#!/usr/bin/env Rscript
# For each perturbation, whether limma_DE.R (before the reference level was set explicitly) used it as the
# reference level, i.e. returned control - treated: the first level of factor() over the perturbation names of its
# model (its own samples and the controls fit_samples() selects). Run it in the same environment (R build, locale)
# as the limma jobs, e.g. inside an sbatch job like launch_limma_DE_generic.sh.
# Usage: Rscript limma_reference_order.R <mdat.csv> <out.tsv>
suppressMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
mdat <- as.data.frame(fread(args[1], header = TRUE))
message(sprintf("LC_COLLATE: %s", Sys.getlocale("LC_COLLATE")))

perts <- unique(mdat$perturbation[mdat$exper_group != "ctrl"])
flipped <- vapply(perts, function(pert) {
    pert_bool <- mdat$perturbation == pert
    timepoint <- unique(mdat$timepoint[pert_bool])
    ctrl_bool <- mdat$exper_group == "ctrl" & mdat$timepoint == timepoint &
        mdat$batch %in% unique(mdat$batch[pert_bool])
    lv <- unique(mdat$perturbation[pert_bool | ctrl_bool])
    length(lv) >= 2 && levels(factor(lv))[1] == pert
}, logical(1))
fwrite(data.table(perturbation = perts, flipped_R = flipped), args[2], sep = "\t")
message(sprintf("%d of %d perturbations were the reference level (%.1f%%): saved %s",
                sum(flipped), length(flipped), 100 * mean(flipped), args[2]))
