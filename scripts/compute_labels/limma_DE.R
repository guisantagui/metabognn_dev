library(limma)
library(data.table)
library(future)
library(future.apply)
library(argparser)

# Functions
################################################################################
write_table_fast <- function (
        df,
        f,
        row.names = T,
        col.names = T,
        sep = ",",
        sep2 = c("", "|", ""),
        quote = "auto"
) {
    if (row.names) {
        rn <- rownames(df)
        df <- data.table(df)
        df[, `:=`(V1, rn)]
        setcolorder(df, c("V1", setdiff(names(df), "V1")))
    }
    else {
        df <- data.table(df)
    }
    fwrite(df, f, col.names = col.names, sep = sep, sep2 = sep2, quote = quote)
}


# Terminal argument parser
################################################################################
parser <- arg_parser("Do limma DE analysis, blocking by batch.")

parser <- add_argument(parser = parser,
                       arg = c("--expr",
                               "--mdat",
                               "--workers",
                               "--out_dir"),
                       help = c("Expression matrix (CSV)",
                                "Metadata file (CSV)",
                                "Number of parallel workers",
                                "Output directory to save the results."),
                       flag = c(F, F, F, F))

parsed <- parse_args(parser)

# Directory stuff
################################################################################

expr_f <- parsed$expr
mdat_f <- parsed$mdat
workers <- as.numeric(parsed$workers)
out_dir <- parsed$out_dir

if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
}

# Load the data
################################################################################

message(sprintf("Loading expression data from %s", expr_f))
# One numeric matrix, filled from the data.table block by block while the
# copied columns are dropped: a data.frame with one column per sample, apply()
# and row subsetting held several full copies at once.
dt <- fread(expr_f, header = TRUE)
feature_ids <- dt[[1]]
sample_ids <- names(dt)[-1]

# Features with no variance (constant or with missing values) are dropped:
# min and max per row, accumulated over the columns.
row_min <- rep(Inf, nrow(dt))
row_max <- rep(-Inf, nrow(dt))
for (s in sample_ids) {
    row_min <- pmin(row_min, dt[[s]])
    row_max <- pmax(row_max, dt[[s]])
}
keep <- !is.na(row_min) & !is.na(row_max) & row_max > row_min
message(sprintf("%d features, %d with no variance or missing values removed",
                length(keep), sum(!keep)))

expr <- matrix(NA_real_, nrow = sum(keep), ncol = length(sample_ids),
               dimnames = list(feature_ids[keep], sample_ids))
block <- 1000
for (start in seq(1, length(sample_ids), by = block)) {
    cols <- sample_ids[start:min(start + block - 1, length(sample_ids))]
    for (s in cols) expr[, s] <- dt[[s]][keep]
    dt[, (cols) := NULL]
}
rm(dt, row_min, row_max, keep)
invisible(gc())

message(sprintf("Loading and parsing metadata from %s", mdat_f))
mdat <- as.data.frame(fread(mdat_f, header = TRUE))

# Samples are selected below by metadata position, so the metadata rows must
# follow the expression columns.
if (!identical(mdat$id, colnames(expr))) {
    if (!setequal(mdat$id, colnames(expr))) {
        stop(sprintf("Expression columns and metadata ids differ (%d only in metadata, %d only in expression)",
                     length(setdiff(mdat$id, colnames(expr))), length(setdiff(colnames(expr), mdat$id))))
    }
    message("Metadata rows are in a different order than the expression columns: reordering them.")
    mdat <- mdat[match(colnames(expr), mdat$id), ]
}
#mdat$batch <- sapply(
#        mdat$rna_plate,
#        function(x) paste(strsplit(x, split = "_")[[1]][1:3], collapse = "_")
#)

perts <- unique(mdat$perturbation[mdat$exper_group != "ctrl"])
genes_all <- rownames(expr)

# Samples of one perturbation's model: its own samples and the controls of the
# same timepoint in its batches.
fit_samples <- function(pert) {
    pert_bool <- mdat$perturbation == pert
    timepoint <- unique(mdat$timepoint[pert_bool])
    ctrl_bool <- (
        mdat$exper_group == "ctrl" &
        mdat$timepoint == timepoint &
        mdat$batch %in% unique(mdat$batch[pert_bool])
    )
    pert_bool | ctrl_bool
}

# limma for one perturbation: genes x (basal, log2FC, pval, qval), or NULL if
# there are no matched controls.
fit_one <- function(pert) {
    sel <- fit_samples(pert)
    mdat_filt <- mdat[sel, ]
    if (length(unique(mdat_filt$perturbation)) < 2) return(NULL)
    dat_filt <- expr[, sel, drop = FALSE]
    design <- model.matrix(~ perturbation, data = mdat_filt)

    if (length(unique(mdat_filt$batch)) > 1) {
        corfit <- duplicateCorrelation(dat_filt, design, block = mdat_filt$batch)
        fit <- lmFit(dat_filt, design, block = mdat_filt$batch,
                     correlation = corfit$consensus.correlation)
    } else {
        fit <- lmFit(dat_filt, design)
    }
    fit <- eBayes(fit)

    # sort.by = "none": keep the input gene order (topTable sorts by B by
    # default, which differs between perturbations)
    results <- topTable(fit, coef = 2, number = Inf, adjust.method = "BH",
                        confint = TRUE, sort.by = "none")
    out <- cbind(fit$coefficients[rownames(results), 1], results$logFC,
                 results$P.Value, results$adj.P.Val)
    if (!identical(rownames(results), genes_all)) {
        out <- out[match(genes_all, rownames(results)), , drop = FALSE]
    }
    dimnames(out) <- NULL
    # Forked workers inherit the parent's garbage-collection threshold, which is
    # large once the expression matrix is loaded: without an explicit gc() each
    # worker accumulates GBs of discarded fits before R collects them. A full
    # gc() costs ~10-100 ms, so it is run after large models and every 8 fits.
    big <- length(dat_filt) > 1e7
    rm(dat_filt, fit, results, mdat_filt, design)
    gc_counter <<- gc_counter + 1
    if (big || gc_counter %% 8 == 0) invisible(gc(verbose = FALSE))
    out
}
gc_counter <- 0   # per worker: each forked worker modifies its own copy

# Do parallelized limma analysis
################################################################################
# Workers are forked and read expr without copying it. Perturbations are run in
# rounds of workers x 4, largest models first (from 2 to thousands of samples,
# so that the slowest ones do not end up last), and each round's results go
# straight into the output matrices: the parent never holds all the per-
# perturbation results next to the final matrices.

options(future.globals.maxSize = 30 * 1024^3)
plan(multicore, workers = workers)

n_samples <- vapply(perts, function(p) sum(fit_samples(p)), integer(1))
order_run <- order(n_samples, decreasing = TRUE)
message(sprintf("%d perturbations, %d to %d samples per model (median %d)", length(perts),
                min(n_samples), max(n_samples), as.integer(median(n_samples))))

empty <- function() matrix(NA_real_, nrow = length(genes_all), ncol = length(perts),
                           dimnames = list(genes_all, perts))
basal <- empty()
log2FCs <- empty()
pvals <- empty()
qvals <- empty()

message("Running limma...")
round_size <- workers * 4
n_rounds <- ceiling(length(perts) / round_size)
for (r in seq_len(n_rounds)) {
    idx <- order_run[((r - 1) * round_size + 1):min(r * round_size, length(perts))]
    res <- future_lapply(perts[idx], fit_one, future.seed = TRUE)
    for (k in seq_along(idx)) {
        m <- res[[k]]
        if (is.null(m)) next
        j <- idx[k]
        basal[, j] <- m[, 1]
        log2FCs[, j] <- m[, 2]
        pvals[, j] <- m[, 3]
        qvals[, j] <- m[, 4]
    }
    rm(res)
    invisible(gc(verbose = FALSE))
    message(sprintf("round %d / %d done (%d perturbations)", r, n_rounds, min(r * round_size, length(perts))))
}
message("limma done.")
message(sprintf("%d perturbations without matched controls (NA columns)",
                sum(colSums(!is.na(log2FCs)) == 0)))

write_table_fast(basal, f = file.path(out_dir, "basal.csv"))
write_table_fast(log2FCs, f = file.path(out_dir, "log2FC.csv"))
write_table_fast(pvals, f = file.path(out_dir, "pvals.csv"))
write_table_fast(qvals, f = file.path(out_dir, "qvals.csv"))

message(sprintf("basal.csv saved at %s", out_dir))
message(sprintf("log2FC.csv saved at %s", out_dir))
message(sprintf("pvals.csv saved at %s", out_dir))
message(sprintf("qvals.csv saved at %s", out_dir))