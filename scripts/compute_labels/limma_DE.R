library(limma)
library(data.table)
library(parallel)
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

# limma::duplicateCorrelation(M, design, block = block), same result to the last
# digit, for matrices without missing values. limma fits each gene with
# statmod::mixedModel2Fit, which recomputes a full SVD of the (samples x
# samples) projected block design for every gene, although it does not depend on
# the gene: ~2 s per gene with 4,000 samples, i.e. ~12 h per perturbation. Here
# the SVD is done once and the genes are projected together, in chunks; the
# per-gene gamma GLM (glmgam.fit) is the same. Other cases go to limma.
dup_cor_fast <- function(M, design, block, trim = 0.15, chunk = 2000L) {
    design <- as.matrix(design)
    QR <- qr(design)
    if (anyNA(M) || max(table(block)) < 2L || QR$rank < ncol(design)) {
        return(duplicateCorrelation(M, design, block = block, trim = trim))
    }
    A <- factor(block)
    nobs <- ncol(M)
    design.block <- model.matrix(~A)[, -1, drop = FALSE]
    if (!(nobs > ncol(design) + 2L && nlevels(A) > 1L && nlevels(A) < nobs - 1L) ||
        max(abs(qr.qty(QR, design.block)[-(1:QR$rank), ])) < 1e-08) {
        return(duplicateCorrelation(M, design, block = block, trim = trim))
    }
    keep <- (QR$rank + 1):nobs
    mq <- length(keep)
    s <- La.svd(qr.qty(QR, model.matrix(~0 + A))[keep, , drop = FALSE], nu = mq, nv = 0)
    d <- rep(0, mq)
    d[seq_along(s$d)] <- s$d^2
    dx <- cbind(Residual = 1, Block = d)
    use_glm <- mq > 2 && sum(abs(d) > 1e-15) > 1 && var(d) > 1e-15
    rho <- rep_len(NA_real_, nrow(M))
    for (start in seq(1L, nrow(M), by = chunk)) {
        g <- start:min(start + chunk - 1L, nrow(M))
        dy <- crossprod(s$u, qr.qty(QR, t(M[g, , drop = FALSE]))[keep, , drop = FALSE])^2
        lf <- lm.fit(dx, dy)
        for (k in seq_along(g)) {
            vc <- lf$coefficients[, k]
            if (use_glm) {
                start_k <- if (all(lf$fitted.values[, k] >= 0)) vc else c(Residual = mean(dy[, k]), Block = 0)
                vc <- tryCatch(suppressWarnings(statmod::glmgam.fit(
                    dx, dy[, k], coef.start = start_k, tol = 1e-6, maxit = 20)$coefficients),
                    error = function(e) NA)
            }
            if (!is.na(vc[1])) rho[g[k]] <- vc[2] / sum(vc)
        }
    }
    rhomin <- 1 / (1 - max(table(block))) + 0.01
    if (min(rho, 0, na.rm = TRUE) < rhomin) rho[rho < rhomin] <- rhomin
    if (max(rho, 0, na.rm = TRUE) > 0.99) rho[rho > 0.99] <- 0.99
    arho <- atanh(rho)
    mrho <- tanh(mean(arho, trim = trim, na.rm = TRUE))
    list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
}

# limma for one perturbation: genes x (basal, log2FC, pval, qval), or NULL if
# there are no matched controls, plus the model's size and correlation.
fit_one <- function(pert) {
    sel <- fit_samples(pert)
    mdat_filt <- mdat[sel, ]
    info <- list(n_samples = sum(sel), n_ctrl = sum(mdat_filt$exper_group == "ctrl"),
                 n_batches = length(unique(mdat_filt$batch)), correlation = NA_real_)
    if (length(unique(mdat_filt$perturbation)) < 2) return(list(out = NULL, info = info))
    dat_filt <- expr[, sel, drop = FALSE]
    design <- model.matrix(~ perturbation, data = mdat_filt)

    if (info$n_batches > 1) {
        corfit <- dup_cor_fast(dat_filt, design, block = mdat_filt$batch)
        info$correlation <- corfit$consensus.correlation
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
    list(out = out, info = info)
}

# Checkpoints
################################################################################
# Each perturbation's result is saved to <out_dir>/checkpoints/<index>.rds as
# soon as it is done (written to a temporary file and renamed, so a job killed
# mid-write leaves no partial result). Rerunning the same command skips the
# perturbations already saved; the final CSVs are built from the checkpoints.
# The manifest makes sure the checkpoints come from the same perturbations and
# genes.

ckpt_dir <- file.path(out_dir, "checkpoints")
dir.create(ckpt_dir, showWarnings = FALSE)
ckpt_f <- function(j) file.path(ckpt_dir, sprintf("%05d.rds", j))
manifest_f <- file.path(ckpt_dir, "manifest.rds")
manifest <- list(perts = perts, genes = genes_all)
if (file.exists(manifest_f)) {
    if (!identical(readRDS(manifest_f), manifest)) {
        stop(sprintf("%s holds results for other perturbations or genes: delete it or use another --out_dir",
                     ckpt_dir))
    }
} else {
    saveRDS(manifest, manifest_f)
}
done <- file.exists(ckpt_f(seq_along(perts)))

# Do parallelized limma analysis
################################################################################
# One forked process per perturbation (mclapply without prescheduling): it reads
# expr without copying it, returns all its memory when it exits, and the next
# perturbation goes to whichever worker is free. Largest models first (from 2 to
# thousands of samples), so that the slowest ones do not end up last.

message(sprintf("BLAS: %s", extSoftVersion()[["BLAS"]]))
n_samples <- vapply(perts, function(p) sum(fit_samples(p)), integer(1))
message(sprintf("%d perturbations, %d to %d samples per model (median %d)", length(perts),
                min(n_samples), max(n_samples), as.integer(median(n_samples))))
message(sprintf("%d already in %s, %d to run with %d workers", sum(done), ckpt_dir, sum(!done), workers))

todo <- which(!done)
todo <- todo[order(n_samples[todo], decreasing = TRUE)]
t_start <- Sys.time()
run_one <- function(j) {
    t0 <- Sys.time()
    res <- tryCatch(fit_one(perts[j]), error = function(e) e)
    if (inherits(res, "error")) {
        message(sprintf("[%s] ERROR %s: %s", format(Sys.time(), "%F %T"), perts[j], conditionMessage(res)))
        return(conditionMessage(res))
    }
    res$pert <- perts[j]
    res$info$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    tmp <- paste0(ckpt_f(j), ".tmp")
    saveRDS(res, tmp, compress = FALSE)
    file.rename(tmp, ckpt_f(j))
    n_done <- length(list.files(ckpt_dir, pattern = "^[0-9]+\\.rds$"))
    elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
    message(sprintf("[%s] %d / %d done (%dh%02dm elapsed) | %s: %d samples, %d batches, %.1f s",
                    format(Sys.time(), "%F %T"), n_done, length(perts), elapsed %/% 60, round(elapsed %% 60),
                    perts[j], res$info$n_samples, res$info$n_batches, res$info$seconds))
    NULL
}

message("Running limma...")
errs <- mclapply(todo, run_one, mc.cores = workers, mc.preschedule = FALSE)
failed <- todo[!vapply(errs, is.null, logical(1))]
if (length(failed)) {
    stop(sprintf("%d perturbations failed (%s%s); rerun to retry them", length(failed),
                 paste(head(perts[failed], 5), collapse = ", "), if (length(failed) > 5) ", ..." else ""))
}
message("limma done.")

# Collect the checkpoints and save
################################################################################

empty <- function() matrix(NA_real_, nrow = length(genes_all), ncol = length(perts),
                           dimnames = list(genes_all, perts))
basal <- empty()
log2FCs <- empty()
pvals <- empty()
qvals <- empty()
info <- vector("list", length(perts))
for (j in seq_along(perts)) {
    res <- readRDS(ckpt_f(j))
    stopifnot(identical(res$pert, perts[j]))
    info[[j]] <- as.data.frame(c(list(perturbation = perts[j]), res$info))
    if (is.null(res$out)) next
    basal[, j] <- res$out[, 1]
    log2FCs[, j] <- res$out[, 2]
    pvals[, j] <- res$out[, 3]
    qvals[, j] <- res$out[, 4]
}
message(sprintf("%d perturbations without matched controls (NA columns)",
                sum(colSums(!is.na(log2FCs)) == 0)))

write_table_fast(basal, f = file.path(out_dir, "basal.csv"))
write_table_fast(log2FCs, f = file.path(out_dir, "log2FC.csv"))
write_table_fast(pvals, f = file.path(out_dir, "pvals.csv"))
write_table_fast(qvals, f = file.path(out_dir, "qvals.csv"))
fwrite(rbindlist(info), file.path(out_dir, "fit_info.csv"))

message(sprintf("basal.csv, log2FC.csv, pvals.csv, qvals.csv and fit_info.csv saved at %s", out_dir))
message(sprintf("The checkpoints in %s can be deleted.", ckpt_dir))
