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

parser <- add_argument(parser = parser, arg = "--reader", default = "fread",
                       help = paste("How the matrix is read: fread (fast) or lines (line by line,",
                                    "for files that crash fread, e.g. ~150k columns of MetChange)."))

parser <- add_argument(parser = parser, arg = "--perts", default = NA,
                       help = "file with one perturbation per line: fit only these (e.g. to check a subset)")

parsed <- parse_args(parser)

# Directory stuff
################################################################################

expr_f <- parsed$expr
mdat_f <- parsed$mdat
workers <- as.numeric(parsed$workers)
out_dir <- parsed$out_dir
reader <- parsed$reader
if (!reader %in% c("fread", "lines")) stop("--reader must be fread or lines")

if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
}

# Load the data
################################################################################

# Features with no variance among their non-missing values (constant, or all
# missing) are dropped. Features with some missing values are kept; limma fits
# them on the samples they have.

# --reader lines: the CSV read in blocks of lines with scan(), each block
# filtered before the next one (memory: the kept rows plus one block).
# data.table's fread crashes (segfault) on the MetChange matrix (~150k columns,
# lines of ~3 MB), with or without the column types given. scan() parses the
# values straight into numbers (strsplit + as.numeric is ~7x slower: one R string
# per value), stops on lines with another number of fields or non-numeric
# values, and reads empty fields as NA.
read_matrix_lines <- function(f, chunk = 50L) {
    con <- file(f, open = "r")
    on.exit(close(con))
    header <- gsub('^"|"$', "", strsplit(readLines(con, n = 1L), ",", fixed = TRUE)[[1]])
    n <- length(header) - 1L
    what <- c(list(""), rep(list(0), n))
    blocks <- list()
    n_read <- 0L
    n_na <- 0L
    repeat {
        r <- scan(con, what = what, sep = ",", nmax = chunk, quiet = TRUE, multi.line = FALSE)
        k <- length(r[[1]])
        if (k == 0L) break
        X <- matrix(unlist(r[-1], use.names = FALSE), nrow = k, dimnames = list(r[[1]], NULL))
        rm(r)
        row_min <- suppressWarnings(apply(X, 1, min, na.rm = TRUE))
        row_max <- suppressWarnings(apply(X, 1, max, na.rm = TRUE))
        keep <- row_max > row_min
        n_na <- n_na + sum(keep & rowSums(is.na(X)) > 0)
        blocks[[length(blocks) + 1L]] <- X[keep, , drop = FALSE]
        n_read <- n_read + k
        if (length(blocks) %% 20L == 0L) message(sprintf("  %d features read", n_read))
    }
    expr <- do.call(rbind, blocks)
    colnames(expr) <- header[-1]
    message(sprintf("%d features, %d with no variance or all values missing removed, %d kept with some missing values",
                    n_read, n_read - nrow(expr), n_na))
    expr
}

message(sprintf("Loading expression data from %s", expr_f))
if (reader == "lines") {
    expr <- read_matrix_lines(expr_f)
} else {
    # One numeric matrix, filled from the data.table block by block while the
    # copied columns are dropped: a data.frame with one column per sample, apply()
    # and row subsetting held several full copies at once. Column types are given,
    # not guessed.
    header <- names(fread(expr_f, header = TRUE, nrows = 0))
    dt <- fread(expr_f, header = TRUE, colClasses = c("character", rep("double", length(header) - 1)))
    feature_ids <- dt[[1]]
    sample_ids <- names(dt)[-1]

    # min and max per row, accumulated over the columns
    row_min <- rep(Inf, nrow(dt))
    row_max <- rep(-Inf, nrow(dt))
    n_missing <- rep(0L, nrow(dt))
    for (k in seq_along(sample_ids) + 1L) {
        x <- dt[[k]]
        row_min <- pmin(row_min, x, na.rm = TRUE)
        row_max <- pmax(row_max, x, na.rm = TRUE)
        n_missing <- n_missing + is.na(x)
    }
    keep <- row_max > row_min
    message(sprintf("%d features, %d with no variance or all values missing removed, %d kept with some missing values",
                    length(keep), sum(!keep), sum(keep & n_missing > 0)))

    expr <- matrix(NA_real_, nrow = sum(keep), ncol = length(sample_ids),
                   dimnames = list(feature_ids[keep], sample_ids))
    block <- 1000
    for (start in seq(1, length(sample_ids), by = block)) {
        idx <- start:min(start + block - 1, length(sample_ids))
        # the columns before this block were dropped: sample idx[1] is column 2
        for (k in idx) expr[, k] <- dt[[k - start + 2L]][keep]
        dt[, (sample_ids[idx]) := NULL]
    }
    rm(dt, row_min, row_max, n_missing, keep, x)
}
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
if (!is.na(parsed$perts)) {
    wanted <- readLines(parsed$perts)
    if (length(setdiff(wanted, perts))) stop(sprintf("%d perturbations of --perts not in the metadata",
                                                     length(setdiff(wanted, perts))))
    perts <- wanted
    message(sprintf("--perts: fitting %d perturbations", length(perts)))
}
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
# per-gene gamma GLM (glmgam.fit) is the same. Genes with missing values go
# through limma (each is fitted on its own samples) and the consensus is taken
# over all genes, as in limma. Other cases go to limma entirely.
dup_cor_fast <- function(M, design, block, trim = 0.15, chunk = 2000L) {
    design <- as.matrix(design)
    QR <- qr(design)
    if (max(table(block)) < 2L || QR$rank < ncol(design)) {
        return(duplicateCorrelation(M, design, block = block, trim = trim))
    }
    incomplete <- rowSums(is.na(M)) > 0
    if (all(incomplete)) return(duplicateCorrelation(M, design, block = block, trim = trim))
    if (any(incomplete)) {
        # limma clips each gene's correlation to the same bounds before the
        # trimmed mean, so the two parts can be computed apart and pooled
        arho <- rep_len(NA_real_, nrow(M))
        arho[!incomplete] <- dup_cor_fast(M[!incomplete, , drop = FALSE], design, block, trim, chunk)$atanh.correlations
        arho[incomplete] <- duplicateCorrelation(M[incomplete, , drop = FALSE], design, block = block,
                                                 trim = trim)$atanh.correlations
        mrho <- tanh(mean(arho, trim = trim, na.rm = TRUE))
        return(list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho))
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

# lmFit, with the genes that have missing values fitted apart: with a single
# missing value limma fits every gene one by one (with block, a Cholesky
# decomposition of samples x samples per gene) instead of all genes together.
# Each gene's fit does not depend on the others, so the per-gene results are
# the same (up to rounding) and are put back in the input order.
lm_fit_split <- function(M, design, ...) {
    incomplete <- rowSums(is.na(M)) > 0
    if (!any(incomplete) || all(incomplete)) return(lmFit(M, design, ...))
    fit <- lmFit(M[!incomplete, , drop = FALSE], design, ...)
    fit_na <- lmFit(M[incomplete, , drop = FALSE], design, ...)
    ord <- order(c(which(!incomplete), which(incomplete)))
    for (f in c("coefficients", "stdev.unscaled")) fit[[f]] <- rbind(fit[[f]], fit_na[[f]])[ord, , drop = FALSE]
    for (f in c("sigma", "df.residual", "Amean")) fit[[f]] <- c(fit[[f]], fit_na[[f]])[ord]
    fit
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
    # Controls are the reference level explicitly: with ~ perturbation the reference is whichever name sorts
    # first, which depends on the locale (en_US: "dacinostat_..." < "DMSO_..." < "vorinostat_..."), and
    # perturbations sorting before the controls got control - treated as log2FC and their own mean as basal.
    group <- factor(ifelse(mdat_filt$exper_group == "ctrl", "ctrl", "pert"), levels = c("ctrl", "pert"))
    design <- model.matrix(~ group)

    if (info$n_batches > 1) {
        corfit <- dup_cor_fast(dat_filt, design, block = mdat_filt$batch)
        info$correlation <- corfit$consensus.correlation
        fit <- lm_fit_split(dat_filt, design, block = mdat_filt$batch,
                            correlation = corfit$consensus.correlation)
    } else {
        fit <- lm_fit_split(dat_filt, design)
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
# marker: these outputs use the controls as the reference level (fix_limma_reference_sign.py refuses to touch them)
writeLines("reference level: ctrl (explicit)", file.path(out_dir, "reference_level.txt"))

message(sprintf("basal.csv, log2FC.csv, pvals.csv, qvals.csv and fit_info.csv saved at %s", out_dir))
message(sprintf("The checkpoints in %s can be deleted.", ckpt_dir))
