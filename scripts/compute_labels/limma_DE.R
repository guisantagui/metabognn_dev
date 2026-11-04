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

# Do parallelized limma analysis
################################################################################

options(future.globals.maxSize = 30 * 1024^3)
plan(multicore, workers = workers)

message("Running limma...")
res <- future_lapply(
    perts,
    function(pert) {
        pert_bool <- mdat$perturbation == pert
        mdat_pert <- mdat[pert_bool, ]

        batches <- unique(mdat_pert$batch)
        timepoint <- unique(mdat_pert$timepoint)

        ctrl_bool <- (
            mdat$exper_group == "ctrl" &
            mdat$timepoint == timepoint &
            mdat$batch %in% batches
        )

        mdat_filt <- mdat[pert_bool | ctrl_bool, ]
        dat_filt <- expr[, pert_bool | ctrl_bool, drop = FALSE]

        if (length(unique(mdat_filt$perturbation)) > 1){
            design <- model.matrix(~ perturbation, data = mdat_filt)

            if (length(unique(mdat_filt$batch)) > 1) {

                corfit <- duplicateCorrelation(
                    dat_filt,
                    design,
                    block = mdat_filt$batch
                )

                fit <- lmFit(
                    dat_filt,
                    design,
                    block = mdat_filt$batch,
                    correlation = corfit$consensus.correlation
                )

            } else {
                fit <- lmFit(dat_filt, design)
            }

            fit <- eBayes(fit)

            # sort.by = "none": keep the input gene order (topTable sorts by B by
            # default, which differs between perturbations)
            results <- topTable(
                fit,
                coef = 2,
                number = Inf,
                adjust.method = "BH",
                confint = TRUE,
                sort.by = "none"
            )
            basal <- fit$coefficients[rownames(results), 1]
            #delta <- fit$coefficients[, 2]

            list(
                basal = basal,
                #delta = delta,
                log2FC = results$logFC,
                pval   = results$P.Value,
                qval   = results$adj.P.Val,
                genes  = rownames(results)
            )
        } else {
            list(
                basal = setNames(rep(NA_real_, nrow(expr)), rownames(expr)),
                #delta = setNames(rep(NA_real_, nrow(expr)), rownames(expr)),
                log2FC = setNames(rep(NA_real_, nrow(expr)), rownames(expr)),
                pval   = setNames(rep(NA_real_, nrow(expr)), rownames(expr)),
                qval   = setNames(rep(NA_real_, nrow(expr)), rownames(expr)),
                genes  = rownames(expr)
            )
        }
    },
    future.seed = TRUE
)
message("limma done.")

names(res) <- perts

# Convert to matrices/data.frames. Each perturbation's values are matched to one
# gene order by name: cbind() of named vectors aligns by position, not by name.
genes_all <- rownames(expr)
to_matrix <- function(field) {
    m <- vapply(res, function(x) x[[field]][match(genes_all, x$genes)], numeric(length(genes_all)))
    rownames(m) <- genes_all
    m
}
basal <- to_matrix("basal")
log2FCs <- to_matrix("log2FC")
pvals <- to_matrix("pval")
qvals <- to_matrix("qval")

colnames(basal) <- perts
colnames(log2FCs) <- perts
colnames(pvals) <- perts
colnames(qvals) <- perts

write_table_fast(basal, f = file.path(out_dir, "basal.csv"))
write_table_fast(log2FCs, f = file.path(out_dir, "log2FC.csv"))
write_table_fast(pvals, f = file.path(out_dir, "pvals.csv"))
write_table_fast(qvals, f = file.path(out_dir, "qvals.csv"))

message(sprintf("basal.csv saved at %s", out_dir))
message(sprintf("log2FC.csv saved at %s", out_dir))
message(sprintf("pvals.csv saved at %s", out_dir))
message(sprintf("qvals.csv saved at %s", out_dir))