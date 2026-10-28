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
expr <- as.data.frame(fread(expr_f, header = TRUE))
rownames(expr) <- expr[, 1]
expr <- expr[, 2:ncol(expr)]

message(sprintf("Loading and parsing metadata from %s", mdat_f))
mdat <- as.data.frame(fread(mdat_f, header = TRUE))
#mdat$batch <- sapply(
#        mdat$rna_plate,
#        function(x) paste(strsplit(x, split = "_")[[1]][1:3], collapse = "_")
#)

perts <- unique(mdat$perturbation[mdat$exper_group != "ctrl"])

# Do parallelized limma analysis
################################################################################

# Remove features with no variance
expr <- expr[apply(expr, 1, var) != 0, ]

options(future.globals.maxSize = 10 * 1024^3)
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

            results <- topTable(
                fit,
                coef = 2,
                number = Inf,
                adjust.method = "BH",
                confint = TRUE
            )
            basal <- fit$coefficients[, 1]
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

# Convert to matrices/data.frames
basal <- do.call(
    cbind,
    lapply(res, function(x) {
        setNames(x$basal, x$genes)
    })
)

log2FCs <- do.call(
    cbind,
    lapply(res, function(x) {
        setNames(x$log2FC, x$genes)
    })
)

pvals <- do.call(
    cbind,
    lapply(res, function(x) {
        setNames(x$pval, x$genes)
    })
)

qvals <- do.call(
    cbind,
    lapply(res, function(x) {
        setNames(x$qval, x$genes)
    })
)

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