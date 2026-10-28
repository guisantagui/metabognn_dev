# Extract perturbation data from A549 lung cancer cell line to assess if 
# the inferred metabolic alterations obtained with LINCS + GEMs constraining
# match the observed in the study by Schuhknecht et al.
# (https://doi.org/10.1038/s41587-024-02524-5)

if(!require(argparser, quietly = T)){
        install.packages("argparser", repos='http://cran.us.r-project.org')
}
library(argparser)
if(!require(BiocManager, quietly = T)){
        install.packages("BiocManager", repos='http://cran.us.r-project.org')
}
if(!require("cmapR", quietly = T)) BiocManager::install("cmapR", update = F)
library(cmapR)
library(plotUtils)

writeCsvFst <- function(df, file, rowNames = T, colNames = T){
        if(rowNames){
                rn <- rownames(df)
                df <- data.table::data.table(df)
                df[, V1 := rn]
                data.table::setcolorder(df, c("V1", setdiff(names(df), "V1")))
        }else{
                df <- data.table::data.table(df)
        }
        data.table::fwrite(df, file, col.names = colNames)
}

# Terminal argument parser
################################################################################
parser <- arg_parser("Extract data from a given cell type into csv files.")

parser <- add_argument(
        parser = parser,
        arg = c(
                "--gctx_pert",
                "--gctx_ctrl",
                "--cell",
                "--ctrls_keep",
                "--drugs",
                "--out_dir"
        ),
        help = c(
                "Perturbation GCTX file.",
                "Perturbation GCTX file.",
                "The cell type to be extracted.",
                "Optional argument, one or more control types, separated by commas, to further filter the controls (e.g. DMSO,noPert)",
                "Optional argument, a CSV file with drugs to be retained (instead of keeping all).",
                "Output directory where placing the results."
        ),
        flag = c(F, F, F, F, F, F)
)

parsed <- parse_args(parser)

# Directory stuff
################################################################################
gctxFile_pert <- parsed$gctx_pert
gctxFile_ctrl <- parsed$gctx_ctrl
cell <- parsed$cell
drugs_f <- parsed$drugs
ctrl_keep <- parsed$ctrls_keep
out_dir <- parsed$out_dir

#gctxFile_pert <- "/scratch/users/gsantamaria/test_large_files/data/perturbation/lincs/cp_predicted_RNAseq_profiles.gctx"
#gctxFile_ctrl <- "/scratch/users/gsantamaria/test_large_files/data/perturbation/lincs/ctl_predicted_RNAseq_profiles.gctx"
#cell <- "A549"
#ctrls_keep <- "DMSO"
#drugs_f <- "/home/users/gsantamaria/projects/metabo_pert/data/LINCS/common_drugs_lincs_info.csv"
#out_dir <- "/home/users/gsantamaria/projects/metabo_pert/data/LINCS"

if (!is.na(ctrls_keep)){
        ctrls_keep <- strsplit(ctrls_keep, split = ",")
}


# Load data
################################################################################

cell_info <- read.table(
    cell_info_f,
    sep = "\t",
    header = TRUE,
    quote = "",
    fill = TRUE,
    comment.char = ""
)

message("Loading perturbation metadata...")
metDat_pert <- read_gctx_meta(gctxFile_pert, dim = "column")
message("Done")

message("Loading control metadata...")
metDat_ctrl <- read_gctx_meta(gctxFile_ctrl, dim = "column")
message("Loading control metadata...")

# Filter data
################################################################################

# Perturbation
if (cell %in% metDat_pert$cell){
        bool_cell <- metDat_pert$cell == cell
}

if (!is.na(drugs_f)){
        drugs_df <- read.csv(drugs_f, row.names = 1)
        drugs <- drugs_df$pert
        bool_drug <- metDat_pert$pertname %in% drugs
}else{
        bool_drug <- rep(TRUE, nrow(metDat_pert))
}

idxs_pert <- which(bool_cell & bool_drug)

gctx_pert_filt <- parse_gctx(gctxFile_pert, cid = idxs_pert)
metDat_pert_filt <- metDat_pert[idxs_pert, ]

# Add batch column, using same criteria as Ma'ayan
metDat_pert_filt$batch <- sapply(
        metDat_pert_filt$id,
        function(x) paste(strsplit(x, split = '_')[[1]][1:3], collapse = '_')
)


# Controls

# Add batch column to ctrl metadata too
metDat_ctrl$batch <- sapply(
        metDat_ctrl$id,
        function(x) paste(strsplit(x, split = '_')[[1]][1:3], collapse = '_')
)

bool_batch <- metDat_ctrl$batch %in% metDat_pert_filt$batch

if (!is.na(ctrls_keep)){
        ctrls_keep_not_in_data <- ctrls_keep[!ctrls_keep %in% metDat_ctrl$pertname]
        if (length(ctrls_keep_not_in_data) > 0){
                ctrls_keep_not_in_data <- paste(ctrls_keep_not_in_data, collapse = ", ")
                warning(
                        sprintf(
                                "%s not found in control metadata's pertname column. They won't be included in the output.",
                                ctrls_keep_not_in_data),
                        call. = FALSE
                )
        }
        ctrls_keep_in_data <- ctrls_keep[ctrls_keep %in% metDat_ctrl$pertname]
        if (length(ctrls_keep_in_data) == 0){
                warning("No ctrls_keep was found in control metadata's pertname column. All controls will be kept.")
                ctrls_keep_bool <- rep(TRUE, nrow(metDat_ctrl))
        }else{
                message(sprintf("Keeping only %s controls", paste(ctrls_keep_in_data, collapse = ", ")))
                ctrls_keep_bool <- metDat_ctrl$pertname %in% ctrls_keep_in_data
        }
}else{
        ctrls_keep_bool <- rep(TRUE, nrow(metDat_ctrl))
}


idxs_ctrl <- which(bool_batch & ctrls_keep_bool)
gctx_ctrl_filt <- parse_gctx(gctxFile_ctrl, cid = idxs_ctrl)

# Generate output files
################################################################################

# Merge controls and perts
mat_pert <- gctx_pert_filt@mat
mdat_pert <- gctx_pert_filt@cdesc

mdat_pert$batch <- sapply(
        gctx_pert_filt@cdesc$id,
        function(x) paste(
                strsplit(
                        x, split = "_"
                )[[1]][1:3],
                collapse = "_"
        )
)

mdat_pert$exper_group <- "pert"

mat_ctrl <- gctx_ctrl_filt@mat
mdat_ctrl <- gctx_ctrl_filt@cdesc
mdat_ctrl$batch <- sapply(
        gctx_ctrl_filt@cdesc$id,
        function(x) paste(
                strsplit(
                        x, split = "_"
                )[[1]][1:3],
                collapse = "_"
        )
)
mdat_ctrl$exper_group <- "ctrl"

pert_batches <- unique(mdat_pert$batch)
ctrl_batches <- unique(mdat_ctrl$batch)

# Remove samples from perturbation data that belong to batches without controls
batches_wo_ctrls <- pert_batches[!pert_batches %in% ctrl_batches]
keep_samps_good_batches <- !mdat_pert$batch %in% batches_wo_ctrls
mat_pert <- mat_pert[, keep_samps_good_batches]
mdat_pert <- mdat_pert[keep_samps_good_batches, ]

mat_all <- cbind(
        mat_pert,
        mat_ctrl
)

cols_2add <- colnames(mdat_pert)[
        !colnames(mdat_pert) %in% colnames(mdat_ctrl)
]

# Add columns missing in ctrl metadata
cols_2add <- data.frame(
        matrix(
                nrow = nrow(mdat_ctrl),
                ncol = length(cols_2add),
                dimnames = list(rownames(mdat_ctrl), cols_2add)
        )
)

mdat_ctrl <- cbind.data.frame(mdat_ctrl, cols_2add)[, colnames(mdat_pert)]

mdat_all <- rbind.data.frame(
        mdat_pert,
        mdat_ctrl
)

mdat_all$cell_type_fine <- sapply(
        mdat_all$batch,
        function(x) strsplit(x, split = "_")[[1]][2]
)

mdat_all$perturbation <- paste(
        mdat_all$pertname,
        mdat_all$dose,
        mdat_all$timepoint,
        sep = "_"
)

# Save results
writeCsvFst(mat_all, file = file.path(out_dir, "expr.csv"))
writeCsvFst(mdat_all, file = file.path(out_dir, "mdat.csv"))