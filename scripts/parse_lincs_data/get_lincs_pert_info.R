# Extract drug metadata from lincs perturbations

if(!require(argparser, quietly = T)){
        install.packages("argparser", repos='http://cran.us.r-project.org')
}
library(argparser)
if(!require(BiocManager, quietly = T)){
        install.packages("BiocManager", repos='http://cran.us.r-project.org')
}
if(!require("cmapR", quietly = T)) BiocManager::install("cmapR", update = F)
library(cmapR)

writeCsvFst <- function(df,
        file, rowNames = T, colNames = T,
        sep=getOption("datatable.fwrite.sep", ","), quote = "auto"
){
        if(rowNames){
                rn <- rownames(df)
                df <- data.table::data.table(df)
                df[, V1 := rn]
                data.table::setcolorder(df, c("V1", setdiff(names(df), "V1")))
        }else{
                df <- data.table::data.table(df)
        }
        data.table::fwrite(df, file, col.names = colNames, sep = sep, quote = quote)
}

# Terminal argument parser
################################################################################
parser <- arg_parser("Extract data from a given cell type into csv files.")

parser <- add_argument(
        parser = parser,
        arg = c(
                "--gctx_pert",
                "--cell",
                "--out_dir"
        ),
        help = c(
                "Perturbation GCTX file.",
                "Cell type, if a filtered version per cell is desired.",
                "Output directory where placing the results."
        ),
        flag = c(F, F, F)
)

parsed <- parse_args(parser)

# Directory stuff
################################################################################
pert_f <- parsed$gctx_pert
cell <- parsed$cell
out_dir <- parsed$out_dir

if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
}

# Load data
################################################################################

message("Loading perturbation metadata...")
metDat_pert <- read_gctx_meta(pert_f, dim = "column")
message("Done")

if(!is.na(cell) & cell %in% metDat_pert$cell){
        metDat_pert <- metDat_pert[metDat_pert$cell == cell, ]
        out_name <- sprintf(file.path(out_dir, "lincs_unique_perts_%s.tsv"), cell)
}else{
        out_name <- file.path(out_dir, "lincs_unique_perts.tsv")
}


out <- unique(
        metDat_pert[, c("pert_id", "pertname", "smiles", "inchi_key")]
)

rownames(out) <- 1:nrow(out)

writeCsvFst(
        out,
        file = out_name,
        sep = '\t',
        quote = TRUE
)

print(sprintf("%s saved in %s.", basename(out_name), out_dir))