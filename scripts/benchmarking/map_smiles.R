# With the perturbation metadata of two datasets and indications of columns,
# canonicalize smiles columns and find intersecting compounds. prt_df is
# intended to be our database, and val_df whatever validation dataset we might
# find (e.g. Attene-Ramos or Hallinger).
#
# Matching, in order:
#   1. salts / mixtures are reduced to their largest fragment ("CC(=O)O.Cl" ->
#      "CC(=O)O"), so a hydrochloride on one side matches the free base on the
#      other. Validation sets are full of these (Attene-Ramos: 2,040 of 10,337
#      SMILES carry a ".", LINCS: 143 of 34,419), and without this they are all
#      lost;
#   2. --stereo exact: canonical SMILES *with* stereochemistry, so enantiomers
#      stay apart;
#      --stereo flat: canonical SMILES without stereochemistry, so they merge;
#      --stereo both (default): pass 1 exact, then the compounds left unmatched
#      are retried without stereochemistry. Every pair is labelled in
#      match_type (exact / stereo_insensitive) and the stereo-insensitive ones
#      carry n_pert_forms / n_val_forms: how many distinct stereo forms that
#      flat structure has on each side (> 1 = the pair is ambiguous, e.g. a
#      racemate matched to one enantiomer). Decide per analysis whether to keep
#      them instead of having the choice made silently.
# All matching pairs are kept: if one structure appears under several ids (the
# same compound under several BRD ids, or several wells / CIDs), every pair is
# written, not an arbitrary one (the previous match() version silently dropped
# the rest). Unparseable SMILES are reported and dropped, never matched to each
# other as NA.
#
# Output in --out_dir: mapped_drugs.csv (the pairs) and mapping_summary.txt
# (counts: parsed, matched, by match type, and how many are ambiguous).

if (!require("rcdk", quietly = TRUE)) install.packages("rcdk")
library(rcdk)
library(data.table)
library(argparser)

# Functions
################################################################################

# Largest fragment of a multi-component SMILES, by number of atoms: drops
# counterions and co-formulants. Ties keep the first.
largest_fragment <- function(smiles) {
        if (is.na(smiles) || !nzchar(smiles)) return(NA_character_)
        parts <- strsplit(smiles, ".", fixed = TRUE)[[1]]
        if (length(parts) == 1) return(parts[1])
        sizes <- vapply(parts, function(p) {
                mol <- tryCatch(parse.smiles(p)[[1]], error = function(e) NULL)
                if (is.null(mol)) 0L else length(get.atoms(mol))
        }, integer(1))
        if (all(sizes == 0)) return(NA_character_)
        parts[which.max(sizes)]
}

# Canonical SMILES of the largest fragment, with or without stereochemistry.
# CDK's Canonical flavor writes canonical SMILES without stereo / isotopes;
# adding Isomeric keeps them.
canonicalize <- function(smiles, stereo) {
        frag <- largest_fragment(smiles)
        if (is.na(frag)) return(NA_character_)
        mol <- tryCatch(parse.smiles(frag)[[1]], error = function(e) NULL)
        if (is.null(mol)) return(NA_character_)
        flavor <- if (stereo) smiles.flavors(c("Canonical", "Isomeric")) else
                smiles.flavors(c("Canonical"))
        out <- tryCatch(get.smiles(mol, flavor), error = function(e) NA_character_)
        if (is.null(out) || !nzchar(out)) NA_character_ else out
}

# Canonicalizes the unique values only (validation tables repeat SMILES over
# wells / concentrations), then maps back.
canonicalize_smiles_vec <- function(smiles_vec, stereo, label) {
        uniq <- unique(smiles_vec[!is.na(smiles_vec) & nzchar(smiles_vec)])
        cat(sprintf("%s: canonicalizing %d unique SMILES (%s stereochemistry)\n",
                    label, length(uniq), if (stereo) "with" else "without"))
        pb <- txtProgressBar(min = 0, max = max(length(uniq), 1), style = 3)
        can <- rep(NA_character_, length(uniq))
        for (i in seq_along(uniq)) {
                setTxtProgressBar(pb, i)
                can[i] <- canonicalize(uniq[i], stereo)
        }
        close(pb)
        unname(setNames(can, uniq)[smiles_vec])
}

# Pairs of rows sharing a key, all combinations, NA keys dropped.
pairs_on <- function(prt, val, key) {
        a <- prt[!is.na(get(key)), .(pert_row = .I, k = get(key))]
        b <- val[!is.na(get(key)), .(val_row = .I, k = get(key))]
        merge(a, b, by = "k", allow.cartesian = TRUE)
}

# Terminal argument parser
################################################################################
parser <- arg_parser("Map a validation dataset to our perturbations by structure.")

parser <- add_argument(
        parser = parser,
        arg = c(
                "--pert",
                "--val",
                "--pert_smiles_col",
                "--pert_pertid_col",
                "--pert_pertname_col",
                "--val_smiles_col",
                "--val_pertid_col",
                "--val_name_col",
                "--stereo",
                "--out_dir"
        ),
        help = c(
                "Training perturbation metadata file (tsv).",
                "Validation perturbation metadata file (csv).",
                "Training metadata column holding SMILES.",
                "Training metadata column holding perturbation ID.",
                "Training metadata column holding perturbation name.",
                "Validation metadata column holding SMILES.",
                "Validation metadata column holding perturbation ID.",
                "Validation metadata column holding a compound name (optional).",
                "exact (keep enantiomers apart), flat (merge them) or both (exact first, then the rest without stereochemistry; default).",
                "Output directory where placing the results."
        ),
        default = list(NA, NA, NA, NA, NA, NA, NA, NA, "both", NA),
        flag = c(F, F, F, F, F, F, F, F, F, F)
)

parsed <- parse_args(parser)

val_df_f <- parsed$val
prt_df_f <- parsed$pert

val_smiles_col <- parsed$val_smiles_col
prt_smiles_col <- parsed$pert_smiles_col
prt_pert_id_col <- parsed$pert_pertid_col
prt_pertname_col <- parsed$pert_pertname_col
val_pertid_col <- parsed$val_pertid_col
val_name_col <- parsed$val_name_col
stereo_mode <- parsed$stereo
out_dir <- parsed$out_dir

if (!stereo_mode %in% c("exact", "flat", "both")) {
        stop("--stereo must be exact, flat or both")
}
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Data loading
################################################################################

val_df <- as.data.table(read.csv(val_df_f, sep = ','))
prt_df <- as.data.table(read.csv(prt_df_f, sep = '\t'))

prt <- data.table(
        pertid_dataset = prt_df[[prt_pert_id_col]],
        pertname_dataset = prt_df[[prt_pertname_col]],
        smiles_dataset = prt_df[[prt_smiles_col]]
)
if ("inchi_key" %in% names(prt_df)) prt[, inchi_key := prt_df[["inchi_key"]]]
val <- data.table(
        val_id = val_df[[val_pertid_col]],
        smiles_val = val_df[[val_smiles_col]]
)
if (!is.na(val_name_col) && val_name_col %in% names(val_df)) {
        val[, val_name := val_df[[val_name_col]]]
}
# one row per compound per side: the validation table repeats wells / doses
prt <- unique(prt, by = c("pertid_dataset", "smiles_dataset"))
val <- unique(val, by = c("val_id", "smiles_val"))

# Canonicalize
################################################################################

need_exact <- stereo_mode %in% c("exact", "both")
need_flat <- stereo_mode %in% c("flat", "both")
if (need_exact) {
        prt[, smi_exact := canonicalize_smiles_vec(smiles_dataset, TRUE, "perturbations")]
        val[, smi_exact := canonicalize_smiles_vec(smiles_val, TRUE, "validation")]
}
if (need_flat) {
        prt[, smi_flat := canonicalize_smiles_vec(smiles_dataset, FALSE, "perturbations")]
        val[, smi_flat := canonicalize_smiles_vec(smiles_val, FALSE, "validation")]
}
parse_col <- if (need_exact) "smi_exact" else "smi_flat"
n_prt_ok <- prt[!is.na(get(parse_col)), .N]
n_val_ok <- val[!is.na(get(parse_col)), .N]

# Match
################################################################################

res <- list()
matched_val <- character(0)
if (need_exact) {
        p <- pairs_on(prt, val, "smi_exact")
        if (nrow(p)) {
                res[["exact"]] <- data.table(
                        pertid_dataset = prt$pertid_dataset[p$pert_row],
                        pertname_dataset = prt$pertname_dataset[p$pert_row],
                        val_id = val$val_id[p$val_row],
                        val_name = if ("val_name" %in% names(val)) val$val_name[p$val_row] else NA,
                        smiles = p$k,
                        inchi_key = if ("inchi_key" %in% names(prt)) prt$inchi_key[p$pert_row] else NA,
                        match_type = "exact",
                        n_pert_forms = 1L,
                        n_val_forms = 1L
                )
                matched_val <- unique(as.character(val$val_id[p$val_row]))
        }
}
if (need_flat) {
        # in "both" mode only the validation compounds still unmatched are retried
        val_left <- if (stereo_mode == "both") val[!as.character(val_id) %in% matched_val] else val
        p <- pairs_on(prt, val_left, "smi_flat")
        if (nrow(p)) {
                # how many distinct stereo forms share this flat structure, per side
                forms_p <- if (need_exact) prt[!is.na(smi_flat), .(n_pert_forms = uniqueN(smi_exact)), by = smi_flat] else
                        prt[!is.na(smi_flat), .(n_pert_forms = 1L), by = smi_flat]
                forms_v <- if (need_exact) val_left[!is.na(smi_flat), .(n_val_forms = uniqueN(smi_exact)), by = smi_flat] else
                        val_left[!is.na(smi_flat), .(n_val_forms = 1L), by = smi_flat]
                d <- data.table(
                        pertid_dataset = prt$pertid_dataset[p$pert_row],
                        pertname_dataset = prt$pertname_dataset[p$pert_row],
                        val_id = val_left$val_id[p$val_row],
                        val_name = if ("val_name" %in% names(val_left)) val_left$val_name[p$val_row] else NA,
                        smiles = p$k,
                        inchi_key = if ("inchi_key" %in% names(prt)) prt$inchi_key[p$pert_row] else NA,
                        match_type = if (stereo_mode == "flat") "flat" else "stereo_insensitive"
                )
                d <- merge(d, forms_p, by.x = "smiles", by.y = "smi_flat", all.x = TRUE)
                d <- merge(d, forms_v, by.x = "smiles", by.y = "smi_flat", all.x = TRUE)
                res[[length(res) + 1L]] <- d
        }
}

map_df <- if (length(res)) rbindlist(res, use.names = TRUE, fill = TRUE) else
        data.table(pertid_dataset = character(0))
if (nrow(map_df)) {
        map_df[, ambiguous_stereo := !is.na(n_pert_forms) & !is.na(n_val_forms) &
                       (n_pert_forms > 1 | n_val_forms > 1) & match_type != "exact"]
        setorder(map_df, match_type, pertname_dataset)
}

write.csv(map_df, file = file.path(out_dir, "mapped_drugs.csv"), row.names = FALSE)

# Summary
################################################################################

lines <- c(
        sprintf("perturbations: %d compounds, %d parsed (%d failed)",
                nrow(prt), n_prt_ok, nrow(prt) - n_prt_ok),
        sprintf("validation:    %d compounds, %d parsed (%d failed)",
                nrow(val), n_val_ok, nrow(val) - n_val_ok),
        sprintf("stereo mode:   %s", stereo_mode),
        sprintf("pairs written: %d (%d perturbations, %d validation compounds)",
                nrow(map_df),
                if (nrow(map_df)) uniqueN(map_df$pertid_dataset) else 0L,
                if (nrow(map_df)) uniqueN(map_df$val_id) else 0L)
)
if (nrow(map_df)) {
        for (mt in unique(map_df$match_type)) {
                s <- map_df[match_type == mt]
                lines <- c(lines, sprintf("  %-18s %d pairs, %d validation compounds%s",
                                          mt, nrow(s), uniqueN(s$val_id),
                                          if (mt != "exact") sprintf(", %d with an ambiguous stereo form",
                                                                     sum(s$ambiguous_stereo, na.rm = TRUE)) else ""))
        }
        n_multi <- map_df[, .N, by = val_id][N > 1, .N]
        lines <- c(lines, sprintf("validation compounds matching more than one perturbation: %d", n_multi))
}
writeLines(lines, file.path(out_dir, "mapping_summary.txt"))
cat(paste(lines, collapse = "\n"), "\n")
