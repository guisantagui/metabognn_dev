################################################################################
# Is the per-drug sign of the MetChange / metabolomics association real?
################################################################################
# Run AFTER metchange_vs_expression_benchmark.R in the same R session. Needs
# from that script: dtp, auc_fn, auc_dir, alpha_sig, out_dir, and the loaded
# data.table / ggplot2 packages.
#
# HYPOTHESIS UNDER TEST
#   Per-drug association between MetChange (production consistency) and measured
#   metabolite pools has a mechanism-dependent sign: for some drugs a pool rises
#   BECAUSE flux is blocked, giving anti-correlation. If true:
#     (1) per-drug rho is more spread out than a mismatched-pairs null,
#     (2) the sign is stable on independent metabolites (split-half),
#     (3) the sign is predictable from drug properties on held-out drugs.
#   (1) and (2) check it is not noise; (3) is what would make it usable.
#
# NOTE: written without R access, not run. Check shapes the first time.
################################################################################
stopifnot(exists("dtp"), exists("auc_fn"), exists("auc_dir"))
set.seed(11)

# A. WHY THE "AUC > 0.5 / AUC < 0.5" SUBGROUP TESTS ARE INVALID ################
# Splitting drugs on the MetChange AUC and then comparing MetChange with
# expression in each group selects on the outcome. With no signal at all, the
# group with mc > 0.5 shows mc > exp and the group with mc < 0.5 shows mc < exp.
sim <- data.table(mc = rnorm(1e4, 0.5, 0.15), exp = rnorm(1e4, 0.5, 0.15))
cat("Null demo, mean(mc - exp) when mc > 0.5:", round(sim[mc > 0.5, mean(mc - exp)], 3),
    "| when mc < 0.5:", round(sim[mc < 0.5, mean(mc - exp)], 3), "\n")
# A valid version classifies drugs on one set of metabolites and tests on
# another (section D).


# B. MISMATCHED-PAIRS NULL ON PER-DRUG rho #####################################
# Matched rho(i,i): mc profile of drug i vs measured profile of drug i.
# Null rho(i,j), j != i: mc profile of drug i vs measured profile of ANOTHER
# drug. Keeps each profile's internal structure (correlated metabolites), unlike
# shuffling metabolites, which would give a null that is too narrow.
# Uses rho over all ~280 metabolites, far less noisy than an AUC computed from
# the 3-10 metabolites that changed.
to_wide <- function(dt, v) {
        w <- dcast(dt, metabolite ~ perturbation, value.var = v)
        m <- as.matrix(w[, -1, with = FALSE])
        rownames(m) <- w$metabolite
        m
}
MET <- to_wide(dtp, "met_c")
MC <- to_wide(dtp, "mc_c")
EX <- to_wide(dtp, "exp_c")
stopifnot(identical(colnames(MET), colnames(MC)), identical(colnames(MET), colnames(EX)))

mismatch_test <- function(S, MET, label) {
        R <- cor(S, MET, method = "spearman", use = "pairwise.complete.obs")  # rows: score of drug i
        n <- nrow(R)
        matched <- diag(R)
        null <- R[row(R) != col(R)]
        p <- sapply(seq_len(n), function(i) (sum(abs(R[i, -i]) >= abs(R[i, i]), na.rm = TRUE) + 1) / n)
        q <- p.adjust(p, "BH")
        list(
                matched = matched, null = null, p = p, q = q,
                summary = data.frame(
                        score = label, n_drugs = n,
                        mean_matched = mean(matched, na.rm = TRUE),
                        mean_abs_matched = mean(abs(matched), na.rm = TRUE),
                        mean_abs_null = mean(abs(null), na.rm = TRUE),
                        sd_matched = sd(matched, na.rm = TRUE), sd_null = sd(null, na.rm = TRUE),
                        n_p05_pos = sum(p < 0.05 & matched > 0, na.rm = TRUE),
                        n_p05_neg = sum(p < 0.05 & matched < 0, na.rm = TRUE),
                        expected_per_side = n * 0.025,
                        n_q10_pos = sum(q < 0.1 & matched > 0, na.rm = TRUE),
                        n_q10_neg = sum(q < 0.1 & matched < 0, na.rm = TRUE)
                )
        )
}
mm <- mismatch_test(MC, MET, "mc")
me <- mismatch_test(EX, MET, "exp")
print(rbind(mm$summary, me$summary))
# Reading: sd_matched >> sd_null and many drugs significant on BOTH sides means
# real heterogeneity. Significant drugs only on the positive side means a
# positive mean effect, not sign heterogeneity.

dens <- rbindlist(list(
        data.table(score = "mc", type = "matched", rho = mm$matched),
        data.table(score = "mc", type = "mismatched", rho = mm$null),
        data.table(score = "exp", type = "matched", rho = me$matched),
        data.table(score = "exp", type = "mismatched", rho = me$null)
))
p_mm <- ggplot(dens, aes(rho, colour = type)) +
        geom_density() + facet_wrap(~score) + geom_vline(xintercept = 0, linetype = 2) +
        theme_minimal() + labs(x = "per-drug rho (metabolite-centered)")
ggsave(file.path(out_dir, "matched_vs_mismatched_rho.pdf"), p_mm, width = 8, height = 4)

# Does the heterogeneity belong to MetChange or to the transcriptome-pool
# relationship? If sign flips for the same drugs under both scores, it is not
# MetChange-specific (both are production-side).
cat("\nSign concordance across drugs, mc vs baseline (Spearman of matched rho):\n")
print(cor.test(mm$matched, me$matched, method = "spearman", exact = FALSE))
cat("\nIs |rho| larger for mc than for the baseline? (paired Wilcoxon)\n")
print(wilcox.test(abs(mm$matched), abs(me$matched), paired = TRUE))


# C. IS THE PER-DRUG AUC BIMODALITY MORE THAN SMALL-SAMPLE NOISE? ##############
# Per-drug direction AUC uses only metabolites with an observed change, often
# 3-10 per class, so its null SD is large. Null: shuffle up/down labels among
# the changed metabolites of each drug (keeps counts per drug).
auc_null <- dtp[obs_label != 0, {
        pos <- obs_label > 0
        if (sum(pos) >= 3 && sum(!pos) >= 3)
                .(null_auc = replicate(200, auc_fn(mc_c, sample(pos))))
}, by = perturbation]
obs_auc <- auc_dir$wide$mc
cat("\nPer-drug direction AUC (mc): observed SD", round(sd(obs_auc), 3),
    "| null SD", round(sd(auc_null$null_auc), 3), "\n")
cat("Fraction |AUC - 0.5| > 0.1: observed", round(mean(abs(obs_auc - 0.5) > 0.1), 3),
    "| null", round(mean(abs(auc_null$null_auc - 0.5) > 0.1), 3), "\n")
p_auc_null <- ggplot() +
        geom_density(data = auc_null, aes(null_auc, colour = "null")) +
        geom_density(data = data.frame(a = obs_auc), aes(a, colour = "observed")) +
        theme_minimal() + labs(x = "per-drug direction AUC (mc)", colour = "")
ggsave(file.path(out_dir, "auc_observed_vs_null.pdf"), p_auc_null, width = 6, height = 4)


# D. SPLIT-HALF STABILITY OF THE SIGN (valid subgroup test) ####################
# Metabolites are split at random into halves A and B. Per-drug rho is computed
# in each half. If the sign is a real property of the drug:
#   * rho in A and rho in B correlate across drugs (reliability),
#   * drugs with negative mc rho in A are also negative in B (anti-correlation
#     reproduced on independent metabolites).
# The exp columns show whether the baseline follows the same drugs.
# Caveat: correlated metabolites can land in different halves, which inflates
# reliability a little. A stricter split puts whole pathways in one half.
split_reliability <- function(dt, B = 200) {
        dd <- copy(dt)
        mets <- unique(dd$metabolite)
        out <- t(replicate(B, {
                A <- sample(mets, floor(length(mets) / 2))
                dd[, grp := fifelse(metabolite %in% A, "A", "B")]
                r <- dd[, .(
                        mc = suppressWarnings(cor(met_c, mc_c, method = "spearman")),
                        exp = suppressWarnings(cor(met_c, exp_c, method = "spearman"))
                ), by = .(perturbation, grp)]
                w <- dcast(r, perturbation ~ grp, value.var = c("mc", "exp"))
                c(rel_mc = cor(w$mc_A, w$mc_B, method = "spearman", use = "complete.obs"),
                  rel_exp = cor(w$exp_A, w$exp_B, method = "spearman", use = "complete.obs"),
                  mcB_if_mcA_pos = mean(w$mc_B[w$mc_A > 0], na.rm = TRUE),
                  mcB_if_mcA_neg = mean(w$mc_B[w$mc_A < 0], na.rm = TRUE),
                  expB_if_mcA_pos = mean(w$exp_B[w$mc_A > 0], na.rm = TRUE),
                  expB_if_mcA_neg = mean(w$exp_B[w$mc_A < 0], na.rm = TRUE))
        }))
        # Intervals reflect the random split only, not sampling uncertainty over drugs.
        data.frame(stat = colnames(out), mean = colMeans(out, na.rm = TRUE),
                   lo = apply(out, 2, quantile, 0.025, na.rm = TRUE),
                   hi = apply(out, 2, quantile, 0.975, na.rm = TRUE), row.names = NULL)
}
cat("\nSplit-half reliability of the sign:\n")
print(split_reliability(dtp))























comp_info <- as.data.frame(data.table::fread(
        "/Users/guillem.santamaria/Documents/postdoc/comput/metabognn_dev/data/LINCS/compoundinfo_beta.txt"
))

drug_info <- unique(comp_info[, c("cmap_name", "moa")])
colnames(drug_info) <- c("perturbation", "moa")
dm_merged <- unique(do.call(rbind.data.frame, lapply(dat_list, function(x) x$dm[, c("compound", "schuhknecht_name")])))

match(drug_info$perturbation, dm_merged$compound)
drug_info$perturbation <- dm_merged$schuhknecht_name[match(drug_info$perturbation, dm_merged$compound)]
drug_info <- drug_info[!is.na(drug_info$perturbation), ]
drug_info$perturbation <- make.names(drug_info$perturbation)
drug_info <- drug_info[drug_info$moa != "", ]

# E. CAN THE SIGN BE PREDICTED FROM DRUG PROPERTIES? ###########################
# Provide `drug_info` with one row per drug: perturbation (names as in
# colnames(MET), i.e. make.names applied) plus whatever annotations you have:
#   targets_gem   TRUE if a DrugBank target is a GEM gene (metabolic enzyme)
#   moa_class     mechanism-of-action class
#   l1000_time, l1000_dose   conditions behind the averaged L1000 delta
#   sch_dose, sch_time       Schuhknecht condition (a mismatch could flip signs)
# Declare the hypothesis before looking (e.g. targets_gem) and validate on
# held-out drugs; testing many annotations on 193 drugs will find something.
if (exists("drug_info")) {
        rho_tab <- merge(
                data.table(perturbation = colnames(MET), rho_mc = mm$matched,
                           rho_exp = me$matched, q_mc = mm$q),
                as.data.table(drug_info), by = "perturbation"
        )
        if ("targets_gem" %in% names(rho_tab)) {
                print(wilcox.test(rho_mc ~ targets_gem, data = rho_tab))
                print(wilcox.test(rho_exp ~ targets_gem, data = rho_tab))
                print(summary(glm(I(rho_mc > 0) ~ targets_gem, data = rho_tab, family = binomial)))
        }
}

library(data.table)

rho_tab[, `:=`(
        delta_abs = abs(rho_mc) - abs(rho_exp),
        delta_signed = rho_mc - rho_exp
)]

rho_tab[, .(
        n = .N,
        mean_mc = mean(rho_mc, na.rm = TRUE),
        mean_exp = mean(rho_exp, na.rm = TRUE),
        mean_abs_mc = mean(abs(rho_mc), na.rm = TRUE),
        mean_abs_exp = mean(abs(rho_exp), na.rm = TRUE),
        mean_delta_abs = mean(delta_abs, na.rm = TRUE),
        median_delta_abs = median(delta_abs, na.rm = TRUE)
)]

moa_summary <- rho_tab[, .(
        n_drugs = uniqueN(perturbation),
        median_delta_abs = median(delta_abs, na.rm = TRUE),
        mean_delta_abs = mean(delta_abs, na.rm = TRUE)
), by = moa][order(-median_delta_abs)]

moa_summary


rho_drug <- unique(rho_tab[, .(
        perturbation, rho_mc, rho_exp, q_mc, delta_abs
)])

# Keep only sufficiently represented MoAs
moa_counts <- rho_tab[, .(
        n_drugs = uniqueN(perturbation)
), by = moa]

moa_keep <- moa_counts[n_drugs >= 5, moa]

moa_data <- unique(rho_tab[moa %in% moa_keep,
                           .(perturbation, moa, delta_abs)])

moa_data <- merge(
        moa_data,
        rho_drug[, .(perturbation, delta_abs)],
        by = "perturbation"
)



moa_data <- unique(
        rho_tab[moa %in% moa_keep, .(perturbation, moa, delta_abs)]
)

moa_data[, .(
        n_drugs = uniqueN(perturbation),
        median_delta = median(delta_abs),
        mean_delta = mean(delta_abs)
), by = moa][order(-median_delta)]

moa_summary <- moa_data[, .(
        n_drugs = uniqueN(perturbation),
        median_delta = median(delta_abs, na.rm = TRUE),
        mean_delta = mean(delta_abs, na.rm = TRUE)
), by = moa][order(-median_delta)]

moa_summary

# F. CONTROL: IS THE SPLIT-HALF STABILITY SPECIFIC TO THE MATCHED PAIRING? #####
# Reliability of per-drug rho between random halves of the metabolites, for
# matched pairs (score of drug i vs metabolomics of drug i) and for mismatched
# pairs (score of drug i vs metabolomics of a random other drug). If mismatched
# pairs are about as reliable, the stability comes from shared profile
# structure (correlated metabolites, common responses), not from a drug-specific
# relationship, and section D does not support the hypothesis.
split_rel_control <- function(S, MET, B = 200) {
        n <- ncol(S)
        mets <- seq_len(nrow(S))
        id <- seq_len(n)
        one <- function(rows, perm) vapply(id, function(i)
                suppressWarnings(cor(S[rows, i], MET[rows, perm[i]],
                                     method = "spearman", use = "pairwise.complete.obs")),
                numeric(1))
        res <- t(replicate(B, {
                A <- sample(mets, floor(length(mets) / 2))
                Bh <- setdiff(mets, A)
                perm <- sample(n)
                while (any(perm == id)) perm <- sample(n)   # no drug paired with itself
                mA <- one(A, id); mB <- one(Bh, id)         # matched, halves A and B
                xA <- one(A, perm); xB <- one(Bh, perm)     # mismatched, same halves
                c(rel_matched = cor(mA, mB, method = "spearman", use = "complete.obs"),
                  rel_mismatched = cor(xA, xB, method = "spearman", use = "complete.obs"),
                  B_if_A_neg_matched = mean(mB[mA < 0], na.rm = TRUE),
                  B_if_A_neg_mismatched = mean(xB[xA < 0], na.rm = TRUE),
                  B_if_A_pos_matched = mean(mB[mA > 0], na.rm = TRUE),
                  B_if_A_pos_mismatched = mean(xB[xA > 0], na.rm = TRUE))
        }))
        data.frame(stat = colnames(res), mean = colMeans(res, na.rm = TRUE),
                   lo = apply(res, 2, quantile, 0.025, na.rm = TRUE),
                   hi = apply(res, 2, quantile, 0.975, na.rm = TRUE), row.names = NULL)
}
cat("\nSplit-half control, mc:\n");  print(split_rel_control(MC, MET))
cat("\nSplit-half control, baseline:\n"); print(split_rel_control(EX, MET))

dd <- dtp[obs_label != 0]
dd[, `:=`(mc_r = frank(mc_c) / .N, exp_r = frank(exp_c) / .N,
          n_up = sum(obs_label > 0), n_dn = sum(obs_label < 0)), by = perturbation]
dd <- dd[n_up >= 8 & n_dn >= 8]
pooled_auc_ci(list(mc = dd$mc_r, exp = dd$exp_r), dd$obs_label > 0, dd$perturbation, B = 1000)


dd[, uniqueN(perturbation)]
