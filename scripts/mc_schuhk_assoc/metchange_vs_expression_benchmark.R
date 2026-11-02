################################################################################
# MetChange (ranks / rank_all / mix_script) vs producers-only expression
# baseline, benchmarked against Schuhknecht A549 metabolomics
################################################################################
# QUESTION
#   Do MetChange-derived metabolite deltas predict the drug-specific metabolite
#   changes measured by Schuhknecht better than a simple expression-derived
#   baseline (mean expression of the reactions that produce each metabolite)?
#
# CONVENTIONS
#   * Unit of inference = the drug. Pooled drug x metabolite rows are not
#     independent, so mixed models are descriptive only; all tests that decide
#     anything are per-drug.
#   * "_c" columns are centered per metabolite across drugs. That removes the
#     metabolite main effect and keeps the drug-specific part of each profile,
#     which is what the tool has to reproduce.
#   * MetChange scores are sign-flipped in get_pooled_deltas (high
#     inconsistency = low production), so positive = more production for all
#     three scores (mc, expression baseline, metabolomics log2FC).
#   * Only ranks, rank_all and mix_script are used. mixed/graded/comb0/comb1
#     had the duplicate-summing bug and are dropped.
#
# NOT TESTED against your data (written without R access): run section by
# section and watch for shape/name mismatches the first time.
################################################################################

suppressPackageStartupMessages({
        library(data.table)
        library(tidyr)
        library(tibble)
        library(dplyr)
        library(lme4)
        library(performance)
        library(ggplot2)
})
set.seed(1)

hdr <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")


# 0. CONFIG ####################################################################
base_dir <- Sys.getenv(
        "MC_BASE_DIR",
        "/Users/guillem.santamaria/Documents/postdoc/comput/metabognn_dev/results/mc_schuhk_assoc"
)
out_dir <- file.path(base_dir, "expression_baseline")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

metMatch_f <- file.path(base_dir, "metabolite_match_table_non_ambig_match.csv")
exp_f <- file.path(
        dirname(base_dir), "expr_met_baseline", "expr_met_products_delta_avg.csv"
)

# One folder per penalty; each holds the same four files.
pen_dirs <- c(
        ranks = "rankPenlty",          # rank penalty, ranking only GEM genes
        rank_all = "rankAllPenlty",    # rank penalty over all genes, then GEM filter
        mix_script = "mixScriptPenlty" # saturating mixture model, bug-free script
)
pen_files <- lapply(pen_dirs, function(d) list(
        mc = file.path(base_dir, d, "mc_delta_avg_matched.csv"),
        dm = file.path(base_dir, d, "drug_match_table.csv"),
        fc = file.path(base_dir, d, "sch_log2FC_matched.csv"),
        qv = file.path(base_dir, d, "sch_qvals_matched.csv")
))

primary <- "rank_all"   # penalty used for the single-penalty analyses (sections 7-12)
alpha_sig <- 0.05    # q-value cutoff defining an observed metabolomics change
perc_keep <- 0.75    # keep drugs in the top 25% by number of significant metabolites
n_perm_spec  <- 1000 # specificity permutation (section 7)
n_perm_calls <- 500  # label permutation for discretised calls (section 12)
B_boot <- 1000       # cluster bootstrap for pooled AUCs
match_status = "all"


# 1. FUNCTIONS #################################################################

## 1a. Alignment (logic unchanged from your script) ----------------------------
# Reads one penalty's files, divides each metabolite's row by its SD across
# drugs (so a few metabolites cannot drive the results), intersects metabolites
# across mc / metabolomics / baseline, matches drug names, and keeps the top
# (1 - perc) fraction of drugs by number of significant metabolites.
get_aligned_dfs <- function(exp, mc_f, dm_f, sch_fc_f, sch_qv_f,
                            perc = 0.75, sd_norm = TRUE, match_status = "all") {
        mc <- read.csv(mc_f, row.names = 1)
        dm <- read.csv(dm_f, row.names = 1)
        sch_fc <- read.csv(sch_fc_f, row.names = 1)
        sch_qv <- read.csv(sch_qv_f, row.names = 1)
        
        if (match_status != "all" &&
            match_status %in% unique(dm$match_status)) {
                
                keep_drugs <- dm$match_status == match_status
                dm <- dm[keep_drugs, ]
                mc <- mc[, keep_drugs]
                sch_fc <- sch_fc[, keep_drugs]
                sch_qv <- sch_qv[, keep_drugs]
        }
        
        if (sd_norm) {
                mc <- t(apply(mc, 1, function(x) x / sd(x)))
                exp <- exp[apply(exp, 1, sd) != 0, ]
                exp <- t(apply(exp, 1, function(x) x / sd(x)))
                sch_fc <- t(apply(sch_fc, 1, function(x) x / sd(x, na.rm = TRUE)))
        }
        comm_mets <- intersect(intersect(rownames(mc), rownames(sch_fc)), rownames(exp))
        mc <- mc[comm_mets, ]
        sch_fc <- sch_fc[comm_mets, ]
        sch_qv <- sch_qv[comm_mets, ]
        exp <- exp[comm_mets, ]
        
        exp <- exp[, match(make.names(dm$matched_column), colnames(exp))]
        colnames(exp) <- make.names(dm$schuhknecht_name)
        
        dm$n_sign_mets_sch <- apply(sch_qv, 2, function(x) sum(x <= 0.05 & !is.na(x)))
        quant <- quantile(dm$n_sign_mets_sch, perc)
        keep <- dm$n_sign_mets_sch >= quant
        list(
                expr = exp[, keep], mc = mc[, keep], sch_fc = sch_fc[, keep],
                sch_qv = sch_qv[, keep], dm = dm[keep, ]
        )
}

build_dat <- function(pen, perc = perc_keep, match_status = "all") {
        f <- pen_files[[pen]]
        get_aligned_dfs(expr, f$mc, f$dm, f$fc, f$qv, perc = perc,
                        match_status = match_status)
}

# Restrict all penalties to the same metabolites and drugs so that differences
# between penalties cannot come from different row/column sets.
harmonise <- function(dat_list) {
        mets <- Reduce(intersect, lapply(dat_list, function(d) rownames(d$mc)))
        drugs <- Reduce(intersect, lapply(dat_list, function(d) colnames(d$mc)))
        lapply(dat_list, function(d) {
                idx <- match(drugs, colnames(d$mc))  # dm rows are positionally aligned with columns
                d$dm <- d$dm[idx, ]
                d$expr <- d$expr[mets, drugs, drop = FALSE]
                d$mc <- d$mc[mets, drugs, drop = FALSE]
                d$sch_fc <- d$sch_fc[mets, drugs, drop = FALSE]
                d$sch_qv <- d$sch_qv[mets, drugs, drop = FALSE]
                d
        })
}

## 1b. Long format --------------------------------------------------------------
get_pooled_deltas <- function(dat) {
        to_long <- function(m, val) {
                as.data.frame(m) |>
                        rownames_to_column("metabolite") |>
                        pivot_longer(-metabolite, names_to = "perturbation", values_to = val)
        }
        pooled <- to_long(dat$expr, "expr_delta") |>
                inner_join(to_long(dat$mc, "mc_delta"), by = c("metabolite", "perturbation")) |>
                inner_join(to_long(dat$sch_fc, "met_delta"), by = c("metabolite", "perturbation"))
        pooled <- as.data.frame(pooled)
        pooled <- pooled[complete.cases(pooled[, c("expr_delta", "mc_delta", "met_delta")]), ]
        pooled$mc_delta <- pooled$mc_delta * -1  # inconsistency -> production
        pooled
}

# Long table used for every per-drug analysis: deltas, metabolite-centered
# deltas (_c), q-values and the observed -1/0/+1 label (the GNN's label format).
get_long <- function(dat) {
        pooled <- get_pooled_deltas(dat)
        q <- as.data.frame(dat$sch_qv) |>
                rownames_to_column("metabolite") |>
                pivot_longer(-metabolite, names_to = "perturbation", values_to = "qval")
        pooled <- left_join(pooled, q, by = c("metabolite", "perturbation"))
        dt <- as.data.table(pooled)
        dt[, `:=`(
                met_c = met_delta - mean(met_delta),
                mc_c = mc_delta - mean(mc_delta),
                exp_c = expr_delta - mean(expr_delta)
        ), by = metabolite]
        # NA q-value is treated as "not significant"
        dt[, obs_label := fifelse(!is.na(qval) & qval <= alpha_sig, sign(met_delta), 0)]
        dt
}

## 1c. Per-drug correlations ----------------------------------------------------
# For each drug: Spearman across metabolites between metabolomics and each score.
# rho = raw deltas; rhoc = metabolite-centered (drug-specific part only).
per_drug_stats <- function(d, min_met = 10) {
        d[, if (.N >= min_met) .(
                n_met = .N,
                rho_mc = suppressWarnings(cor(met_delta, mc_delta, method = "spearman")),
                rho_exp = suppressWarnings(cor(met_delta, expr_delta, method = "spearman")),
                rhoc_mc = suppressWarnings(cor(met_c, mc_c, method = "spearman")),
                rhoc_exp = suppressWarnings(cor(met_c, exp_c, method = "spearman"))
        ), by = perturbation]
}
get_pd <- function(dt) {
        pd <- per_drug_stats(dt)
        pd[, `:=`(delta = rho_mc - rho_exp, delta_c = rhoc_mc - rhoc_exp)]
        pd
}

## 1d. Paired-comparison helpers -------------------------------------------------
boot_mean_ci <- function(x, B = 5000) {
        x <- x[!is.na(x)]
        b <- replicate(B, mean(sample(x, replace = TRUE)))
        c(mean = mean(x), lo = unname(quantile(b, 0.025)), hi = unname(quantile(b, 0.975)))
}
# Paired sign-flip permutation test on the mean difference. This is the proper
# null for "A beats B": it keeps both scores intact. (Shuffling only mc and
# comparing with an intact expression baseline just destroys mc, so its null
# for the difference is trivially negative.)
signflip_p <- function(d, B = 10000) {
        d <- d[!is.na(d)]
        obs <- abs(mean(d))
        null <- replicate(B, abs(mean(d * sample(c(-1, 1), length(d), replace = TRUE))))
        (sum(null >= obs) + 1) / (B + 1)
}
# Compares two per-drug vectors a and b (a - b): mean difference with bootstrap
# CI, fraction of drugs where a > b, paired Wilcoxon and sign-flip p-values.
summarise_paired <- function(a, b, label) {
        d <- a - b
        ci <- boot_mean_ci(d)
        data.frame(
                comparison = label, n = sum(!is.na(d)),
                mean_diff = unname(ci["mean"]), ci_lo = unname(ci["lo"]), ci_hi = unname(ci["hi"]),
                win_frac = mean(d > 0, na.rm = TRUE),
                p_wilcox = suppressWarnings(wilcox.test(a, b, paired = TRUE)$p.value),
                p_signflip = signflip_p(d), row.names = NULL
        )
}

## 1e. AUC / ROC helpers ---------------------------------------------------------
# Rank-based AUC (Mann-Whitney), ties get average rank.
auc_fn <- function(score, pos) {
        ok <- !is.na(score) & !is.na(pos)
        score <- score[ok]; pos <- pos[ok]
        n1 <- sum(pos); n0 <- sum(!pos)
        if (n1 < 1 || n0 < 1) return(NA_real_)
        r <- rank(score)
        (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
roc_points <- function(score, pos) {
        ok <- !is.na(score) & !is.na(pos)
        score <- score[ok]; pos <- pos[ok]
        o <- order(score, decreasing = TRUE)
        pos <- pos[o]
        data.frame(fpr = c(0, cumsum(!pos) / sum(!pos)), tpr = c(0, cumsum(pos) / sum(pos)))
}
# Per-drug AUC.
#  direction: among metabolites with an observed change, does the score separate
#             up from down? (does the SIGN of the label come out right)
#  detection: does |score| separate changed from unchanged metabolites?
auc_by_drug <- function(dt, task, score_col, min_n = 3) {
        dt[, {
                s <- get(score_col)
                if (task == "direction") {
                        keep <- obs_label != 0
                        pos <- obs_label[keep] > 0
                        s <- s[keep]
                } else {
                        pos <- obs_label != 0
                        s <- abs(s)
                }
                if (sum(pos) >= min_n && sum(!pos) >= min_n)
                        .(auc = auc_fn(s, pos), n_pos = sum(pos), n_neg = sum(!pos))
        }, by = perturbation]
}
auc_summary <- function(dt, task, cols, min_n = 3) {
        wide <- NULL
        for (nm in names(cols)) {
                x <- auc_by_drug(dt, task, cols[[nm]], min_n)[, .(perturbation, auc)]
                setnames(x, "auc", nm)
                wide <- if (is.null(wide)) x else merge(wide, x, by = "perturbation")
        }
        per_score <- rbindlist(lapply(names(cols), function(nm) {
                ci <- boot_mean_ci(wide[[nm]])
                data.table(
                        task = task, score = nm, n_drugs = nrow(wide),
                        mean_auc = unname(ci["mean"]), lo = unname(ci["lo"]), hi = unname(ci["hi"]),
                        frac_above_half = mean(wide[[nm]] > 0.5),
                        p_vs_half = suppressWarnings(wilcox.test(wide[[nm]], mu = 0.5)$p.value)
                )
        }))
        paired <- rbindlist(lapply(
                list(c("mc", "exp"), c("mc_raw", "exp_raw"), c("comb", "mc"), c("comb", "exp")),
                function(p) as.data.table(summarise_paired(
                        wide[[p[1]]], wide[[p[2]]], paste(task, p[1], "-", p[2])
                ))
        ))
        list(wide = wide, per_score = per_score, paired = paired)
}
# Pooled AUC with a cluster bootstrap over drugs (drugs are the independent unit).
pooled_auc_ci <- function(S, pos, drug, B = 1000) {
        idx_by_drug <- split(seq_along(drug), drug)
        drugs <- names(idx_by_drug)
        est <- sapply(S, auc_fn, pos = pos)
        bs <- replicate(B, {
                idx <- unlist(idx_by_drug[sample(drugs, replace = TRUE)], use.names = FALSE)
                sapply(S, function(s) auc_fn(s[idx], pos[idx]))
        })
        d <- bs["mc", ] - bs["exp", ]
        list(
                auc = data.frame(
                        score = names(S), auc = unname(est),
                        lo = unname(apply(bs, 1, quantile, 0.025)),
                        hi = unname(apply(bs, 1, quantile, 0.975))
                ),
                diff_mc_exp = c(
                        est = unname(est["mc"] - est["exp"]),
                        lo = unname(quantile(d, 0.025)), hi = unname(quantile(d, 0.975))
                )
        )
}

## 1f. Discretised calls ---------------------------------------------------------
# The GNN labels are -1/0/+1. Here a score is discretised per drug (top frac of
# the score -> +1, bottom frac -> -1, rest 0) and compared with the observed label.
call_metrics <- function(obs, called) {
        idx <- !is.na(called) & called != 0
        s <- idx & obs != 0
        c(
                n_called = sum(idx),
                frac_sig = mean(obs[idx] != 0),               # called metabolites that truly changed
                sign_acc = mean(sign(obs[s]) == called[s]),   # of those, sign correct
                exact = mean(obs[idx] == called[idx])         # call equals observed label
        )
}
# Observed metrics vs a null where observed labels are shuffled WITHIN each drug.
# The same shuffles are applied to every score, so differences between scores
# (mc - exp) get a proper paired null.
perm_call_test <- function(dt, cols, B = 500) {
        stats_of <- function(obs) sapply(cols, function(cc) call_metrics(obs, dt[[cc]]))
        obs <- stats_of(dt$obs_label)
        null <- replicate(B, {
                dt[, obs_perm := sample(obs_label), by = perturbation]
                stats_of(dt$obs_perm)
        })  # array: metric x score x permutation
        dt[, obs_perm := NULL]
        rbindlist(lapply(c("frac_sig", "sign_acc", "exact"), function(m) {
                rbindlist(list(
                        data.table(
                                metric = m, what = names(cols), observed = obs[m, ],
                                null_mean = apply(null[m, , , drop = FALSE], 2, mean, na.rm = TRUE),
                                p_perm = sapply(seq_along(cols), function(j)
                                        (sum(null[m, j, ] >= obs[m, j], na.rm = TRUE) + 1) / (B + 1))
                        ),
                        data.table(
                                metric = m, what = "mc - exp",
                                observed = obs[m, "mc"] - obs[m, "exp"],
                                null_mean = mean(null[m, "mc", ] - null[m, "exp", ], na.rm = TRUE),
                                p_perm = (sum(abs(null[m, "mc", ] - null[m, "exp", ]) >=
                                                      abs(obs[m, "mc"] - obs[m, "exp"]), na.rm = TRUE) + 1) / (B + 1)
                        )
                ))
        }))
}


# 2. LOAD BASELINE #############################################################
# Producers-only baseline: per sample and metabolite, mean expression of the
# (gene-backed) reactions producing it, as control-subtracted deltas averaged
# per perturbation. Rows are GEM metabolite ids -> mapped to metabolomics names.
metMatch <- read.csv(metMatch_f, row.names = 1)
expr <- data.frame(data.table::fread(exp_f))
rownames(expr) <- expr$V1
expr <- expr[, 2:ncol(expr)]
rownames(expr) <- metMatch$name[match(rownames(expr), metMatch$id)]


# 3. ALIGN ALL PENALTIES #######################################################
hdr("3. Alignment")
dat_list <- lapply(names(pen_files), build_dat, match_status = match_status)
names(dat_list) <- names(pen_files)
dat_list <- harmonise(dat_list)
print(sapply(dat_list, function(d) c(n_met = nrow(d$mc), n_drug = ncol(d$mc))))

# Global (unstratified) Spearman between each score and metabolomics. Expect
# values near 0: global structure is dominated by metabolite main effects.
print(sapply(dat_list, function(d) c(
        mc = cor(as.vector(as.matrix(d$mc)), as.vector(as.matrix(d$sch_fc)),
                 method = "spearman", use = "complete.obs"),
        exp = cor(as.vector(as.matrix(d$expr)), as.vector(as.matrix(d$sch_fc)),
                  method = "spearman", use = "complete.obs")
)))

dt_list <- lapply(dat_list, get_long)
pd_list <- lapply(dt_list, get_pd)

# Primary table: adds a combined score (z(mc) + z(exp) within drug) used for the
# complementarity and ROC analyses.
zs <- function(x) as.numeric(scale(x))
dtp <- copy(dt_list[[primary]])
dtp[, comb_c := zs(mc_c) + zs(exp_c), by = perturbation]


# 4. POOLED MIXED MODELS (descriptive) #########################################
# What they check: whether the pooled slope of metabolomics on each score is
# positive, whether it is within-drug or between-drug, and how much variance it
# explains. p-values here are inflated (rows within drug/metabolite are not
# independent); the per-drug tests in section 5 are the ones to rely on.
fit_pooled <- function(pooled) {
        m_exp <- lmer(met_delta ~ expr_delta + (1 | perturbation) + (1 | metabolite), data = pooled)
        m_mc <- lmer(met_delta ~ mc_delta + (1 | perturbation) + (1 | metabolite), data = pooled)
        m_both <- lmer(met_delta ~ expr_delta + mc_delta + (1 | perturbation) + (1 | metabolite), data = pooled)
        m_exp_s <- lmer(met_delta ~ expr_delta + (1 + expr_delta || perturbation) + (1 | metabolite), data = pooled)
        m_mc_s <- lmer(met_delta ~ mc_delta + (1 + mc_delta || perturbation) + (1 | metabolite), data = pooled)
        
        slope_sd <- function(m, v) {
                vc <- as.data.frame(VarCorr(m))
                #vc$sdcor[vc$var1 == v & is.na(vc$var2)]
                vc$sdcor[which(vc$var1 == v & is.na(vc$var2))]
        }
        # Fraction of drugs with a positive slope implied by the random-slope model.
        # A bigger slope SD (and better AIC) does NOT mean more consistent: compare
        # these two numbers, not the AIC.
        pos_mc <- pnorm(fixef(m_mc_s)["mc_delta"] / slope_sd(m_mc_s, "mc_delta"))
        pos_exp <- pnorm(fixef(m_exp_s)["expr_delta"] / slope_sd(m_exp_s, "expr_delta"))
        
        a_mc_adds <- suppressMessages(anova(m_exp, m_both))   # does mc add to expression?
        a_exp_adds <- suppressMessages(anova(m_mc, m_both))   # does expression add to mc?
        
        wb <- function(v) {  # within- vs between-drug slope
                x <- pooled[[v]]
                pooled$x_mean <- ave(x, pooled$perturbation)
                pooled$x_within <- x - pooled$x_mean
                m <- lmer(met_delta ~ x_within + x_mean + (1 + x_within || perturbation) + (1 | metabolite), data = pooled)
                ci <- suppressMessages(confint(m, method = "Wald"))
                data.frame(
                        score = v, within = unname(fixef(m)["x_within"]),
                        within_lo = ci["x_within", 1], within_hi = ci["x_within", 2],
                        between = unname(fixef(m)["x_mean"]),
                        between_lo = ci["x_mean", 1], between_hi = ci["x_mean", 2]
                )
        }
        list(
                summary = data.frame(
                        slope_mc = unname(fixef(m_mc)["mc_delta"]), slope_exp = unname(fixef(m_exp)["expr_delta"]),
                        R2m_mc = as.numeric(r2_nakagawa(m_mc)$R2_marginal),
                        R2m_exp = as.numeric(r2_nakagawa(m_exp)$R2_marginal),
                        R2m_both = as.numeric(r2_nakagawa(m_both)$R2_marginal),
                        pos_frac_mc = unname(pos_mc), pos_frac_exp = unname(pos_exp),
                        dAIC_slope_mc_minus_exp = AIC(m_mc_s) - AIC(m_exp_s),
                        p_mc_adds_to_exp = a_mc_adds$`Pr(>Chisq)`[2],
                        p_exp_adds_to_mc = a_exp_adds$`Pr(>Chisq)`[2]
                ),
                wb = rbind(wb("mc_delta"), wb("expr_delta"))
        )
}
hdr("4. Pooled mixed models")
pooled_res <- lapply(dat_list, function(d) fit_pooled(get_pooled_deltas(d)))
print(rbindlist(lapply(names(pooled_res), function(p)
        data.table(penalty = p, pooled_res[[p]]$summary))))
print(rbindlist(lapply(names(pooled_res), function(p)
        data.table(penalty = p, pooled_res[[p]]$wb))))


# 5. PER-DRUG COMPARISONS ######################################################
# What they check: for each drug, is the within-drug metabolite profile of a
# score correlated with the measured profile, and does MetChange beat the
# baseline across drugs? Holm correction across the three penalties.
hdr("5. Per-drug comparisons")
common <- Reduce(intersect, lapply(pd_list, function(x) x$perturbation))
M <- sapply(pd_list, function(x) x$rhoc_mc[match(common, x$perturbation)])
base <- pd_list[[primary]]$rhoc_exp[match(common, pd_list[[primary]]$perturbation)]

# 5a. Each penalty vs the baseline (the central claim)
vs_base <- rbindlist(lapply(colnames(M), function(p)
        as.data.table(summarise_paired(M[, p], base, paste(p, "- baseline")))))
vs_base[, p_wilcox_holm := p.adjust(p_wilcox, "holm")]
vs_base[, p_signflip_holm := p.adjust(p_signflip, "holm")]
print(vs_base)

# 5b. Is each score above zero on its own?
print(data.table(
        score = c(colnames(M), "baseline"),
        mean_rhoc = c(colMeans(M, na.rm = TRUE), mean(base, na.rm = TRUE)),
        p_vs_zero = c(apply(M, 2, function(x) suppressWarnings(wilcox.test(x, mu = 0)$p.value)),
                      suppressWarnings(wilcox.test(base, mu = 0)$p.value))
))

# 5c. Penalties against each other (only the three valid ones)
Mc <- M[complete.cases(M), ]
print(friedman.test(Mc))
pen_pairs <- combn(colnames(M), 2, simplify = FALSE)
vs_pen <- rbindlist(lapply(pen_pairs, function(p)
        as.data.table(summarise_paired(M[, p[1]], M[, p[2]], paste(p[1], "-", p[2])))))
vs_pen[, p_wilcox_holm := p.adjust(p_wilcox, "holm")]
print(vs_pen)
cat("\nSpearman between penalties' per-drug rhoc (are they the same predictions?)\n")
print(round(cor(M, method = "spearman", use = "pairwise.complete.obs"), 3))


# 6. PENALTY SELECTION WITHOUT WINNER'S CURSE ##################################
# What it checks: if you choose the best penalty on half of the drugs, how much
# does it beat the baseline on the other half? Repeated random splits.
hdr("6. Split-half selection")
split_half <- function(M, base, B = 500) {
        n <- nrow(M)
        res <- replicate(B, {
                a <- sample(n, floor(n / 2)); b <- setdiff(seq_len(n), a)
                best <- which.max(colMeans(M[a, , drop = FALSE], na.rm = TRUE))
                c(best = unname(best), adv = mean(M[b, best] - base[b], na.rm = TRUE))
        })
        list(
                chosen = table(colnames(M)[res["best", ]]),
                heldout_adv = c(mean = mean(res["adv", ]),
                                lo = unname(quantile(res["adv", ], 0.025)),
                                hi = unname(quantile(res["adv", ], 0.975)),
                                frac_pos = mean(res["adv", ] > 0))
        )
}
print(split_half(M, base))


# 7. SPECIFICITY PERMUTATION (primary penalty) #################################
# What it checks: is the MetChange signal drug-specific? mc is shuffled across
# drugs within each metabolite (keeps each metabolite's distribution, breaks the
# drug-to-metabolomics link). It tests mc vs shuffled mc only; "mc vs expression"
# is tested by the paired sign-flip test in section 5.
hdr("7. Specificity permutation")
obs_mc <- mean(pd_list[[primary]]$rhoc_mc, na.rm = TRUE)
null_mc <- vapply(seq_len(n_perm_spec), function(i) {
        d <- copy(dt_list[[primary]])
        d[, mc_delta := sample(mc_delta), by = metabolite]
        d[, mc_c := mc_delta - mean(mc_delta), by = metabolite]
        s <- d[, if (.N >= 10) .(r = suppressWarnings(cor(met_c, mc_c, method = "spearman"))),
               by = perturbation]
        mean(s$r, na.rm = TRUE)
}, numeric(1))
print(c(obs = obs_mc, null_mean = mean(null_mc),
        p = (sum(null_mc >= obs_mc) + 1) / (n_perm_spec + 1)))
p_null <- ggplot(data.frame(null = null_mc), aes(null)) +
        geom_histogram(bins = 40) + geom_vline(xintercept = obs_mc, colour = "red") +
        theme_minimal() + ggtitle("Specificity null: mean per-drug rho (mc)")
ggsave(file.path(out_dir, "specificity_null.pdf"), p_null, width = 6, height = 4)


# 8. DRUG-SELECTION CUTOFF AND RESPONSE STRENGTH ###############################
# What it checks: (a) whether the 0.75 cutoff was lucky, (b) whether the
# MetChange advantage grows with the number of significant metabolites.
hdr("8. perc sweep and n_sign dependence")
sweep <- rbindlist(lapply(c(0, 0.5, 0.75, 0.9), function(pc) {
        d <- build_dat(primary, pc)
        pd <- get_pd(get_long(d))
        s <- as.data.table(summarise_paired(pd$rhoc_mc, pd$rhoc_exp, paste0("perc=", pc)))
        s[, `:=`(mean_rhoc_mc = mean(pd$rhoc_mc, na.rm = TRUE),
                 mean_rhoc_exp = mean(pd$rhoc_exp, na.rm = TRUE))]
        s
}))
print(sweep)

d0 <- build_dat(primary, 0)
pd0 <- get_pd(get_long(d0))
pd0[, n_sign := d0$dm$n_sign_mets_sch[match(perturbation, colnames(d0$mc))]]
print(cor.test(pd0$delta_c, pd0$n_sign, method = "spearman", exact = FALSE))
pd0[, quartile := ntile(n_sign, 4)]
print(pd0[, .(n = .N, median_nsign = median(n_sign),
              rhoc_mc = mean(rhoc_mc, na.rm = TRUE),
              rhoc_exp = mean(rhoc_exp, na.rm = TRUE),
              delta_c = mean(delta_c, na.rm = TRUE)), by = quartile][order(quartile)])
p_ns <- ggplot(pd0, aes(n_sign, delta_c)) +
        geom_point(alpha = 0.5) + geom_smooth(method = "loess", formula = y ~ x) +
        geom_hline(yintercept = 0, linetype = 2) + theme_minimal() +
        labs(x = "n significant metabolites (drug)", y = "rhoc(mc) - rhoc(baseline)")
ggsave(file.path(out_dir, "delta_vs_nsign.pdf"), p_ns, width = 6, height = 4)


# 9. COMPLEMENTARITY AT THE DRUG LEVEL #########################################
# What it checks: is MetChange just a rescaled version of the expression
# baseline (high rho_mc_vs_exp), and does z(mc)+z(exp) beat each alone?
hdr("9. Complementarity")
comp <- dtp[, if (.N >= 10) .(
        rho_mc = suppressWarnings(cor(met_c, mc_c, method = "spearman")),
        rho_exp = suppressWarnings(cor(met_c, exp_c, method = "spearman")),
        rho_comb = suppressWarnings(cor(met_c, comb_c, method = "spearman")),
        rho_mc_vs_exp = suppressWarnings(cor(mc_c, exp_c, method = "spearman"))
), by = perturbation]
cat("mean within-drug correlation between mc and baseline across metabolites:",
    round(mean(comp$rho_mc_vs_exp, na.rm = TRUE), 3), "\n")
print(rbind(
        summarise_paired(comp$rho_comb, comp$rho_mc, "comb - mc"),
        summarise_paired(comp$rho_comb, comp$rho_exp, "comb - exp")
))


# 10. STRATIFICATION BY METABOLITE CONNECTIVITY (optional) #####################
# What it checks: does MetChange win where network propagation should matter
# (highly connected metabolites) or also where plain GPR logic would suffice?
# Needs metabolite_degree.csv exported from Python:
#   deg = pd.DataFrame({"id": [m.id for m in gem.metabolites],
#     "n_prod": [len(met2rxnsProd.get(m.id, [])) for m in gem.metabolites],
#     "n_reac": [len(met2rxnsReac.get(m.id, [])) for m in gem.metabolites]})
#   deg.to_csv(".../mc_schuhk_assoc/metabolite_degree.csv", index=False)
hdr("10. Connectivity strata")
deg_f <- file.path(base_dir, "metabolite_degree.csv")
if (file.exists(deg_f)) {
        deg <- read.csv(deg_f)
        deg$metabolite <- metMatch$name[match(deg$id, metMatch$id)]
        deg <- deg |>
                filter(!is.na(metabolite)) |>
                group_by(metabolite) |>
                summarise(degree = sum(n_prod + n_reac), .groups = "drop")
        deg$stratum <- factor(ntile(deg$degree, 3), labels = c("low", "mid", "high"))
        dtd <- merge(dtp, as.data.table(deg), by = "metabolite")
        strat <- dtd[, if (.N >= 10) .(
                n_met = .N,
                rhoc_mc = suppressWarnings(cor(met_c, mc_c, method = "spearman")),
                rhoc_exp = suppressWarnings(cor(met_c, exp_c, method = "spearman"))
        ), by = .(stratum, perturbation)]
        strat[, delta_c := rhoc_mc - rhoc_exp]
        print(strat[, c(.(n_drugs = .N, mean_rhoc_mc = mean(rhoc_mc, na.rm = TRUE),
                          mean_rhoc_exp = mean(rhoc_exp, na.rm = TRUE)),
                        as.list(boot_mean_ci(delta_c))), by = stratum][order(stratum)])
} else {
        cat("metabolite_degree.csv not found; skipping.\n")
}


# 11. AUROC AND ROC CURVES #####################################################
# What it checks, in the format the GNN will use (observed label = q <= 0.05 and
# sign of log2FC):
#   direction: among metabolites that truly changed, does the score rank the ups
#              above the downs?
#   detection: does |score| rank truly changed metabolites above unchanged?
# Primary inference = per-drug AUC (mean vs 0.5, and mc vs baseline paired).
# Pooled ROC is for the figure; its CI is a cluster bootstrap over drugs.
hdr("11. AUROC")
score_cols <- c(mc = "mc_c", exp = "exp_c", comb = "comb_c",
                mc_raw = "mc_delta", exp_raw = "expr_delta")  # *_raw = sensitivity
auc_dir <- auc_summary(dtp, "direction", score_cols)
auc_det <- auc_summary(dtp, "detection", score_cols)
print(rbind(auc_dir$per_score, auc_det$per_score))
print(rbind(auc_dir$paired, auc_det$paired))

main_scores <- c(mc = "mc_c", exp = "exp_c", comb = "comb_c")
pooled_out <- lapply(c("direction", "detection"), function(task) {
        d <- if (task == "direction") dtp[obs_label != 0] else dtp
        pos <- if (task == "direction") d$obs_label > 0 else d$obs_label != 0
        tf <- if (task == "direction") identity else abs
        S <- lapply(main_scores, function(cn) tf(d[[cn]]))
        ci <- pooled_auc_ci(S, pos, d$perturbation, B = B_boot)
        roc <- rbindlist(lapply(names(S), function(s) {
                p <- roc_points(S[[s]], pos); p$score <- s; p
        }))
        roc$task <- task
        ci$auc$task <- task
        list(roc = roc, auc = ci$auc, diff = ci$diff_mc_exp)
})
names(pooled_out) <- c("direction", "detection")
print(rbindlist(lapply(pooled_out, `[[`, "auc")))
print(lapply(pooled_out, `[[`, "diff"))

roc_df <- rbindlist(lapply(pooled_out, `[[`, "roc"))
auc_lab <- rbindlist(lapply(pooled_out, `[[`, "auc"))
auc_lab[, `:=`(x = 0.6, y = 0.25 - 0.07 * (seq_len(.N) - 1),
               lab = sprintf("%s: AUC %.3f [%.3f-%.3f]", score, auc, lo, hi)), by = task]
p_roc <- ggplot(roc_df, aes(fpr, tpr, colour = score)) +
        geom_line() + geom_abline(linetype = 2) +
        geom_text(data = auc_lab, aes(x, y, label = lab, colour = score),
                  hjust = 0.5, size = 3, show.legend = FALSE) +
        facet_wrap(~task) + coord_equal() + theme_minimal() +
        labs(x = "False positive rate", y = "True positive rate")
ggsave(file.path(out_dir, "roc_pooled.pdf"), p_roc, width = 9, height = 4.8)

# Per-drug AUC distributions (the primary evidence)
auc_long <- rbindlist(list(
        melt(auc_dir$wide, id.vars = "perturbation", variable.name = "score", value.name = "auc")[, task := "direction"],
        melt(auc_det$wide, id.vars = "perturbation", variable.name = "score", value.name = "auc")[, task := "detection"]
))[score %in% names(main_scores)]
p_auc <- ggplot(auc_long, aes(score, auc)) +
        geom_violin() + geom_boxplot(width = 0.15, outlier.shape = NA) +
        geom_hline(yintercept = 0.5, linetype = 2) + facet_wrap(~task) + theme_minimal() +
        labs(y = "per-drug AUC")
ggsave(file.path(out_dir, "auc_per_drug.pdf"), p_auc, width = 8, height = 4.5)


x <- auc_dir$wide$mc

summary(x)
median(x, na.rm = TRUE)
mean(x, na.rm = TRUE)

plot(density(x, na.rm = TRUE))
abline(v = median(x, na.rm = TRUE), col = "red", lwd = 2)
abline(v = mean(x, na.rm = TRUE), col = "blue", lwd = 2)


plot_roc_drug <- function(drug, task = c("direction", "detection"), dt = dtp) {
        task <- match.arg(task)
        
        d <- copy(dt[perturbation == drug])
        
        if (task == "direction") {
                # Among significantly changed metabolites: increases vs decreases
                d <- d[obs_label != 0]
                pos <- d$obs_label > 0
                scores <- list(
                        Expression = d$exp_c,
                        MetChange = d$mc_c,
                        Combined = d$comb_c
                )
        } else {
                # Among all metabolites: changed vs unchanged
                pos <- d$obs_label != 0
                scores <- list(
                        Expression = abs(d$exp_c),
                        MetChange = abs(d$mc_c),
                        Combined = abs(d$comb_c)
                )
        }
        
        if (sum(pos) < 3 || sum(!pos) < 3) {
                stop("Not enough positive or negative observations for this drug and task.")
        }
        
        roc <- rbindlist(lapply(names(scores), function(nm) {
                s <- scores[[nm]]
                p <- roc_points(s, pos)
                p$score <- sprintf("%s (AUC = %.2f)", nm, auc_fn(s, pos))
                p
        }))
        
        ggplot(roc, aes(x = fpr, y = tpr, colour = score, group = score)) +
                geom_abline(
                        slope = 1, intercept = 0,
                        linetype = 2, colour = "grey60"
                ) +
                geom_path(linewidth = 1) +
                coord_equal() +
                labs(
                        title = paste(drug, "-", task),
                        x = "False positive rate",
                        y = "True positive rate",
                        colour = NULL
                ) +
                theme_bw() +
                theme(legend.position = "bottom")
}

# 12. DISCRETISED CALLS (-1 / 0 / +1) ##########################################
# What it checks: if you turned each score into GNN-style labels by calling the
# top/bottom `frac` of metabolites per drug, how often is a called metabolite
# truly changed (frac_sig), how often is the sign right among the changed
# (sign_acc), and does the call equal the observed label (exact)? Compared with
# shuffled observed labels within drug; "mc - exp" is the paired difference.
hdr("12. Discretised calls")
calls_out <- rbindlist(lapply(c(0.05, 0.10, 0.20), function(frac) {
        dd <- copy(dtp)
        cn <- character()
        for (nm in names(main_scores)) {
                cc <- paste0("call_", nm)
                sc <- main_scores[[nm]]
                dd[, (cc) := {
                        s <- get(sc)
                        hi <- quantile(s, 1 - frac, na.rm = TRUE, names = FALSE)
                        lo <- quantile(s, frac, na.rm = TRUE, names = FALSE)
                        fifelse(s >= hi, 1, fifelse(s <= lo, -1, 0))
                }, by = perturbation]
                cn[nm] <- cc
        }
        res <- perm_call_test(dd, cn, B = n_perm_calls)
        res[, frac := frac]
        res
}))
print(calls_out[order(frac, metric, what)])


# 13. SAVE ######################################################################
saveRDS(list(
        pooled = pooled_res, vs_base = vs_base, vs_pen = vs_pen, sweep = sweep,
        comp = comp, auc_dir = auc_dir, auc_det = auc_det, pooled_auc = pooled_out,
        calls = calls_out, pd_list = pd_list
), file.path(out_dir, "benchmark_results.rds"))
fwrite(vs_base, file.path(out_dir, "vs_baseline.csv"))
fwrite(calls_out, file.path(out_dir, "discretised_calls.csv"))
cat("\nDone. Outputs in", out_dir, "\n")
