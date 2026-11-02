################################################################################
# Additional analyses: timing, pool balance, mechanism, pathway level, retrieval
################################################################################
# Run AFTER metchange_vs_expression_benchmark.R and sign_heterogeneity_analysis.R
# in the same R session. Needs from them: dtp, dt_list, build_dat, get_long,
# get_pd, get_aligned_dfs, harmonise, summarise_paired, boot_mean_ci, signflip_p,
# auc_fn, pen_files, primary, perc_keep, alpha_sig, metMatch, base_dir, out_dir,
# MET, MC, EX, mismatch_test, to_wide.
# Needs drug_mechanism_annotation.csv from annotate_drugs_metabolites.py.
#
# A. TIMING      Does concordance fall with the time gap between the LINCS profile
#                (6 / 24 / 48 h) and the Schuhknecht measurement (48 h)?
# B. POOL BALANCE A pool changes with production MINUS consumption. Does a
#                consumer-side or a net (producers - consumers) expression score
#                explain pools better than the producers-only baseline, and does
#                adding consumption to MetChange help? A total (producers +
#                consumers) score tests the alternative: pools follow overall
#                enzyme expression around the metabolite, not flux direction.
# C. MECHANISM   Pre-declared: (H1) drugs acting directly on a GEM enzyme or
#                transporter have LOWER concordance (a block raises the substrate
#                pool while lowering flux); (H2) drugs whose target reaches the
#                transcriptional route have HIGHER concordance.
# D. PATHWAYS    Concordance after averaging metabolites per GEM subsystem, with
#                the mismatched-pairs null.
# E. RETRIEVAL   The tool's use case: given a measured metabolomics profile, is the
#                drug that produced it ranked above the other drugs?
# F. DETECTION   Why is the detection AUC below 0.5?
#
# Unit of inference is the drug throughout. All drugs (perc = 0) are used where
# the analysis adjusts for response strength (log number of significant
# metabolites); the top-25% responders (dtp) where the earlier sections did.
################################################################################
stopifnot(exists("dtp"), exists("MET"), exists("mismatch_test"), exists("build_dat"))
set.seed(21)
add_dir <- file.path(out_dir, "additional")
dir.create(add_dir, showWarnings = FALSE)

# Per-drug metabolite-centered Spearman, for any pair of score columns
per_drug_rho <- function(dt, score_cols, min_met = 10) {
        dt[, if (.N >= min_met) lapply(score_cols, function(cc)
                suppressWarnings(cor(met_c, get(cc), method = "spearman"))), by = perturbation]
}
ci_row <- function(x, label) {
        ci <- boot_mean_ci(x)
        data.table(what = label, n = sum(!is.na(x)), mean = ci[["mean"]], lo = ci[["lo"]], hi = ci[["hi"]],
                   p_vs_zero = suppressWarnings(wilcox.test(x, mu = 0)$p.value))
}

# All drugs, with their condition match and response strength
d_all <- build_dat(primary, 0)
dt_all <- get_long(d_all)
pd_all <- get_pd(dt_all)
dm_all <- as.data.table(d_all$dm)
dm_all[, `:=`(perturbation = colnames(d_all$mc), n_sign = n_sign_mets_sch)]
pd_all <- merge(pd_all, dm_all[, .(perturbation, compound, time_distance_h, match_status,
                                   validated_metabolic, n_sign)], by = "perturbation")


# A. TIMING ####################################################################
hdr("A. Concordance by LINCS-Schuhknecht time gap")
timing <- rbindlist(lapply(sort(unique(pd_all$time_distance_h)), function(t) {
        x <- pd_all[time_distance_h == t]
        rbind(cbind(ci_row(x$rhoc_mc, "mc"), gap_h = t, median_nsign = median(x$n_sign)),
              cbind(ci_row(x$rhoc_exp, "exp"), gap_h = t, median_nsign = median(x$n_sign)))
}))
print(timing)
# Adjusted for response strength: does the gap still matter?
for (v in c("rhoc_mc", "rhoc_exp")) {
        cat("\n", v, "~ time gap + log1p(n_sign), all drugs\n")
        print(round(coef(summary(lm(reformulate(c("factor(time_distance_h)", "log1p(n_sign)"), v),
                                    data = pd_all))), 4))
}
# Strong responders only (top 25% by n_sign within all drugs)
strong <- pd_all[n_sign >= quantile(n_sign, perc_keep)]
cat("\nStrong responders by gap (n, mean rhoc_mc, mean rhoc_exp):\n")
print(strong[, .(n = .N, rhoc_mc = mean(rhoc_mc, na.rm = TRUE), rhoc_exp = mean(rhoc_exp, na.rm = TRUE)),
             by = time_distance_h][order(time_distance_h)])
fwrite(timing, file.path(add_dir, "A_timing.csv"))


# B. POOL BALANCE ##############################################################
hdr("B. Producers, consumers and net expression")
read_baseline <- function(f) {
        x <- data.frame(data.table::fread(f))
        rownames(x) <- x$V1
        x <- x[, -1]
        rownames(x) <- metMatch$name[match(rownames(x), metMatch$id)]
        x[!is.na(rownames(x)), ]
}
prod <- read_baseline(file.path(dirname(base_dir), "expr_met_baseline", "expr_met_products_delta_avg.csv"))
reac <- read_baseline(file.path(dirname(base_dir), "expr_met_baseline", "expr_met_reactants_delta_avg.csv"))
rows <- intersect(rownames(prod), rownames(reac))
cols <- intersect(colnames(prod), colnames(reac))
net <- prod[rows, cols] - reac[rows, cols]   # more production than consumption -> pool should rise
total <- prod[rows, cols] + reac[rows, cols] # all enzymes around the metabolite

f <- pen_files[[primary]]
variants <- list(prod = prod, reac = reac, net = net, total = total)
dat_v <- harmonise(lapply(variants, function(e) get_aligned_dfs(e, f$mc, f$dm, f$fc, f$qv, perc = perc_keep)))
long_v <- lapply(dat_v, get_long)
pb <- long_v$prod[, .(metabolite, perturbation, met_c, mc_c, prod_c = exp_c)]
pb <- merge(pb, long_v$reac[, .(metabolite, perturbation, reac_c = exp_c)], by = c("metabolite", "perturbation"))
pb <- merge(pb, long_v$net[, .(metabolite, perturbation, net_c = exp_c)], by = c("metabolite", "perturbation"))
pb <- merge(pb, long_v$total[, .(metabolite, perturbation, total_c = exp_c)], by = c("metabolite", "perturbation"))
zs <- function(x) as.numeric(scale(x))
pb[, `:=`(neg_reac_c = -reac_c, mc_minus_reac = zs(mc_c) - zs(reac_c)), by = perturbation]
rb <- per_drug_rho(pb, c(mc = "mc_c", prod = "prod_c", neg_reac = "neg_reac_c", net = "net_c",
                         total = "total_c", mc_minus_reac = "mc_minus_reac"))
cat(nrow(rb), "drugs,", uniqueN(pb$metabolite), "metabolites (top-25% responders)\n")
print(rbindlist(lapply(c("mc", "prod", "neg_reac", "net", "total", "mc_minus_reac"), function(v) ci_row(rb[[v]], v))))
print(rbind(
        summarise_paired(rb$net, rb$prod, "net - producers"),
        summarise_paired(rb$neg_reac, rb$prod, "-consumers - producers"),
        summarise_paired(rb$total, rb$prod, "total - producers"),
        summarise_paired(rb$mc, rb$total, "mc - total"),
        summarise_paired(rb$mc, rb$net, "mc - net"),
        summarise_paired(rb$mc_minus_reac, rb$mc, "(mc - consumers) - mc")
))
# What does MetChange track? Within-drug correlation of mc with each expression side
cat("\nWithin-drug Spearman of mc with producer, consumer and total expression (mean over drugs):\n")
print(pb[, .(prod = cor(mc_c, prod_c, method = "spearman"), reac = cor(mc_c, reac_c, method = "spearman"),
             total = cor(mc_c, total_c, method = "spearman")), by = perturbation][
        , lapply(.SD, mean, na.rm = TRUE), .SDcols = c("prod", "reac", "total")])
fwrite(rb, file.path(add_dir, "B_pool_balance_per_drug.csv"))


# C. MECHANISM #################################################################
hdr("C. Mechanism strata (pre-declared H1, H2)")
ann <- fread(file.path(base_dir, "drug_mechanism_annotation.csv"))
ann <- ann[, .(compound, targets_gem, direct_gem_inhibitor, on_route,
               mech_targets_gem, mech_direct_gem_inhibitor, mech_on_route)]
pm <- merge(pd_all, ann, by = "compound", all.x = TRUE)
mech_test <- function(var, hypothesis, data = pm) {
        d <- data[!is.na(get(var))]
        d[, g := as.logical(get(var))]
        rbindlist(lapply(c("rhoc_mc", "rhoc_exp"), function(s) {
                m <- lm(reformulate(c("g", "log1p(n_sign)"), s), data = d)
                co <- coef(summary(m))["gTRUE", ]
                ci <- confint(m)["gTRUE", ]
                st <- d[n_sign >= quantile(pd_all$n_sign, perc_keep)]
                data.table(variable = var, hypothesis = hypothesis, score = s,
                           n_true = sum(d$g), n_false = sum(!d$g),
                           diff_adj = co[["Estimate"]], lo = ci[[1]], hi = ci[[2]], p_adj = co[["Pr(>|t|)"]],
                           strong_n_true = sum(st$g), strong_n_false = sum(!st$g),
                           strong_mean_true = mean(st[g == TRUE, get(s)], na.rm = TRUE),
                           strong_mean_false = mean(st[g == FALSE, get(s)], na.rm = TRUE),
                           strong_p_wilcox = tryCatch(suppressWarnings(
                                   wilcox.test(st[g == TRUE, get(s)], st[g == FALSE, get(s)])$p.value),
                                   error = function(e) NA_real_))
        }))
}
mech <- rbindlist(list(
        mech_test("mech_targets_gem", "H1: lower"),
        mech_test("mech_on_route", "H2: higher"),
        mech_test("targets_gem", "H1 (all targets): lower"),
        mech_test("on_route", "H2 (all targets): higher"),
        mech_test("validated_metabolic", "exploratory")
))
# Holm over the two pre-declared primary tests on MetChange
prim <- mech$variable %in% c("mech_targets_gem", "mech_on_route") & mech$score == "rhoc_mc"
mech[prim, p_adj_holm := p.adjust(p_adj, "holm")]
print(mech)
fwrite(mech, file.path(add_dir, "C_mechanism.csv"))


# D. PATHWAY LEVEL #############################################################
hdr("D. Subsystem-level concordance")
generic <- c("Transport reactions", "Pool reactions", "Exchange/demand reactions", "Miscellaneous",
             "Isolated", "Artificial reactions", "")
met_sub <- as.data.table(metMatch)[, .(metabolite = name, subsystem = subsystems)]
met_sub <- met_sub[, .(subsystem = unlist(strsplit(subsystem, "|", fixed = TRUE))), by = metabolite]
met_sub <- unique(met_sub[!subsystem %in% generic])
pathway_long <- function(dt, min_mets = 4) {
        x <- merge(dt[, .(metabolite, perturbation, met_c, mc_c, exp_c)], met_sub, by = "metabolite",
                   allow.cartesian = TRUE)
        x <- x[, .(n = .N, met = mean(met_c), mc = mean(mc_c), exp = mean(exp_c)), by = .(subsystem, perturbation)]
        x <- x[n >= min_mets]
        x[, `:=`(met_c = met - mean(met), mc_c = mc - mean(mc), exp_c = exp - mean(exp)), by = subsystem]
        x
}
pw_strong <- pathway_long(dtp)
pw_all <- pathway_long(dt_all)
cat(uniqueN(pw_strong$subsystem), "subsystems with >= 4 measured metabolites\n")
pw_wide <- function(x, v) {
        w <- dcast(x, subsystem ~ perturbation, value.var = v)
        m <- as.matrix(w[, -1, with = FALSE]); rownames(m) <- w$subsystem; m
}
pw_res <- rbindlist(lapply(list(strong = pw_strong, all = pw_all), function(x) {
        PM <- pw_wide(x, "met_c")
        rbind(mismatch_test(pw_wide(x, "mc_c"), PM, "mc")$summary,
              mismatch_test(pw_wide(x, "exp_c"), PM, "exp")$summary)
}), idcol = "drugs")
print(pw_res)
rp <- per_drug_rho(pw_strong[, .(perturbation, met_c, mc_c, exp_c)], c(mc = "mc_c", exp = "exp_c"), min_met = 8)
print(summarise_paired(rp$mc, rp$exp, "pathway level, strong responders: mc - exp"))
fwrite(pw_res, file.path(add_dir, "D_pathway_level.csv"))


# E. RETRIEVAL #################################################################
hdr("E. Retrieval: is the right drug ranked first for its measured profile?")
# R[i, j] = similarity of drug i's predicted profile to drug j's measured profile.
# For each measured profile j (the query), the percentile of the true drug j among
# all predictions (1 = top, 0.5 = chance). Null: shuffle which prediction belongs
# to which drug. The same shuffles are used for mc and exp (paired null).
retrieval <- function(S, M, B = 1000) {
        stopifnot(identical(colnames(S), colnames(M)))
        R <- cor(S, M, method = "spearman", use = "pairwise.complete.obs")
        pct_of <- function(R, truth) vapply(seq_len(ncol(R)), function(j)
                mean(R[-truth[j], j] < R[truth[j], j], na.rm = TRUE), numeric(1))
        n <- ncol(R)
        obs <- pct_of(R, seq_len(n))
        null <- replicate(B, mean(pct_of(R, sample(n))))
        list(pct = obs, summary = data.table(
                n_drugs = n, mean_percentile = mean(obs), median_percentile = median(obs),
                top10pct = mean(obs >= 0.9), top5pct = mean(obs >= 0.95),
                null_mean = mean(null), p_perm = (sum(null >= mean(obs)) + 1) / (B + 1)))
}
ret <- list()
wide_all <- list(MET = to_wide(dt_all, "met_c"), MC = to_wide(dt_all, "mc_c"), EX = to_wide(dt_all, "exp_c"))
sets <- list(
        `metabolites, strong` = list(M = MET, mc = MC, exp = EX),
        `metabolites, all` = list(M = wide_all$MET, mc = wide_all$MC, exp = wide_all$EX),
        `subsystems, strong` = list(M = pw_wide(pw_strong, "met_c"), mc = pw_wide(pw_strong, "mc_c"),
                                    exp = pw_wide(pw_strong, "exp_c"))
)
ret_tab <- rbindlist(lapply(names(sets), function(nm) {
        s <- sets[[nm]]
        r_mc <- retrieval(s$mc, s$M); r_ex <- retrieval(s$exp, s$M)
        ret[[nm]] <<- list(mc = r_mc$pct, exp = r_ex$pct)
        rbind(cbind(set = nm, score = "mc", r_mc$summary), cbind(set = nm, score = "exp", r_ex$summary))
}))
print(ret_tab)
print(rbindlist(lapply(names(ret), function(nm)
        as.data.table(summarise_paired(ret[[nm]]$mc, ret[[nm]]$exp, paste(nm, ": mc - exp percentile"))))))
fwrite(ret_tab, file.path(add_dir, "E_retrieval.csv"))
p_ret <- ggplot(rbindlist(lapply(names(ret), function(nm) rbind(
        data.table(set = nm, score = "mc", pct = ret[[nm]]$mc),
        data.table(set = nm, score = "exp", pct = ret[[nm]]$exp)))), aes(pct, colour = score)) +
        geom_density() + geom_vline(xintercept = 0.5, linetype = 2) + facet_wrap(~set) + theme_minimal() +
        labs(x = "percentile of the true drug among all predictions (1 = ranked first)")
ggsave(file.path(add_dir, "E_retrieval.pdf"), p_ret, width = 10, height = 3.8)


# F. DETECTION BELOW CHANCE ####################################################
hdr("F. Why |score| is lower for metabolites that change")
# If metabolites that change often in the metabolomics are the ones MetChange
# rarely moves (or the opposite), |score| ranks them below unchanged ones.
met_level <- dtp[, .(frac_changed = mean(obs_label != 0), mean_abs_mc = mean(abs(mc_c)),
                     mean_abs_exp = mean(abs(exp_c)), sd_met = sd(met_delta)), by = metabolite]
print(cor.test(met_level$frac_changed, met_level$mean_abs_mc, method = "spearman", exact = FALSE))
print(cor.test(met_level$frac_changed, met_level$mean_abs_exp, method = "spearman", exact = FALSE))
# Detection AUC within drug after ranking each metabolite against itself across drugs
# (removes the metabolite-level confound): is a metabolite's |score| higher in the
# drugs where it changed?
dtp[, `:=`(mc_rank = frank(abs(mc_c)) / .N, exp_rank = frank(abs(exp_c)) / .N), by = metabolite]
det_w <- dtp[, if (sum(obs_label != 0) >= 3 && sum(obs_label == 0) >= 3) .(
        mc = auc_fn(mc_rank, obs_label != 0), exp = auc_fn(exp_rank, obs_label != 0)), by = perturbation]
print(rbind(ci_row(det_w$mc, "detection AUC, metabolite-ranked |mc|"),
            ci_row(det_w$exp, "detection AUC, metabolite-ranked |exp|")))
cat("(p_vs_zero above tests AUC against 0, not 0.5; vs 0.5:",
    signif(wilcox.test(det_w$mc, mu = 0.5)$p.value, 3), "mc,",
    signif(wilcox.test(det_w$exp, mu = 0.5)$p.value, 3), "exp)\n")
fwrite(met_level, file.path(add_dir, "F_metabolite_level.csv"))

cat("\nDone. Outputs in", add_dir, "\n")
