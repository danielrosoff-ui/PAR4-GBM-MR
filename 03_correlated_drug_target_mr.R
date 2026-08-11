# =============================================================================
# 03_correlated_drug_target_mr.R
# Correlated-instrument drug-target Mendelian randomization (MR) of five
# platelet-pathway targets on age-stratified glioblastoma risk.
#
# For each gene x outcome-stratum pair:
#   1. Harmonise exposure and outcome (TwoSampleMR, action = 1) and apply
#      Steiger filtering.
#   2. Compute a local LD correlation matrix for the instrument set from the
#      1000G EUR reference panel (variant server-independent).
#   3. Run correlated MR (MendelianRandomization package):
#        >= 3 SNPs : IVW, Egger, maximum-likelihood
#        == 2 SNPs : IVW, maximum-likelihood
#        == 1 SNP  : Wald ratio (single-instrument, uncorrelated)
#
# Output: results/correlated_mr_results.csv  (tidy, one row per method/estimate)
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(TwoSampleMR)
  library(MendelianRandomization)
  library(ieugwasr)
  library(genetics.binaRies)
})
source("config.R")

plink_bin <- genetics.binaRies::get_plink_binary()

# ---- Load harmonised inputs -------------------------------------------------
exposure <- as.data.frame(fread(file.path(WORK_DIR, "cis_instruments.csv")))
outcome  <- as.data.frame(fread(file.path(WORK_DIR, "gbm_outcomes_cis.csv")))

# Harmonise all gene x outcome data at once, then Steiger-filter.
harmonised <- harmonise_data(exposure, outcome, action = 1)
harmonised <- steiger_filtering(harmonised)
fwrite(harmonised, file.path(WORK_DIR, "harmonised_data.csv"))

# -----------------------------------------------------------------------------
# Local LD correlation matrix + allele alignment.
# ieugwasr::ld_matrix returns SNPs in reference-panel order and codes the
# correlation to its own reference alleles, so we (a) subset to SNPs present in
# the panel, (b) recover the panel's coding alleles, and (c) flip exposure and
# outcome betas so both are expressed relative to the panel's A1. The LD matrix
# and the beta vectors are finally reordered to a single common SNP order.
# -----------------------------------------------------------------------------
build_ld_and_align <- function(d) {
  snps <- unique(d$SNP)

  ld   <- ieugwasr::ld_matrix(variants = snps, with_alleles = FALSE,
                              plink_bin = plink_bin, bfile = LD_REFERENCE)
  ld_a <- ieugwasr::ld_matrix(variants = snps, with_alleles = TRUE,
                              plink_bin = plink_bin, bfile = LD_REFERENCE)

  panel_snps <- colnames(ld)                 # SNPs actually in the reference
  d <- d[d$SNP %in% panel_snps, , drop = FALSE]

  # coding alleles used by the LD matrix, e.g. "rs123_A_G" -> SNP, A1, A2
  tags <- colnames(ld_a)
  allele_map <- data.frame(
    SNP = sub("_.*", "", tags),
    A1  = sub(".*_(.)_.*", "\\1", tags),
    A2  = sub(".*_._(.)$", "\\1", tags),
    stringsAsFactors = FALSE
  )
  m <- merge(d, allele_map, by = "SNP")

  # Align both betas to the LD-panel A1 allele.
  flip <- m$effect_allele.exposure != m$A1
  m$beta.exposure_aligned <- ifelse(flip, -m$beta.exposure, m$beta.exposure)
  m$beta.outcome_aligned  <- ifelse(flip, -m$beta.outcome,  m$beta.outcome)

  # Put the LD matrix and the aligned betas in one identical SNP order.
  ord <- m$SNP
  list(m = m, corr = ld[ord, ord, drop = FALSE])
}

# tidy one-row-per-estimate extractor for MendelianRandomization objects -------
tidy_ivw <- function(o, gene, outcome, method) data.frame(
  gene = gene, outcome = outcome, method = method, nsnp = o@SNPs,
  estimate = o@Estimate, se = o@StdError,
  ci_lower = o@CILower, ci_upper = o@CIUpper, pval = o@Pvalue,
  egger_intercept = NA_real_, egger_intercept_pval = NA_real_,
  stringsAsFactors = FALSE)

tidy_egger <- function(o, gene, outcome) data.frame(
  gene = gene, outcome = outcome, method = "MR-Egger (correlated)", nsnp = o@SNPs,
  estimate = o@Estimate, se = o@StdError.Est,
  ci_lower = o@CILower.Est, ci_upper = o@CIUpper.Est, pval = o@Pvalue.Est,
  egger_intercept = o@Intercept, egger_intercept_pval = o@Pvalue.Int,
  stringsAsFactors = FALSE)

# ---- Correlated MR for one gene x outcome -----------------------------------
run_correlated_mr <- function(d, gene, outcome) {
  n <- nrow(d)
  if (n == 0) return(NULL)

  # Single instrument (e.g. P2RY12): standard uncorrelated Wald ratio.
  if (n == 1) {
    d$mr_keep <- TRUE
    w <- TwoSampleMR::mr(d, method_list = "mr_wald_ratio")
    return(data.frame(
      gene = gene, outcome = outcome, method = "Wald ratio", nsnp = 1,
      estimate = w$b, se = w$se, ci_lower = w$b - 1.96 * w$se,
      ci_upper = w$b + 1.96 * w$se, pval = w$pval,
      egger_intercept = NA_real_, egger_intercept_pval = NA_real_,
      stringsAsFactors = FALSE))
  }

  a  <- build_ld_and_align(d)
  bx <- a$m$beta.exposure_aligned; bxse <- a$m$se.exposure
  by <- a$m$beta.outcome_aligned;  byse <- a$m$se.outcome
  input <- MendelianRandomization::mr_input(bx, bxse, by, byse, corr = a$corr)

  ivw <- MendelianRandomization::mr_ivw(input)
  mxl <- MendelianRandomization::mr_maxlik(input)
  res <- rbind(
    tidy_ivw(ivw, gene, outcome, "IVW (correlated)"),
    tidy_ivw(mxl, gene, outcome, "Maximum likelihood (correlated)")
  )
  # Egger requires >= 3 instruments.
  if (nrow(a$m) >= 3) {
    egg <- MendelianRandomization::mr_egger(input)
    res <- rbind(res, tidy_egger(egg, gene, outcome))
  }
  res
}

# ---- Loop over gene x outcome -----------------------------------------------
results <- list()
for (oc in unique(harmonised$id.outcome)) {
  for (g in unique(harmonised$id.exposure)) {
    d <- harmonised[harmonised$id.exposure == g &
                    harmonised$id.outcome  == oc, , drop = FALSE]
    r <- tryCatch(run_correlated_mr(d, g, oc),
                  error = function(e) { message(sprintf(
                    "  error %s x %s: %s", g, oc, conditionMessage(e))); NULL })
    if (!is.null(r)) results[[paste(g, oc, sep = "_")]] <- r
  }
}

results_tbl <- bind_rows(results)
# Report odds ratios (GBM outcomes are binary case/control strata).
results_tbl <- results_tbl %>%
  mutate(OR = exp(estimate), OR_lower = exp(ci_lower), OR_upper = exp(ci_upper))

out_file <- file.path(RESULTS_DIR, "correlated_mr_results.csv")
fwrite(results_tbl, out_file)
message("Wrote -> ", out_file)
print(results_tbl[, c("gene", "outcome", "method", "nsnp",
                      "OR", "OR_lower", "OR_upper", "pval")])
