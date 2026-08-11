# =============================================================================
# 02_clean_gbm_outcomes.R
# Clean the age-stratified glioblastoma (GBM) GWAS meta-analysis summary
# statistics and recode to TwoSampleMR outcome format, restricted to the cis
# instrument variants produced by 01_extract_cis_instruments.R.
#
# For each age stratum (age18_53, age54_63, age64plus):
#   1. Read the per-chromosome meta-analysis files.
#   2. Standardise variant IDs to rsIDs using a 1000G EUR reference
#      (rows already carrying an rsID are matched by rsID; rows with a
#      chr:pos-style ID are matched by chromosome + position).
#   3. Recode to TwoSampleMR outcome format and keep only cis instrument SNPs.
#
# Output: work/gbm_outcomes_cis.csv  (all strata stacked; id.outcome = stratum).
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})
source("config.R")

# ---- cis instruments from step 01 -------------------------------------------
cis <- fread(file.path(WORK_DIR, "cis_instruments.csv"))
cis_snps <- unique(cis$SNP)
message(length(cis_snps), " cis instrument SNPs to extract from outcomes.")

# ---- rsID / position reference ----------------------------------------------
# Columns: V1 = chr, V2 = pos (GRCh37), V3 = rsID, V6 = allele frequency
ref <- fread(RSID_REFERENCE)
ref[, V1 := as.character(V1)]

O <- OUTCOME_COLS

clean_one_chromosome <- function(file_path) {
  df <- as.data.frame(fread(file_path))
  if (!nrow(df)) return(NULL)

  # Split rows that already have an rsID from those that do not.
  has_rs <- grepl("^rs", df[[O$snp]])
  with_rs    <- df[has_rs, , drop = FALSE]
  without_rs <- df[!has_rs, , drop = FALSE]

  # (a) rows with rsID: attach reference chr/AF by rsID
  if (nrow(with_rs)) {
    with_rs <- merge(with_rs, ref, by.x = O$snp, by.y = "V3")
    with_rs$SNP <- with_rs[[O$snp]]
    with_rs$chr <- as.character(with_rs$V1)
  }
  # (b) rows without rsID: parse chr from the ID and match by chr + position
  if (nrow(without_rs)) {
    without_rs$chr <- sapply(strsplit(without_rs[[O$snp]], ":"), `[`, 1)
    without_rs <- merge(without_rs, ref,
                        by.x = c("chr", O$pos), by.y = c("V1", "V2"))
    without_rs$SNP <- without_rs$V3
  }

  temp <- bind_rows(with_rs, without_rs)
  if (!nrow(temp)) return(NULL)

  # Outcome effect on the log-odds (beta) scale required by MR. GBM strata are
  # binary; if the meta-analysis reports odds ratios, convert beta = log(OR).
  # (SE is assumed to be on the log-odds scale, as produced by logistic
  # meta-analysis; adjust here if your file stores it differently.)
  beta_out <- if (identical(OUTCOME_EFFECT, "or")) log(temp[[O$or]]) else temp[[O$beta]]

  # Recode to TwoSampleMR outcome format.
  data.frame(
    SNP                   = temp$SNP,
    chr.outcome           = as.character(temp$chr),
    pos.outcome           = temp[[O$pos]],
    effect_allele.outcome = temp[[O$effect]],
    other_allele.outcome  = temp[[O$other]],
    beta.outcome          = beta_out,
    se.outcome            = temp[[O$se]],
    pval.outcome          = temp[[O$pval]],
    eaf.outcome           = temp$V6,
    stringsAsFactors      = FALSE
  )
}

clean_stratum <- function(stratum, dir_path) {
  files <- list.files(dir_path,
                      pattern = "meta_results_\\d+_gbm_fixed_all\\.out\\.gz$",
                      full.names = TRUE)
  if (!length(files)) {
    message("  no per-chromosome files found in ", dir_path); return(NULL)
  }
  cleaned <- lapply(files, clean_one_chromosome) %>% bind_rows()
  cleaned <- cleaned[cleaned$SNP %in% cis_snps, , drop = FALSE]
  if (nrow(cleaned)) {
    cleaned$outcome    <- stratum
    cleaned$id.outcome <- stratum
  }
  message(sprintf("  %-9s -> %d cis SNPs extracted", stratum, nrow(cleaned)))
  cleaned
}

message("Cleaning GBM outcome strata...")
outcomes <- Map(clean_stratum, names(OUTCOME_STRATA), OUTCOME_STRATA) %>%
  bind_rows()

out_file <- file.path(WORK_DIR, "gbm_outcomes_cis.csv")
fwrite(outcomes, out_file)
message("Wrote -> ", out_file)
