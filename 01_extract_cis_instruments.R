# =============================================================================
# 01_extract_cis_instruments.R
# Build cis genetic instruments for five platelet-pathway drug targets from the
# platelet-count GWAS (exposure).
#
# Steps:
#   1. Read the platelet-count GWAS and recode to TwoSampleMR exposure format.
#   2. For each target gene, subset variants within +/- CIS_WINDOW_BP of the
#      gene body.
#   3. Threshold at PVAL_THRESHOLD.
#   4. Correlated clumping (r2 = CLUMP_R2, kb = CLUMP_KB) against the 1000G EUR
#      reference panel.
#
# Output: work/cis_instruments.csv  (one row per retained instrument, tagged by
#         gene in id.exposure).
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(TwoSampleMR)
  library(ieugwasr)
  library(genetics.binaRies)
})
source("config.R")

plink_bin <- genetics.binaRies::get_plink_binary()

# ---- 1. Read exposure GWAS, recode to TwoSampleMR exposure format ------------
message("Reading exposure GWAS: ", EXPOSURE_GWAS)
raw <- fread(EXPOSURE_GWAS)
C   <- EXPOSURE_COLS
exposure <- data.frame(
  SNP                    = raw[[C$snp]],
  chr.exposure           = as.integer(raw[[C$chr]]),
  pos.exposure           = as.integer(raw[[C$pos]]),
  effect_allele.exposure = raw[[C$effect]],
  other_allele.exposure  = raw[[C$other]],
  beta.exposure          = raw[[C$beta]],
  se.exposure            = raw[[C$se]],
  pval.exposure          = raw[[C$pval]],
  eaf.exposure           = raw[[C$eaf]],
  samplesize.exposure    = EXPOSURE_N,
  stringsAsFactors       = FALSE
)
exposure <- exposure[!is.na(exposure$SNP) & exposure$SNP != "", ]
message("  ", nrow(exposure), " variants read.")

# ---- 2. cis window per gene -------------------------------------------------
extract_cis <- function(g) {
  lo <- g$start - CIS_WINDOW_BP
  hi <- g$end   + CIS_WINDOW_BP
  cis <- subset(exposure,
                chr.exposure == g$chr &
                pos.exposure >= lo &
                pos.exposure <= hi)
  if (nrow(cis)) {
    cis$id.exposure <- g$gene
    cis$exposure    <- g$label
  }
  message(sprintf("  %-7s chr%s:%d-%d  ->  %d cis variants",
                  g$gene, g$chr, lo, hi, nrow(cis)))
  cis
}
cis_all <- do.call(rbind, lapply(seq_len(nrow(TARGET_GENES)),
                                 function(i) extract_cis(TARGET_GENES[i, ])))

# ---- 3. p-value threshold ---------------------------------------------------
cis_sig <- subset(cis_all, pval.exposure < PVAL_THRESHOLD)
message(sprintf("cis variants at p < %.0e: %d", PVAL_THRESHOLD, nrow(cis_sig)))

# ---- 4. correlated clumping (per gene) --------------------------------------
# Uses ieugwasr::ld_clump with a local PLINK binary and 1000G EUR reference.
# clump_p = 1 keeps all variants below PVAL_THRESHOLD as clumping candidates;
# the correlated r2 = CLUMP_R2 retains a set of partially correlated cis
# variants suitable for correlated MR (as opposed to the usual r2 = 0.001).
clump_gene <- function(df) {
  cand <- df %>%
    distinct(SNP, .keep_all = TRUE) %>%
    mutate(rsid = SNP, pval = pval.exposure, id = id.exposure)
  keep <- tryCatch(
    ieugwasr::ld_clump(dat = cand,
                       clump_kb = CLUMP_KB, clump_r2 = CLUMP_R2, clump_p = 1,
                       plink_bin = plink_bin, bfile = LD_REFERENCE),
    error = function(e) { message("  clump error (", df$id.exposure[1], "): ",
                                  conditionMessage(e)); cand[0, ] })
  df[df$SNP %in% keep$SNP, ]
}
cis_clumped <- cis_sig %>%
  group_split(id.exposure) %>%
  lapply(clump_gene) %>%
  bind_rows()

message("Retained instruments after clumping:")
print(as.data.frame(table(cis_clumped$id.exposure)))

out_file <- file.path(WORK_DIR, "cis_instruments.csv")
fwrite(cis_clumped, out_file)
message("Wrote -> ", out_file)
