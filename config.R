# =============================================================================
# config.R
# Shared configuration for the correlated drug-target Mendelian randomization
# (MR) analysis of platelet-pathway drug targets on age-stratified
# glioblastoma (GBM) risk.
#
# Edit the paths in this file to point at your local copies of the input data
# and the 1000 Genomes LD reference panel, then run the numbered scripts in
# order (or source run_all.R).
# =============================================================================

# ---- Directory layout -------------------------------------------------------
# All inputs/outputs are addressed relative to PROJECT_DIR so the pipeline is
# portable. Set this to wherever you unpack the release.
PROJECT_DIR <- normalizePath(".", mustWork = FALSE)

DATA_DIR    <- file.path(PROJECT_DIR, "data")        # input GWAS summary stats
WORK_DIR    <- file.path(PROJECT_DIR, "work")        # intermediate files
RESULTS_DIR <- file.path(PROJECT_DIR, "results")     # final MR result tables
for (d in c(DATA_DIR, WORK_DIR, RESULTS_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ---- Input data -------------------------------------------------------------
# 1) EXPOSURE: platelet-count GWAS (Vuckovic et al. 2020, PMID 32888494),
#    EBI GWAS Catalog accession GCST90002402, N = 408,112, European ancestry,
#    genome build GRCh37/hg19. Download the harmonised summary statistics from
#    https://www.ebi.ac.uk/gwas/ and place the (optionally gzipped) TSV here.
EXPOSURE_GWAS <- file.path(DATA_DIR, "GCST90002402_platelet_count_GRCh37.tsv.gz")

# Column names in the exposure file. Defaults match the EBI GWAS Catalog
# harmonised format; change the right-hand side if your download differs.
EXPOSURE_COLS <- list(
  snp         = "variant_id",
  chr         = "chromosome",
  pos         = "base_pair_location",
  effect      = "effect_allele",
  other       = "other_allele",
  beta        = "beta",
  se          = "standard_error",
  pval        = "p_value",
  eaf         = "effect_allele_frequency"
)
EXPOSURE_N <- 408112

# 2) OUTCOME: age-stratified glioblastoma GWAS meta-analysis summary
#    statistics (study-internal; fixed-effect meta-analysis across contributing
#    cohorts), one set of per-chromosome files per age-at-onset stratum.
#    Provide, for each stratum, a directory of per-chromosome gzipped files
#    named "meta_results_<CHR>_gbm_fixed_all.out.gz".
OUTCOME_STRATA <- list(
  age18_53  = file.path(DATA_DIR, "gbm", "age18_53"),
  age54_63  = file.path(DATA_DIR, "gbm", "age54_63"),
  age64plus = file.path(DATA_DIR, "gbm", "age64plus")
)
# Column names in the raw per-chromosome outcome files.
# The outcome is a binary (case/control) GBM stratum. If the meta-analysis
# reports the per-allele effect as an ODDS RATIO, set OUTCOME_EFFECT = "or"
# and give the OR column name in `or`; the cleaning step then converts to the
# log-odds (beta) scale required by MR (beta = log(OR)). If the effect is
# already on the log-odds/beta scale, set OUTCOME_EFFECT = "beta".
OUTCOME_EFFECT <- "or"       # "or"  -> beta = log(OR);  "beta" -> use as-is
OUTCOME_COLS <- list(
  snp    = "rsid",
  pos    = "pos",
  effect = "allele_A",     # coded/effect allele
  other  = "allele_B",
  beta   = "beta",         # used when OUTCOME_EFFECT == "beta"
  or     = "OR",           # used when OUTCOME_EFFECT == "or"
  se     = "se",           # standard error on the log-odds (beta) scale
  pval   = "P_value",
  eaf    = "coded_af"       # coded-allele frequency
)

# rsID / position reference used to standardise outcome variant IDs. A plain
# table of 1000 Genomes phase 3 EUR variants with columns:
#   V1 = chromosome, V2 = position (GRCh37), V3 = rsID, V6 = allele frequency
RSID_REFERENCE <- file.path(DATA_DIR, "1000G_phase3_EUR_variants.txt.gz")

# 3) LD REFERENCE: 1000 Genomes phase 3 European PLINK bfile (GRCh37), given as
#    a path prefix (without .bed/.bim/.fam). Used for correlated clumping and
#    for the local LD correlation matrix.
LD_REFERENCE <- file.path(DATA_DIR, "1kg_v3_EUR", "EUR")

# ---- Analysis parameters ----------------------------------------------------
CIS_WINDOW_BP  <- 500000   # +/- window around each gene body (base pairs)
PVAL_THRESHOLD <- 5e-06    # cis-variant significance threshold (primary)
CLUMP_R2       <- 0.2      # correlated clumping r^2 (NOT the usual 0.001)
CLUMP_KB       <- 10000    # clumping window (kb)

# ---- Target genes -----------------------------------------------------------
# Five platelet-pathway drug targets. Coordinates are gene body start/end on
# GRCh37/hg19 (Ensembl). PAR4 is the protein; its gene is F2RL3.
TARGET_GENES <- data.frame(
  gene  = c("F2",       "F2R",      "F2RL3",    "TBXA2R",   "P2RY12"),
  label = c("F2",       "F2R",      "PAR4",     "TBXA2R",   "P2RY12"),
  chr   = c(11,          5,          19,         19,          3),
  start = c(46740730,    76011868,   16999671,   3594504,     151055168),
  end   = c(46761056,    76031606,   17003417,   3606838,     151102600),
  stringsAsFactors = FALSE
)

message("config.R loaded: ", nrow(TARGET_GENES), " target genes, ",
        length(OUTCOME_STRATA), " outcome strata.")
