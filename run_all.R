# =============================================================================
# run_all.R
# Run the full correlated drug-target MR pipeline end to end.
# Edit paths/parameters in config.R first, then:  Rscript run_all.R
# =============================================================================
source("01_extract_cis_instruments.R")
source("02_clean_gbm_outcomes.R")
source("03_correlated_drug_target_mr.R")
message("\nPipeline complete. See results/correlated_mr_results.csv")
