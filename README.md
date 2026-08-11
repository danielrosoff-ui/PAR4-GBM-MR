# Correlated drug-target Mendelian randomization of platelet-pathway targets on age-stratified glioblastoma risk

Reproducible code for the correlated-instrument drug-target Mendelian
randomization (MR) analysis reported in the manuscript. The pipeline tests
whether genetically proxied perturbation of five platelet-pathway drug targets
is associated with glioblastoma (GBM) risk within three age-at-onset strata.

The analysis uses published, open-source R packages and no individual-level
data. The exposure GWAS and LD reference panel are publicly available; the
age-stratified glioblastoma outcome summary statistics are not publicly
available and are held by the study team.

---

## Analysis overview

| Component | Source |
|-----------|--------|
| **Exposure** | Platelet-count GWAS (Vuckovic et al., *Cell* 2020; PMID 32888494), EBI GWAS Catalog accession **GCST90002402**; N = 408,112; European ancestry; GRCh37/hg19 |
| **Targets** | F2, F2R, F2RL3 (PAR4), TBXA2R, P2RY12 — cis windows = gene body ± 500 kb |
| **Outcome** | Age-stratified glioblastoma GWAS meta-analysis summary statistics (fixed-effect meta-analysis), three strata: `age18_53`, `age54_63`, `age64plus` |
| **LD reference** | 1000 Genomes phase 3 European panel (GRCh37) |

**Pipeline**

1. **Cis-instrument selection** (`01_extract_cis_instruments.R`)
   - Recode the exposure GWAS to TwoSampleMR format.
   - For each target gene, keep variants within ±500 kb of the gene body.
   - Threshold at **p < 5×10⁻⁶** (primary; 5×10⁻⁸ and 5×10⁻⁴ were examined as
     sensitivity analyses — change `PVAL_THRESHOLD` in `config.R`).
   - Correlated clumping at **r² = 0.2, 10,000 kb** against the 1000G EUR panel
     (retaining partially correlated cis variants for correlated MR).

2. **Outcome preparation** (`02_clean_gbm_outcomes.R`)
   - Read per-chromosome GBM meta-analysis files for each stratum.
   - Standardise variant IDs to rsIDs using a 1000G EUR reference.
   - Convert per-allele odds ratios to log-odds (beta = log(OR)) for MR.
   - Recode to TwoSampleMR outcome format; restrict to the cis instruments.

3. **Correlated MR** (`03_correlated_drug_target_mr.R`)
   - Harmonise exposure/outcome (`action = 1`) and apply Steiger filtering.
   - Compute a **local** LD correlation matrix per instrument set from the
     1000G EUR panel (avoids reliance on the remote LD server).
   - Run correlated MR (`MendelianRandomization`):
     - ≥ 3 instruments: IVW, MR-Egger, maximum-likelihood
     - 2 instruments: IVW, maximum-likelihood
     - 1 instrument (e.g. P2RY12): Wald ratio

Effect estimates are reported on the odds-ratio scale (binary GBM strata), transformed as needed.

---

## Requirements

- R (≥ 4.1)
- R packages: `TwoSampleMR`, `MendelianRandomization`, `ieugwasr`,
  `genetics.binaRies`, `data.table`, `dplyr`
- A local PLINK 1.9 binary (provided via `genetics.binaRies::get_plink_binary()`)

```r
install.packages(c("data.table", "dplyr", "MendelianRandomization"))
# remotes::install_github("MRCIEU/TwoSampleMR")
# remotes::install_github("MRCIEU/ieugwasr")
# remotes::install_github("explodecomputer/genetics.binaRies")
```

---

## Input data (not distributed here)

Place the following under `data/` (paths are configurable in `config.R`):

- `GCST90002402_platelet_count_GRCh37.tsv(.gz)` — exposure GWAS
  (download from the EBI GWAS Catalog).
- `gbm/<stratum>/meta_results_<CHR>_gbm_fixed_all.out.gz` — per-chromosome
  outcome meta-analysis files for each stratum (not publicly available). Per-allele effects are on the odds-ratio scale
  and are converted to log-odds (beta) for MR; set `OUTCOME_EFFECT` in
  `config.R` (`"or"` or `"beta"`) to match your files.
- `1000G_phase3_EUR_variants.txt(.gz)` — rsID/position reference
  (columns: chr, position [GRCh37], rsID, allele frequency).
- `1kg_v3_EUR/EUR.{bed,bim,fam}` — 1000G phase 3 EUR PLINK reference panel.

---

## Running

Edit the paths and parameters at the top of `config.R`, then:

```bash
Rscript run_all.R
```

or source the numbered scripts in order from an R session. Final results are
written to `results/correlated_mr_results.csv`; intermediate files are written
to `work/`.

---

## Files

```
config.R                          paths, parameters, target-gene coordinates
01_extract_cis_instruments.R      exposure cis-instrument selection + clumping
02_clean_gbm_outcomes.R           outcome cleaning + cis extraction
03_correlated_drug_target_mr.R    harmonisation, local LD matrix, correlated MR
run_all.R                         runs the three steps in order
```
