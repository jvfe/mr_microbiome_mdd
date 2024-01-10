library(TwoSampleMR)
library(stringr)
library(dplyr)
library(vroom)
library(tidyr)

read_depression <- function(path) {
  exp <- read_exposure_data(
    path,
    sep = "\t",
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    pval_col = "P",
    other_allele_col = "A1",
    effect_allele_col = "A2",
    samplesize_col = "N"
  )

  exposure_filtered <- exp[which(exp$pval.exposure < 0.0000001), ]

  clumped <- clump_data(exposure_filtered, clump_r2 = 0.3)

  return(clumped)
}

read_bacteria <- function(path) {
  read_outcome_data(
    path,
    sep = "\t",
    snp_col = "rsID",
    beta_col = "beta",
    se_col = "SE",
    pval_col = "P.weightedSumZ",
    other_allele_col = "ref.allele",
    effect_allele_col = "eff.allele",
    samplesize_col = "N",
    phenotype_col = "bac"
  )
}

run_MR <- function(exposure_data, outcome_data) {

  message(paste("Running", exposure_data, "against", outcome_data))

  mdd <- read_depression(exposure_data)

  bacteria <- read_bacteria(outcome_data)

  dat <-
    harmonise_data(exposure_dat = mdd, outcome_dat = bacteria)

  df <- add_rsq(dat)

  outcome <- bacteria$outcome[1]

  exposure_name <- paste(
    str_extract(exposure_data, "\\d{4}"),
    str_extract(exposure_data, "(\\w*male|both_sexes)"),
    "MDD",
    sep = "_"
  )

  result_name <- paste0(exposure_name, "-", gsub("\\.", "_", outcome))

  vroom_write(df, file = paste0("results/harmonised/", result_name, ".txt.gz"))

  mr_report(
    dat,
    output_path = paste0("results/reports/", result_name),
    output_type = "md"
  )
}

depression_exps <- fs::dir_ls("results/clean_sumstats/")

bacs <-
  fs::dir_ls("data/mibiogen/", glob = "*txt.gz", recurse = TRUE)

combs <- tidyr::expand_grid(depression_exps, bacs)
colnames(combs) <- c("exposure_data", "outcome_data")

purrr::pmap(combs, run_MR)

# run_MR(
#   exposure_data = "results/clean_sumstats/2090.gwas.imputed_v3.both_sexes.tsv.gz",
#   outcome_data = "data/mibiogen/order/order.Actinomycetales.id.420.summary.txt.gz",
#   exposure_name = "2090_both_MDD"
# )
