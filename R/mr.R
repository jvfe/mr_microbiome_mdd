library(TwoSampleMR)
library(dplyr)
library(vroom)
library(tidyr)

actinobacteria <- vroom("data/phylum.Actinobacteria.id.400.summary.txt.gz")


mdd_female <- vroom("data/2090.gwas.imputed_v3.female.tsv")


mdd_female <- mdd_female %>%
  separate(variant, into = c("chr", "pos", "other_allele", "effect_allele"), remove = F)

# get rsIDs
# calculate effect allele frequency




mdd_female <- read_exposure_data(
  "data/2090.gwas.imputed_v3.female.tsv",
  sep = "\t",
  snp_col = "",
  beta_col = "beta",
  se_col = "se",
  pval_col = "pval",
  effect_allele_col = "",
  other_allele_col = "",

  samplesize_col = "n_complete_samples",
)