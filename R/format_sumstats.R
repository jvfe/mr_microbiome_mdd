BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh37")
BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")
install.packages("vroom")

library(MungeSumstats)
library(vroom)
library(dplyr)

phenotype <- vroom("data/head_female.tsv")

phenotype_rsids <- vroom("data/variants.tsv.bgz") %>%
  filter(variant %in% phenotype$variant) %>%
  right_join(phenotype, by = "variant")

gc()

reformatted <-
  MungeSumstats::format_sumstats(path=phenotype_rsids,
                                 ref_genome="GRCh37",
                                 nThread = 8,
                                 save_path = "results/2090_clean_female.tsv.gz")

