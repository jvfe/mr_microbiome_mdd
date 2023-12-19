options(timeout = 10000)
BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh37")
BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")
install.packages("vroom")

library(MungeSumstats)
library(vroom)
library(dplyr)
library(fs)

read_ukbb_sumstat <-
  function(filepath,
           refgenome = "GRCh37",
           threads = 8,
           ...) {
    phenotype_sumstat <- vroom(filepath)

    # Get filepath, result path will be the same but .gz
    filename <- path_ext_set(path_file(filepath), ext = "gz")

    # Join in with rsIDs and etc
    phenotype_rsids <- vroom("data/variants.tsv.bgz") %>%
      filter(variant %in% phenotype_sumstat$variant) %>%
      right_join(phenotype_sumstat, by = "variant")

    gc()

    reformatted <-
      MungeSumstats::format_sumstats(
        path = phenotype_rsids,
        ref_genome = refgenome,
        nThread = threads,
        save_path = paste0("results/clean_sumstats/", filename),
        ...
      )

    return(reformatted)
  }

sumstats_list <- dir_ls("data/neal_depression/howard_phenotypes/", glob = "*tsv.bgz")

purrr::walk(sumstats_list, read_ukbb_sumstat, threads = 16)
