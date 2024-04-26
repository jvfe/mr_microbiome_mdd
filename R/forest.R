library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)

make_forest <- function(data) {

  ggplot(data = data,
         aes(
           y = p1,
           x = estimate,
           xmin = estimate - stderr,
           xmax = estimate + stderr
         )) +
    geom_point(aes(size = psize)) +
    scale_size_continuous(range = c(0, 2)) +
    geom_errorbarh(height = 0.5) +
    labs(x = 'Effect Size', y = NULL) +
    geom_vline(
      xintercept = 0,
      color = 'black',
      linetype = 'dashed',
      alpha = .5
    ) +
    guides(size = FALSE) +
    theme_classic()

}


full <- vroom::vroom("results/full.tsv") |>
  filter(nsnp > 0)

col_for_plot <- c("p1", "p2", "estimate", "stderr", "p", "NSNPs")

filtered <- full |>
  setNames(col_for_plot) |>
  mutate(Q=p.adjust(p,method="BH")) |>
  filter(Q<0.05) |>
  mutate(
    upper = estimate + stderr,
    lower = estimate - stderr
  )

ci.lb <-  filtered$lower
ci.ub <-  filtered$upper

#taken from forest.default source code
level = 95
alpha <- ifelse(level > 1, (100 - level) / 100, 1 - level)
vi <- ((ci.ub - ci.lb) / (2 * qnorm(alpha / 2, lower.tail = FALSE))) ^
  2
wi <- 1 / sqrt(vi)
psize <- wi / sum(wi, na.rm = TRUE)
psize <- (psize - min(psize, na.rm = TRUE)) / (max(psize,
                                                   na.rm = TRUE) - min(psize, na.rm = TRUE))
filtered$psize <- (psize * 1) + 0.5

get_specific_taxon <- function(table, taxon) {
  if (str_detect(taxon, "s__\\w+")) {
    return(taxon)
  } else {
    maxed <- table |>
      filter(str_detect(p1, taxon)) |>
      filter(str_length(p1) == max(str_length(p1))) |>
      pull(p1)

    return(maxed)
  }
}

pathways <- filtered |>
  filter(str_starts(p1, "GWAS_"))

taxa <- filtered |>
  filter(!(p1 %in% pathways$p1)) %>%
  mutate(specific = purrr::map_chr(unique(p1), get_specific_taxon, table = .)) |>
  filter(p1 == specific)

make_forest(pathways |> mutate(p1 = str_remove_all(p1, "GWAS_")))

ggsave(
  "results/pathways_plot.png",
  units = "in",
  width = 9,
  dpi = 300
)

make_forest(taxa |> mutate(p1 = str_split_i(p1, "_\\w__", i = -1)))

ggsave(
  "results/taxa_plot.png",
  units = "in",
  dpi = 300
)




