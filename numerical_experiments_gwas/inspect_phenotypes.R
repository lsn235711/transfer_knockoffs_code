#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

pop.list <- c("African", "Asian", "European", "Indian")
n.causal <- 100
spec.list <- c(0, 0.5, 1)

suppressMessages(library(tidyverse))

oak <- "/oak/stanford/groups/candes/matteo/transfer_knockoffs"
scratch <- "/scratch/groups/candes/matteo/transfer_knockoffs"

in.dir <- "/oak/stanford/groups/candes/matteo/transfer_knockoffs/data/phenotypes"

Phenotypes <- do.call("rbind", lapply(spec.list, function(spec) {
  do.call("rbind", lapply(pop.list, function(pop) {
    in.file <- sprintf("%s/phenotypes_n%s_s%s_%s.txt", in.dir, n.causal, round(100*spec), pop)
    read_delim(in.file, delim=" ", col_types=cols()) %>% mutate(Population = pop, Specificity = spec)
  }))
}))


Phenotypes %>%
  ggplot(aes(x=Population, y=Y1, color=Population)) +
  geom_boxplot() +
  facet_grid(Specificity ~ .) +
  theme_bw()
