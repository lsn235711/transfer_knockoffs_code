suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(kableExtra))

## Directories
oak_dir <- "/oak/stanford/groups/candes"
scratch_dir <- "/scratch/groups/candes"

## Parameters
alpha <- 0.1
offset <- 0
pheno_name <- "platelet"
relatedness <- "related"
resolution_list <- 1:6

## Read the results 
all_res <- data.frame()
for(resolution in resolution_list){
  file_dir <- sprintf("%s/transfer_gwas/gwas/replicate/%s_res%d_fdr%.1f_offset%d_%s.txt", 
                      oak_dir, pheno_name, resolution, alpha, offset, relatedness)
  res <- read_delim(file_dir, delim = " ", col_types = cols()) 
  all_res <- rbind(all_res, res)
}

method_list <- c("Vanilla Knockoffs", "Invariant Knockoffs", "Adaptive Knockoffs", "Weighted Lasso")
new_method_list <- c("Vanilla Knockoffs", "Knockoffs-LRO", "Adaptive Knockoffs", "Weighted Lasso")
resolution_list <- 1:6
new_resolution_list <- c(3, 20, 41, 81, 208, 425)



pop_list <- c("whitenonbritish", "britishindia", "black", "asian")
new_pop_list <- c("European", "Indian", "African", "Asian")

all_res <- all_res %>% 
  filter(method != "Adaptive Weighted Knockoffs") %>%
  mutate(method = factor(method, method_list, new_method_list)) %>%
  mutate(pop = factor(pop, pop_list, new_pop_list)) %>%
  mutate(res = factor(res, resolution_list, new_resolution_list)) %>%
  select(-pheno) %>%
  spread(method, rej) %>%
  rename(Population = pop, Resolution = res)
all_res
#    arrange(res) %>%
#    mutate()

all_res %>% kbl(booktabs = TRUE, "latex") %>%
  collapse_rows(columns = 1:2)
    

