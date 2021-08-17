## Arguments
args <- commandArgs(trailingOnly=TRUE)
resolution <- as.integer(args[1])
phenotype <- args[2]

## Load packages
suppressMessages(library(tidyverse))
suppressMessages(library(kableExtra))
suppressMessages(library(SNPknock))
suppressMessages(library(adaptiveKnockoff))
suppressMessages(library(randomForest))
suppressMessages(library(gam))

## Directories
oak_dir <- "/oak/stanford/groups/candes"
scratch_dir <- "/scratch/groups/candes"

## Load partition information
read_partitions <- function(chr, resolution){
  in.file <- sprintf("/oak/stanford/groups/candes/popstruct/data/partitions/hap_chr%d.txt", chr)
  Partitions <- read_delim(in.file, delim=" ", col_types = cols())
  res.name <- sprintf("res_%d", resolution)
  groups <- Partitions %>% mutate(CHR = chr) %>%
  rename(Group = res.name) %>% 
  select(CHR, SNP, Group)
  return(groups)
}

## Function to generate results
filter <- function(pheno_name, population, resolution, relatedness, prop,
                   alpha = 0.1, offset = 1, diagnostic = FALSE){
  ## Initializa the output
  output <- data.frame()

  ## Load W
  stat_file <- sprintf("%s/popstruct/analysis/stats/lasso_%s_%s_%s_res%d.txt", oak_dir, pheno_name, population, relatedness, resolution) 
  stat <- read_delim(stat_file, delim = " ", col_types = cols())

  ## Load the side information
  u_file <- sprintf("%s/popstruct/analysis/stats/lasso_%s_british_%s_res%d.txt", oak_dir, pheno_name, relatedness, resolution)
  U <- read_delim(u_file, delim = " ", col_types = cols()) %>%
    rename(U = W) %>%
    mutate(U = abs(U))
  

  ## Match W with the side information
  stat <- stat %>% left_join(U, by = c("CHR", "Group"))
  stat$U[is.na(stat$U)] <- 0
  
  
  ############################
  ## Run adaptive knockoffs ##
  ############################
  res <- filter_gam(stat$W, stat$U, alpha = alpha, offset = offset, mute = TRUE, df = 2,
                                      reveal_prop = prop)
  output <- rbind(output, data.frame(rej = length(res$rejs[[1]]), method = "Adaptive Knockoffs"))

  ########################
  ## Run mean knockoffs ##
  ########################
  stat <- stat %>% mutate(W1 = sign(W) * (0.9 * abs(W) + 0.1 * abs(U)))
  tau <- knockoff::knockoff.threshold(stat$W1, fdr = alpha, offset = offset)
  res <- which(stat$W1 >= tau)
  output <- rbind(output, data.frame(rej = length(res), method = "Invariant Knockoffs"))

  ###########################
  ## Run vanilla knockoffs ##
  ###########################
  tau <- knockoff::knockoff.threshold(stat$W, fdr = alpha, offset = offset)
  res <- which(stat$W >= tau)
  output <- rbind(output, data.frame(rej = length(res), method = "Vanilla Knockoffs"))

  ##################################################
  ## Run knockoffs with weighted lasso statistics ##
  ##################################################
  new_population <- factor(population, population_list, new_population_list) %>% as.character()
  wstat_file <- sprintf("%s/transfer_gwas/gwas/stats/lasso_transfer_res%d_%s_%s.txt", 
                        oak_dir, resolution, new_population, pheno_name) 
  wstat <- read_delim(wstat_file, delim = " ", col_types = cols())
  if(resolution !=0){
    Groups <- lapply(1:22, read_partitions, resolution = resolution)
    Groups <- do.call(rbind, Groups)
  }

  if(resolution != 0){
    wstat <- wstat %>% 
      mutate(Knockoff = endsWith(SNP, ".k"), SNP = str_replace(SNP, ".k", "")) %>%
      left_join(Groups, by=c("SNP", "CHR"))  %>%
      group_by(CHR, Group, Knockoff) %>%
      summarise(SNP.lead=SNP[which.max(abs(Z))], Z=sum(abs(Z)), Size = n()) %>%
      group_by(CHR, Group) %>%
      summarise(SNP.lead=SNP.lead[1], W=sum(Z[!Knockoff])-sum(Z[Knockoff]), Size = sum(Size)) %>%
      ungroup() %>%
      arrange(desc(abs(W)))
  }else{
    wstat <- wstat %>% 
      mutate(Knockoff = endsWith(SNP, ".k"), SNP = str_replace(SNP, ".k", "")) %>%
      group_by(CHR, SNP) %>%
      summarise(W=sum(Z[!Knockoff])-sum(Z[Knockoff])) %>%
      ungroup() %>%
      arrange(desc(abs(W)))
  }

  tau <- knockoff::knockoff.threshold(wstat$W, fdr = alpha, offset = offset)
  res <- which(wstat$W >= tau)
  output <- rbind(output, data.frame(rej = length(res), method = "Weighted Lasso"))


  ###################
  ## Collet output ##
  ###################

  output$pheno <- pheno_name
  output$pop <- population
  output$res <- resolution 
  
  return(output)
}

## parameters
population_list <- c("whitenonbritish", "britishindia", "asian", "black")
new_population_list <- c("European", "Indian", "Asian", "African")
alpha <- 0.1
offset <- 0
prop <- 0.5
relatedness <- "related"

all_res <- data.frame()
for (population in population_list){
    cat(sprintf("Computing results for %s in the %s population...", phenotype, population))
    res <- lapply(resolution:resolution, filter, pheno_name = phenotype, population = population, 
                  relatedness = relatedness, prop = prop, alpha = alpha, offset = offset)
    res <- do.call(rbind, res)
    all_res <- rbind(all_res, res)
    cat("done.\n")
}
#all_res
out_file <- sprintf("%s/transfer_gwas/gwas/replicate/%s_res%d_fdr%s_offset%d_%s.txt", 
                    oak_dir, phenotype, resolution, alpha, offset, relatedness)
write_delim(all_res, out_file, delim = " ")

#all_res %>% kbl(booktabs = TRUE, "latex")





