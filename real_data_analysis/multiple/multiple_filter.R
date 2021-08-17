#!/usr/bin/env Rscript

## Arguments
args <- commandArgs(trailingOnly=TRUE)
resolution <- as.integer(args[1])

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
source("filter_gam.R")

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
  if(pheno_name %in% c("bmi", "diabetes") | population == "whitenonbritish"){
    stat_file <- sprintf("%s/popstruct/analysis/stats/lasso_%s_%s_%s_res%d.txt", oak_dir, pheno_name, population, relatedness, resolution) 
    stat <- read_delim(stat_file, delim = " ", col_types = cols())
  }else{
    new_population <- factor(population, population_list, new_population_list) %>% as.character()
    stat_file <- sprintf("%s/transfer_gwas/gwas/stats/lasso_res%s_%s_%s.txt", oak_dir, resolution, new_population, pheno_name) 
    stat <- read_delim(stat_file, delim = " ", col_types = cols())
    print(stat)

    if(resolution !=0){
      Groups <- lapply(1:22, read_partitions, resolution = resolution)
      Groups <- do.call(rbind, Groups)
    }

    if(resolution != 0){
      stat <- stat %>% 
        mutate(Knockoff = endsWith(SNP, ".k"), SNP = str_replace(SNP, ".k", "")) %>%
        left_join(Groups, by=c("SNP", "CHR"))  %>%
        group_by(CHR, Group, Knockoff) %>%
        summarise(SNP.lead=SNP[which.max(abs(Z))], Z=sum(abs(Z)), Size = n()) %>%
        group_by(CHR, Group) %>%
        summarise(SNP.lead=SNP.lead[1], W=sum(Z[!Knockoff])-sum(Z[Knockoff]), Size = sum(Size)) %>%
        ungroup() %>%
        arrange(desc(abs(W)))
    }else{
      stat <- stat %>% 
        mutate(Knockoff = endsWith(SNP, ".k"), SNP = str_replace(SNP, ".k", "")) %>%
        group_by(CHR, SNP) %>%
        summarise(W=sum(Z[!Knockoff])-sum(Z[Knockoff])) %>%
        ungroup() %>%
        arrange(desc(abs(W)))
    }
  }
  
  ## Load the side information (bmi, cvd, sbp, diabetes)
  bmi_file <- sprintf("%s/popstruct/analysis/stats/lasso_%s_british_%s_res%d.txt", oak_dir, "bmi", relatedness, resolution)
  bmi <- read_delim(bmi_file, delim = " ", col_types = cols()) %>%
    rename(BMI = W) %>%
    mutate(BMI = abs(BMI))

  sbp_file <- sprintf("%s/popstruct/analysis/stats/lasso_%s_british_%s_res%d.txt", oak_dir, "sbp", relatedness, resolution)
  sbp <- read_delim(sbp_file, delim = " ", col_types = cols()) %>%
    rename(SBP = W) %>%
    mutate(SBP = abs(SBP))
  
  cvd_file <- sprintf("%s/popstruct/analysis/stats/lasso_%s_british_%s_res%d.txt", oak_dir, "cvd", relatedness, resolution)
  cvd <- read_delim(cvd_file, delim = " ", col_types = cols()) %>%
    rename(CVD = W) %>%
    mutate(CVD = abs(CVD))

  diabetes_file <- sprintf("%s/popstruct/analysis/stats/lasso_%s_british_%s_res%d.txt", oak_dir, "diabetes", relatedness, resolution)
  diabetes <- read_delim(diabetes_file, delim = " ", col_types = cols()) %>%
    rename(DIABETES = W) %>%
    mutate(DIABETES = abs(DIABETES))

  ## Match W with the side information
  stat <- stat %>% 
      ## the bmi prior info
      left_join(bmi, by = c("CHR", "Group")) %>%
      mutate(bmi = ifelse(is.na(BMI), 0 ,BMI)) %>%
      ## the sbp prior info 
      left_join(sbp, by = c("CHR", "Group")) %>%
      mutate(sbp = ifelse(is.na(SBP), 0, SBP)) %>%
      ## the diabetes prior info
      left_join(diabetes, by = c("CHR", "Group")) %>%
      mutate(diabetes = ifelse(is.na(DIABETES), 0, DIABETES)) %>%
      ## the cvd prior info
      left_join(cvd, by = c("CHR", "Group")) %>%
      mutate(cvd = ifelse(is.na(CVD), 0, CVD)) %>%
      ## Select the relavant features
      select(CHR, Group, W, diabetes, bmi, sbp, cvd) 

      print(stat)
  
  ## Diagnostics
  if(diagnostic){
    order_u <- order(stat$U, decreasing = TRUE)
    print(stat[order_u[1:10],])
    pp <- ggplot(stat, aes(x = 1 : dim(stat)[1], y = W[order_u])) + 
      geom_bar(stat = "identity", show.legend = FALSE)+
      theme_bw()
    ggsave("../figs/ordering.pdf", plot = pp, width = 6, height = 5)
  }
  
  ############################
  ## Run adaptive knockoffs ##
  ############################
  res <- filter_gam(stat$W, stat[[pheno_name]], alpha = alpha, offset = offset, mute = TRUE, df = 2,
                                      reveal_prop = prop)
  output <- rbind(output, data.frame(rej = length(res$rejs[[1]]), method = "Adaptive Knockoffs"))
  #cat(sprintf("The number of SNPs discovered by adaptive knockoffs is %d.\n", length(res$rejs[[1]])))
  
  ##############################################
  ## Run adaptive knockoffs W/ multiple prios ##
  ##############################################
  res <- filter_gam(stat$W, cbind(stat$diabetes, stat$bmi, stat$sbp, stat$cvd), 
                    alpha = alpha, offset = offset, mute = TRUE, reveal_prop = prop, df = 2)
  output <- rbind(output, data.frame(rej = length(res$rejs[[1]]), method = "Adaptive Knockoffs (multiple)"))
  #cat(sprintf("The number of SNPs discovered by adaptive knockoffs is %d.\n", length(res$rejs[[1]])))

  ########################
  ## Run mean knockoffs ##
  ########################
  W1 <- sign(stat$W) * (0.9 * abs(stat$W) + 0.1 * abs(stat[[pheno_name]]))
  tau <- knockoff::knockoff.threshold(W1, fdr = alpha, offset = offset)
  res <- which(W1 >= tau)
  output <- rbind(output, data.frame(rej = length(res), method = "Invariant Knockoffs"))
  #cat(sprintf("The number of SNPs discovered by mean knockoffs is %d.\n", length(res)))

  ###########################
  ## Run vanilla knockoffs ##
  ###########################
  tau <- knockoff::knockoff.threshold(stat$W, fdr = alpha, offset = offset)
  res <- which(stat$W >= tau)
  output <- rbind(output, data.frame(rej = length(res), method = "Vanilla Knockoffs"))
  #cat(sprintf("The number of SNPs discovered by vanilla knockoffs is %d.\n", length(res)))

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
  #print(wstat)

  tau <- knockoff::knockoff.threshold(wstat$W, fdr = alpha, offset = offset)
  res <- which(wstat$W >= tau)
  output <- rbind(output, data.frame(rej = length(res), method = "Weighted Lasso"))
  
  #####################################
  ## Run adaptive weighted knockoffs ##
  #####################################
  wstat <- wstat %>% left_join(diabetes, by = c("CHR", "Group")) %>%
      mutate(DIABETES = ifelse(is.na(DIABETES), 0, DIABETES) ) %>%
      left_join(bmi, by = c("CHR", "Group")) %>%
      mutate(BMI = ifelse(is.na(BMI), 0 ,BMI)) %>%
      left_join(sbp, by = c("CHR", "Group")) %>%
      mutate(SBP = ifelse(is.na(SBP), 0, SBP)) %>%
      left_join(cvd, by = c("CHR", "Group")) %>%
      mutate(CVD = ifelse(is.na(CVD), 0, CVD)) 

  res <-  filter_gam(wstat$W, cbind(wstat$DIABETES, wstat$BMI, wstat$SBP, wstat$CVD), alpha = alpha, offset = offset, 
                      mute = TRUE, df = 2, reveal_prop = prop)
  output <- rbind(output, data.frame(rej = length(res$rejs[[1]]), method = "Adaptive Weighted Knockoffs"))

  ####################################################################
  ## Run knockoffs with weighted lasso statistics (multiple priors) ##
  ####################################################################
  new_population <- factor(population, population_list, new_population_list) %>% as.character()
  wstat_file <- sprintf("%s/transfer_gwas/gwas/stats/lasso_multiple_res%d_%s_%s.txt", 
                        oak_dir, resolution, new_population, pheno_name) 
  wstat <- read_delim(wstat_file, delim = " ", col_types = cols())

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
  output <- rbind(output, data.frame(rej = length(res), method = "Weighted Lasso (Multiple Prior)"))



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
#population_list <- c("whitenonbritish")
#new_population_list <- ("European")
phenotype <- "cvd"
alpha <- 0.1
offset <- 0
prop <- 0.5
relatedness <- "related"

all_res <- data.frame()
for (population in population_list){
#  for (phenotype in pheno_list){
    cat(sprintf("Computing results for %s in the %s population...", phenotype, population))
    res <- lapply(resolution:resolution, filter, pheno_name = phenotype, population = population, 
                  relatedness = relatedness, prop = prop, alpha = alpha, offset = offset)
    res <- do.call(rbind, res)
    all_res <- rbind(all_res, res)
    cat("done.\n")
#  }
}
#all_res
out_file <- sprintf("%s/transfer_gwas/gwas/results/%s_res%d_fdr%s_offset%d_%s_multi.txt", 
                    oak_dir, phenotype, resolution, alpha, offset, relatedness)
write_delim(all_res, out_file, delim = " ")

#all_res %>% kbl(booktabs = TRUE, "latex")





