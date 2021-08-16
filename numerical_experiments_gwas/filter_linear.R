#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

pop <- "Indian"
n.causal <- 100
specificity <- 0
snr <- 0.05
fold <- 1

pop <- args[1]
n.causal <- as.integer(args[2])
specificity <- as.integer(args[3])/100
snr <- as.integer(args[4])/100
fold <- as.integer(args[5])

## Load libraries
suppressMessages(library(tidyverse))
source("utils.R")

pop.label <- pop
if(length(str_split(pop, "-small")[[1]])==2) {
  pop = str_split(pop, "-small")[[1]][1]
}

oak <- "/oak/stanford/groups/candes/matteo/transfer_knockoffs"
scratch <- "/scratch/groups/candes/matteo/transfer_knockoffs"

## Load partition information
in.file <- "/oak/stanford/groups/candes/popstruct/data/partitions/hap_chr1.txt"
Partitions <- read_delim(in.file, delim=" ")
Groups <- Partitions %>% mutate(Group = `res_5`) %>% select(SNP, Group)

## Load prior information from British samples
prior.file <- sprintf("%s/stats/lasso_res5_n%s_snr%s_s%s_%s_fold%d.txt", oak, n.causal, round(100*snr), round(100*specificity), "Everyone", fold)
Beta.prior <- read_delim(prior.file, delim=" ", col_types=cols()) %>%
  mutate(Knockoff=endsWith(SNP, ".k"), SNP=str_replace(SNP, ".k", "")) %>%
  full_join(Groups) %>%
  mutate(Z = ifelse(is.na(Z), 0, Z), Beta = ifelse(is.na(Beta), 0, Beta)) %>%
  group_by(CHR, Group, SNP) %>%
  summarise(Z = sum(abs(Z)))

## Compute knockoff statistics
Stats.prior <- Beta.prior %>%
  mutate(Knockoff=endsWith(SNP, ".k"), SNP=str_replace(SNP, ".k", "")) %>%
  left_join(Groups) %>%
  group_by(CHR, Group, Knockoff) %>%
  summarise(SNP.lead=SNP[which.max(abs(Z))], Z=sum(abs(Z))) %>%
  group_by(CHR, Group) %>%
  summarise(SNP.lead=SNP.lead[1], W.prior=sum(Z[!Knockoff])-sum(Z[Knockoff])) %>%
  ungroup() %>%
  arrange(desc(abs(W.prior))) %>%
  select(-SNP.lead)

## Load importance measures
stats.file <- sprintf("%s/stats/lasso_res5_n%s_snr%s_s%s_%s_fold%d.txt", oak, n.causal, round(100*snr), round(100*specificity), pop.label, fold)
Beta <- read_delim(stats.file, delim=" ")

## Compute knockoff statistics
Stats.new <- Beta %>%
  mutate(Knockoff=endsWith(SNP, ".k"), SNP=str_replace(SNP, ".k", "")) %>%
  left_join(Groups, by="SNP") %>%
  group_by(CHR, Group, Knockoff) %>%
  summarise(SNP.lead=SNP[which.max(abs(Z))], Z=sum(abs(Z))) %>%
  group_by(CHR, Group) %>%
  summarise(SNP=SNP.lead[1], W.new=sum(Z[!Knockoff])-sum(Z[Knockoff])) %>%
  ungroup() %>%
  arrange(desc(abs(W.new)))

## Load list of causal variants
in.file <- sprintf("%s/data/causal/causal_n%s_s%s.txt", oak, n.causal, round(100*specificity))
Causal <- read_delim(in.file, delim=" ") %>%
  mutate(Causal = TRUE) %>%
  mutate(SNP.c = SNP) %>%
  select(CHR, Group, SNP.c, Causal)

##
theta.list <- seq(0.1,0.9,by=0.1)

Results <- do.call("rbind", lapply(theta.list, function(theta) {
  ## Combine stats
  Stats <- left_join(Stats.new, Stats.prior)%>%
    mutate(W = (1-theta)*abs(W.new) + theta*abs(W.prior), W = W*sign(W.new),
           W = ifelse(is.na(W), 0, W)) %>%
    arrange(desc(abs(W)))
  ## Apply knockoff filter
  Discoveries <- Stats %>%
    mutate(Threshold = knockoff.threshold(W, fdr=0.1, offset=1)) %>%
    filter(W>Threshold) 
  ## Evaluate discoveries
  df <- Causal %>%
    full_join(Discoveries) %>%
    mutate(Causal = ifelse(is.na(Causal), FALSE, Causal)) %>%
    mutate(Discovered = ifelse(is.na(W), FALSE, TRUE)) %>%
    mutate(Theta = theta)
  res.fdp <- df %>%
    filter(Discovered) %>%  
    summarise(FDP = mean(!Causal)) %>%
    mutate(FDP = ifelse(is.nan(FDP), 0, FDP))
  res.pow <- df %>%
    filter(Causal) %>%  
    summarise(Power = mean(Discovered))
  cat(sprintf("Theta = %.2f, Power: %.2f, FDP: %.2f\n", theta, res.pow, res.fdp))
  return(df)
}))

## Save list of discoveries
out.file <- sprintf("%s/discoveries_others/discoveries_transfer_2_res5_n%s_snr%s_s%s_%s_fold%d.txt",
                    oak, n.causal, round(100*snr), round(100*specificity), pop.label, fold)
Results %>% write_delim(out.file, delim=" ")
cat(sprintf("List of discoveries written to:\n  %s\n", out.file))
