#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

pop <- "British"
n.causal <- 100
specificity <- 0.5
snr <- 0.1
fold <- 1
transfer <- 0

pop <- args[1]
n.causal <- as.integer(args[2])
specificity <- as.integer(args[3])/100
snr <- as.integer(args[4])/100
fold <- as.integer(args[5])
transfer <- as.integer(args[6])

## Load libraries
suppressMessages(library(tidyverse))
suppressMessages(library(bigstatsr))
suppressMessages(library(bigsnpr))

pop.label <- pop
if(length(str_split(pop, "-small")[[1]])==2) {
  pop = str_split(pop, "-small")[[1]][1]
}


oak <- "/oak/stanford/groups/candes/matteo/transfer_knockoffs"
scratch <- "/scratch/groups/candes/matteo/transfer_knockoffs"

## Load list of ethnicities
in.file <- sprintf("%s/data/qc/ethnicity.txt", oak)
Ethnicities <- read_delim(in.file, delim=" ") %>% mutate(FID=ID, IID=ID) %>% select(FID, IID, Group)

## Load list of subjects
fam.file <- sprintf("%s/data/knockoffs/ukb_gen_chr1_res5_%s.fam", oak, pop)
Subjects <- read_delim(fam.file, delim=" ",
                       col_names = c("FID", "IID", "X1", "X2", "Sex", "X3"), col_types=cols()) %>%
  left_join(Ethnicities)

## Load list of causal variants
in.file <- sprintf("%s/data/causal/causal_n%s_s%s.txt", oak, n.causal, round(100*specificity))
Causal <- read_delim(in.file, delim=" ") %>%
  mutate(Causal = TRUE) %>%
  select(CHR, Group, BP, SNP, Causal)

## Load importance measures
if(transfer) {
  stats.file <- sprintf("%s/stats/lasso_transfer_res5_n%s_snr%s_s%s_%s_fold%d.txt",
                        oak, n.causal, round(100*snr), round(100*specificity), pop.label, fold)
} else {
  stats.file <- sprintf("%s/stats/lasso_res5_n%s_snr%s_s%s_%s_fold%d.txt",
                        oak, n.causal, round(100*snr), round(100*specificity), pop.label, fold)
}
Beta <- read_delim(stats.file, delim=" ")

## Load partition information
in.file <- "/oak/stanford/groups/candes/popstruct/data/partitions/hap_chr1.txt"
Partitions <- read_delim(in.file, delim=" ")
Groups <- Partitions %>% mutate(Group = `res_5`) %>% select(SNP, Group)

## Compute knockoff statistics
Stats <- Beta %>%
  mutate(Knockoff=endsWith(SNP, ".k"), SNP=str_replace(SNP, ".k", "")) %>%
  left_join(Groups, by="SNP") %>%
  group_by(CHR, Group, Knockoff) %>%
  summarise(SNP.lead=SNP[which.max(abs(Z))], Z=sum(abs(Z))) %>%
  group_by(CHR, Group) %>%
  summarise(SNP.lead=SNP.lead[1], W=sum(Z[!Knockoff])-sum(Z[Knockoff])) %>%
  ungroup() %>%
  arrange(desc(abs(W)))

## Apply knockoff filter
source("utils.R")
Discoveries <- Stats %>%
  mutate(Threshold = knockoff.threshold(W, fdr=0.1, offset=1)) %>%
  filter(W>Threshold)
  
## Evaluate discoveries
df <- Causal %>%
  full_join(Discoveries) %>%
  mutate(Causal = ifelse(is.na(Causal), FALSE, Causal)) %>%
  mutate(Discovered = ifelse(is.na(W), FALSE, TRUE))

res.fdp <- df %>%
  filter(Discovered) %>%  
  summarise(FDP = mean(!Causal))
res.pow <- df %>%
  filter(Causal) %>%  
  summarise(Power = mean(Discovered))
cat(sprintf("Power: %.2f, FDP: %.2f\n", res.pow, res.fdp))

## Save list of discoveries
if(transfer) {
  out.file <- sprintf("%s/discoveries/discoveries_transfer_res5_n%s_snr%s_s%s_%s_fold%d.txt",
                      oak, n.causal, round(100*snr), round(100*specificity), pop.label, fold)
} else {
  out.file <- sprintf("%s/discoveries/discoveries_res5_n%s_snr%s_s%s_%s_fold%d.txt",
                      oak, n.causal, round(100*snr), round(100*specificity), pop.label, fold)    
}
df %>% write_delim(out.file, delim=" ")
cat(sprintf("List of discoveries written to:\n  %s\n", out.file))
