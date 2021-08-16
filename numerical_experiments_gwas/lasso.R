#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

pop <- "Everyone"
n.causal <- 100
specificity <- 0.5
snr <- 0.1
fold <- 1

pop <- args[1]
n.causal <- as.integer(args[2])
specificity <- as.integer(args[3])/100
snr <- as.integer(args[4])/100
fold <- as.integer(args[5])

## Load libraries
suppressMessages(library(tidyverse))
suppressMessages(library(bigstatsr))
suppressMessages(library(bigsnpr))

pop.label <- pop
if(length(str_split(pop, "-small")[[1]])==2) {
  pop = str_split(pop, "-small")[[1]][1]
  full.sample <- FALSE
} else {
  full.sample <- TRUE
}

oak <- "/oak/stanford/groups/candes/matteo/transfer_knockoffs"
scratch <- "/scratch/groups/candes/matteo/transfer_knockoffs"

## Parameters
dfmax <- 1000
n.cores <- 5

## Load list of ethnicities
in.file <- sprintf("%s/data/qc/ethnicity.txt", oak)
Ethnicities <- read_delim(in.file, delim=" ") %>% mutate(FID=ID, IID=ID) %>% select(FID, IID, Group)

## Load list of subjects
fam.file <- sprintf("%s/data/knockoffs/ukb_gen_chr1_res5_%s.fam", oak, pop)
Subjects <- read_delim(fam.file, delim=" ",
                       col_names = c("FID", "IID", "X1", "X2", "Sex", "X3"), col_types=cols()) %>%
  left_join(Ethnicities)

## Load data file
fbm.file <- sprintf("%s/fbm/ukb_gen_chr1_res5_%s", scratch, pop)
rds.file <- sprintf("%s.rds", fbm.file)
if(file.exists(rds.file)){
  cat(sprintf("Found FBM in %s.\n", rds.file))
} else {
  cat(sprintf("Could not find FBM in %s.\n", rds.file))
  cat(sprintf("Converting BED data to FBM... "))
  ## Load data file
  if(pop=="Everyone") {
    ind.row <- 1:nrow(Subjects)
  } else {
    ind.row <- which(Subjects$Group == pop)
  }
  bed.file <- sprintf("%s/data/knockoffs/ukb_gen_chr1_res5_%s.bed", oak, pop)
  rds.file <- snp_readBed2(bed.file, backingfile=fbm.file, ind.row=ind.row)    
  cat("OK\n")
}
obj.bigSNP <- snp_attach(rds.file)

## Get alias for genotype data
G <- obj.bigSNP$genotypes
cat(sprintf("Loaded genotypes: %d by %d.\n", nrow(G), ncol(G)))

## Get list of variants
BIM <- obj.bigSNP$map %>% as_tibble() %>% transmute(CHR=chromosome, SNP=marker.ID, BP=physical.pos)

## Get list of samples
FAM <- obj.bigSNP$fam %>% as_tibble() %>% transmute(FID=family.ID, IID=sample.ID)

## Load phenotypes
pheno.dir <- "/oak/stanford/groups/candes/matteo/transfer_knockoffs/data/phenotypes"
if(pop=="Everyone") {
  pop.list <- c("African", "Asian", "British", "European", "Indian")
  Phenotypes <- do.call("rbind", lapply(pop.list, function(pop.tmp) {
    pheno.file <- sprintf("%s/phenotypes_n%s_snr%s_s%s_%s.txt", pheno.dir, n.causal, round(100*snr), round(100*specificity), pop.tmp)
    tmp <- read_delim(pheno.file, delim=" ") %>% mutate(IID=FID)
    return(tmp)
  }))
} else {
  pheno.file <- sprintf("%s/phenotypes_n%s_snr%s_s%s_%s.txt", pheno.dir, n.causal, round(100*snr), round(100*specificity), pop)
  Phenotypes <- read_delim(pheno.file, delim=" ") %>% mutate(IID=FID)
}
Phenotypes <- FAM %>% inner_join(Phenotypes, by = c("FID", "IID"))

## Extract phenotypes
pheno.name <- sprintf("Y%d", fold)
cat(sprintf("- Extracting phenotype %s... ", pheno.name))
y <- Phenotypes[[pheno.name]]
cat("OK.\n")

## Compute scaling factor for the genotypes
scale.file <- sprintf("%s/data/knockoffs/scale_chr1_res5_%s.txt", oak, pop)
if(file.exists(scale.file)) {
  cat("- Loading scaling factors for all variants... ")
  scaling.factors <- read_delim(scale.file, delim=" ", col_types=cols(),progress=FALSE)$Scale
  cat("OK.\n")
} else {
  ## Compute scaling factor for the genotypes
  cat("- Computing scaling factors for all variants... ")
  scaler <- big_scale()
  G.scale <- scaler(G)
  scaling.factors <- G.scale$scale
  cat("OK.\n")
  ## Save scaling factors
  BIM %>% mutate(Scale=scaling.factors) %>% write_delim(scale.file, delim=" ")
  cat("- Saved scaling factors for all variants.\n")
}

## Determine which samples to use
set.seed(2021)
if(full.sample) {
  ind.train <- 1:nrow(G)
} else {
  ind.train <- sort(sample(nrow(G), 3284))
}

## Fit lasso
cat(sprintf("- Fitting lasso on (%d x %d) data... ", nrow(G), ncol(G)))
start_time <- Sys.time()
lasso.fit <- big_spLinReg(G, y.train=y[ind.train], dfmax=dfmax, ncores=n.cores, ind.train=ind.train)
end_time <- Sys.time()
cat("OK.\n")
print(end_time - start_time)

## Extracting lasso coefficients
cat("- Extracting model coefficients... ")
betas <- sapply(1:10, function(k) lasso.fit[[1]][k][[1]]$beta)
beta.extracted <- get_beta(betas, method = c("mean"))
beta.nonzero <- beta.extracted[1:(length(beta.extracted))]
lasso.indices <- attr(lasso.fit, "ind.col")
beta.variants <- rep(0, ncol(G))
beta.variants[lasso.indices] <- beta.nonzero
model.support <- which(beta.variants!=0)
cat("OK.\n")
cat(sprintf("- Model support size = %d.\n", length(model.support)))

## Compute importance measures
cat("- Computing importance measures... ")
Beta <- BIM %>% mutate(Beta = beta.variants,
                       Scale=scaling.factors,
                       Z=Beta*Scale) %>%
  filter(Z!=0) %>%
  arrange(desc(abs(Z)))
cat("OK.\n")

## Save importance measures
out.file <- sprintf("%s/stats/lasso_res5_n%s_snr%s_s%s_%s_fold%d.txt",
                    oak, n.causal, round(100*snr), round(100*specificity), pop.label, fold)
Beta %>% write_delim(out.file, delim=" ")
cat(sprintf("Lasso importance measures written to:\n  %s\n", out.file))
