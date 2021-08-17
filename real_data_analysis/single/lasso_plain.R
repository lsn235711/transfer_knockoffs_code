#!/usr/bin/env Rscript

## Arguments
args <- commandArgs(trailingOnly=TRUE)
pop <- args[1]
res <- as.integer(args[2])
pheno.name <- args[3]

## Load libraries
suppressMessages(library(tidyverse))
suppressMessages(library(bigstatsr))
suppressMessages(library(bigsnpr))

oak <- "/oak/stanford/groups/candes"
scratch <- "/scratch/groups/candes"

## Parameters
dfmax <- 1000
n.cores <- 10
sum_stat <- TRUE
family <- ifelse(pheno.name %in% c("diabetes", "cvd"), "binomial", "gaussian")


## Load list of ethnicities
in.file <- sprintf("%s/transfer_gwas/gwas/qc/samples_%s.txt", oak, pop)
Ethnicities <- read_delim(in.file, delim = " ", col_names = c("FID", "IID"))

## Load list of subjects
fam.file <- sprintf("%s/popstruct/analysis/knockoffs_merged/ukb_gen_merged_res%d.fam", oak, res)
Subjects <- read_delim(fam.file, delim=" ",
                       col_names = c("FID", "IID", "X1", "X2", "Sex", "X3"), col_types=cols()) %>%
              left_join(Ethnicities)
ind.row <- which(Subjects$FID %in% Ethnicities$FID)
cat(sprintf("%d samples are loaded for training.\n", length(ind.row)))

## Load data file
fbm.file <- sprintf("%s/transfer_gwas/gwas/fbm/ukb_gen_merged_res%d_%s", scratch, res, pop)
rds.file <- sprintf("%s.rds", fbm.file)

if(file.exists(rds.file)){
  cat(sprintf("Found FBM in %s.\n", rds.file))
} else {
  cat(sprintf("Could not find FBM in %s.\n", rds.file))
  cat(sprintf("Converting BED data to FBM... "))
  ## Load data file
  bed.file <- sprintf("%s/popstruct/analysis/knockoffs_merged/ukb_gen_merged_res%d.bed", oak, res)
  rds.file <- snp_readBed2(bed.file, backingfile = fbm.file, ind.row = ind.row)
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
pheno.file <- "/oak/stanford/groups/candes/popstruct/data/phenotypes/phenotypes.tab"
Phenotypes <- read_delim(pheno.file, delim="\t", col_names = TRUE) %>% mutate(IID=FID)
Phenotypes <- FAM %>% inner_join(Phenotypes, by = c("FID", "IID"))

## Extract phenotypes
cat(sprintf("- Extracting phenotype %s... ", pheno.name))
y <- Phenotypes[[pheno.name]]
id.notna <- which(!is.na(y))
y <- y[id.notna]
cat(sprintf("There are %d responses with non-missing value.\n", length(y)))
cat("OK.\n")

## Compute scaling factor for the genotypes
#pop.list <- c("African", "Asian", "European", "Indian")
#new.pop.list <- c("black", "asian", "whitenonbritish", "indian")
#new.pop <- factor(pop, levels = pop.list, labels = new.pop.list)
#new.pop <- as.character(new.pop)
scale.file <- sprintf("%s/transfer_gwas/gwas/scale/scale_%s_related_res%d.txt", oak, pop, res)
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

## Load partition information
read_partitions <- function(chr){
  in.file <- sprintf("/oak/stanford/groups/candes/popstruct/data/partitions/hap_chr%d.txt", chr)
  Partitions <- read_delim(in.file, delim=" ")
  res.name <- sprintf("res_%d", res)
  groups <- Partitions %>% mutate(CHR = chr) %>%
  rename(Group = res.name) %>% 
  select(CHR, SNP, Group)
  return(groups)
}
if(res !=0){
  Groups <- lapply(1:22, read_partitions)
  Groups <- do.call(rbind, Groups)
}
## print(Groups)

## Load covariates
cat("- Loading covariate data... ")
covariate.file <- sprintf("%s/popstruct/data/phenotypes/covariates.tab", oak)
Analysis <- read_tsv(covariate.file, col_types=cols(), progress=FALSE)
covariate.names <- Analysis %>% filter(Name==pheno.name) %>% select(Covariates) %>%
    as.character() %>% strsplit(",")
covariate.names <- c(covariate.names[[1]], paste("PC", seq(5), sep="."))
Covariates <- Phenotypes %>% select(covariate.names)
cat("OK.\n")

## Determine which samples to use
ind.train <- 1:nrow(G)
ind.train <- ind.train[id.notna]
print(length(ind.train))
covar.train <- as.matrix(Covariates)[ind.train,]

## Find the class of response (numeric or factor)
cat("- Checking phenotype family: ")
if(family=="binomial") {
  y <- factor(y, levels=c(1,2), labels=c(0,1))
  y <- as.numeric(levels(y))[y]
}
cat(sprintf("%s.\n", family))


## Fit lasso
if(family == "binomial"){
  cat(sprintf("- Fitting lasso on %d x (%d + %d) data... ", length(ind.train), ncol(G), ncol(covar.train)))
  start_time <- Sys.time()
  lasso.fit <- big_spLogReg(G, y01.train=y, dfmax=dfmax, ncores=n.cores, 
                           ind.train=ind.train, covar.train = covar.train)
  end_time <- Sys.time()
  cat("OK.\n")
  print(end_time - start_time)
}else{
  cat(sprintf("- Fitting lasso on %d x (%d + %d) data... ", length(ind.train), ncol(G), ncol(covar.train)))
  start_time <- Sys.time()
  lasso.fit <- big_spLinReg(G, y.train=y, dfmax=dfmax, ncores=n.cores, 
                           ind.train=ind.train, covar.train = covar.train)
  end_time <- Sys.time()
  cat("OK.\n")
  print(end_time - start_time)
}

## Extracting lasso coefficients
cat("- Extracting model coefficients... ")
betas <- sapply(1:10, function(k) lasso.fit[[1]][k][[1]]$beta)
beta.extracted <- get_beta(betas, method = c("mean"))
beta.nonzero <- beta.extracted[1:(length(beta.extracted) - ncol(covar.train))]
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
out.file <- sprintf("%s/transfer_gwas/gwas/stats/lasso_res%d_%s_%s.txt",
                    oak, res, pop, pheno.name) 
Beta %>% write_delim(out.file, delim=" ")
cat(sprintf("Lasso importance measures written to:\n  %s\n", out.file))
