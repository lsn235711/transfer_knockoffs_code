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

## Load prior information from pooled samples
if(sum_stat){
  prior.file <- sprintf("%s/popstruct/analysis/importance/lasso_%s_british_related_res%d.txt",
                        oak, pheno.name, res)
  if(res != 0){
    Beta.prior <- read_delim(prior.file, delim=" ", col_types=cols()) %>%
      mutate(Knockoff=endsWith(SNP, ".k"), SNP=str_replace(SNP, ".k", "")) %>%
      full_join(Groups, by=c("CHR", "SNP")) %>%
      mutate(Z = ifelse(is.na(Z), 0, Z), Beta = ifelse(is.na(Beta), 0, Beta)) %>%
      group_by(CHR, Group, SNP) %>%
      summarise(Z = sum(abs(Z)))

    ## Compute prior weights
    Prior <- Beta.prior %>%
      group_by(CHR, Group, SNP) %>%
      mutate(Z = sum(abs(Z)), Pi = 1/(0.05+Z)) %>%
      ungroup() %>%
      mutate(Pi = Pi / max(Pi))

    ## Re-order prior
    Prior <- BIM %>%
      mutate(Knockoff=endsWith(SNP, ".k"), SNP=str_replace(SNP, ".k", "")) %>%
      filter(!Knockoff) %>%
      left_join(Prior, by=c("CHR", "SNP"))
  }else{
    
    Beta.prior <- read_delim(prior.file, delim=" ", col_types=cols()) %>%
      mutate(Knockoff=endsWith(SNP, ".k"), SNP=str_replace(SNP, ".k", "")) %>%
      group_by(CHR, SNP) %>%
      summarise(Z = sum(abs(Z)))

    ## Compute prior weights
    #Prior <- Beta.prior %>%
    #  group_by(CHR, SNP) %>%
    #  mutate(Z = sum(abs(Z)), Pi = 1/(0.05+Z)) %>%
    #  ungroup() %>%
    #  mutate(Pi = Pi / max(Pi))
    #print(Prior)

    ## Re-order prior
    Prior <- BIM %>%
      mutate(Knockoff=endsWith(SNP, ".k"), SNP=str_replace(SNP, ".k", "")) %>%
      filter(!Knockoff) %>%
      left_join(Beta.prior, by=c("CHR", "SNP")) %>%
      mutate(Z = ifelse(is.na(Z), 0, Z)) %>%
      mutate(Pi = 1 / (0.05 + Z)) %>%
      mutate(Pi = Pi / max(Pi))
  }
    pi.values <- Prior$Pi

}else{
  prior.file <- sprintf("%s/popstruct/analysis/stats/lasso_%s_everyone_related_res%d.txt",
                      oak, pheno.name, res)

  ## Compute prior weights
  Prior <- read_delim(prior.file, delim = " ", col_types = cols()) 
  if(res != 0){
    Prior <- Prior %>% right_join(Groups, by = c("CHR", "Group")) %>%
            mutate(W = ifelse(is.na(W), 0, W), Pi = 1 / (0.05 + abs(W))) %>%
            mutate(Pi = Pi / max(Pi))
    ## Re-order prior
    Prior <- BIM %>%
      mutate(Knockoff=endsWith(SNP, ".k"), SNP=str_replace(SNP, ".k", "")) %>%
      filter(!Knockoff) %>%
      left_join(Prior, by=c("CHR", "SNP"))
  
    pi.values <- Prior$Pi
  }else{
  
    Prior <- Prior %>% rename(SNP = SNP.lead)
  
    Prior <-BIM %>%
          mutate(Knockoff=endsWith(SNP, ".k"), SNP=str_replace(SNP, ".k", "")) %>%
          filter(!Knockoff) %>%
          left_join(Prior, by=c("CHR", "SNP")) %>%
          mutate(W = ifelse(is.na(W), 0 ,W)) %>%
          mutate(Pi = 1 / (0.05 + abs(W))) %>%
          mutate(Pi = Pi / max(Pi))
  
          pi.values <- Prior$Pi

  }
}
#Prior <- Beta.prior %>%
#  group_by(CHR, Group, SNP) %>%
#  mutate(Z = sum(abs(Z)), Pi = 1/(0.05+Z)) %>%
#  ungroup() %>%
#  mutate(Pi = Pi / max(Pi))

cat(sprintf("The prior is of length %d.\n", length(pi.values)))

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

## Outer hyperparameter
gamma.seq <- seq(0, 1, length.out=10)

ind.sets <- sample(rep(1:10, length.out=length(ind.train)))

fit_lasso <- function(gamma) {
  pf.X <- rep((1-gamma) + gamma * pi.values, each=2)
  cat(sprintf("- The length of the response is %d; the length of the prior information is %d; the length of
              the indsets is %d.\n", length(y), length(pf.X), length(ind.sets)))
  cat(sprintf("- Fitting lasso on %d x (%d + %d) data with gamma = %.2f...\n", length(ind.train), ncol(G), ncol(covar.train), gamma))
  start_time <- Sys.time()
  if(family == "binomial"){
    lasso.fit <- big_spLogReg(G, y01.train=y, dfmax=dfmax, ncores=n.cores, pf.X=pf.X, 
                            ind.sets=ind.sets, ind.train=ind.train, covar.train = covar.train)
  }else{
    lasso.fit <- big_spLinReg(G, y.train=y, dfmax=dfmax, ncores=n.cores, pf.X=pf.X, 
                            ind.sets=ind.sets, ind.train=ind.train, covar.train = covar.train)
  }
  end_time <- Sys.time()
  cat("OK.\n")
  print(end_time - start_time)

  ## Extract validation loss
  #loss <- summary(lasso.fit, best.only = TRUE)$validation_loss
  lambda.seq <- do.call("c", lapply(1:10, function(k) lasso.fit[[1]][[k]]$lambda))
  lambda.seq <- seq(min(lambda.seq), max(lambda.seq), length.out=100)
  Losses <- do.call("cbind", lapply(1:10, function(k) {
    spl <- smooth.spline(lasso.fit[[1]][[k]]$lambda, lasso.fit[[1]][[k]]$loss.val)
    predict(spl, lambda.seq)$y
  }))
  loss.avg.seq <- apply(Losses, 1, "mean")
  loss.sd.seq <- apply(Losses, 1, "sd")
  idx.min <- which.min(loss.avg.seq)
  loss <- loss.avg.seq[idx.min]
  loss.sd <- loss.sd.seq[idx.min]
  cat(sprintf("Loss = %.2f (%.2f, %.2f)\n", loss, loss-loss.sd, loss+loss.sd))
  return(c(loss, loss.sd))
}

loss.seq <- sapply(gamma.seq, function(gamma) fit_lasso(gamma))

if(FALSE) {

  p1 <- tibble(gamma=gamma.seq, loss=loss.seq[1,], loss.upp=loss.seq[1,]+loss.seq[2,], loss.low=loss.seq[1,]-loss.seq[2,]) %>%
    gather(loss, loss.upp, loss.low, key="key", value="value") %>%
    ggplot(aes(x=gamma, y=value, color=key)) +
    geom_point() +
    geom_line() +
    theme_bw()
  p1

  tibble(gamma=gamma.seq, loss=loss.seq[1,], loss.upp=loss.seq[1,]+loss.seq[2,], loss.low=loss.seq[1,]-loss.seq[2,]) %>%
    gather(loss, loss.upp, loss.low, key="key", value="value") %>%
    ggplot(aes(x=gamma, y=value, color=key)) +
    geom_point() +
    geom_line() +
    theme_bw()

}

## Find the best gamma value
idx.best <- which.min(loss.seq[1,])
gamma.best <- gamma.seq[idx.best]

## Fit lasso
pf.X <- rep((1-gamma.best) + gamma.best * pi.values, each=2)
cat(sprintf("- Fitting lasso on %d x (%d + %d) data with gamma = %.2f... ", length(ind.train), ncol(G), ncol(covar.train), gamma.best))
start_time <- Sys.time()
if(family == "binomial"){
  lasso.fit <- big_spLogReg(G, y01.train=y, dfmax=dfmax, ncores=n.cores, pf.X=pf.X, ind.sets=ind.sets, ind.train=ind.train, covar.train = covar.train)
}else{
  lasso.fit <- big_spLinReg(G, y.train=y, dfmax=dfmax, ncores=n.cores, pf.X=pf.X, ind.sets=ind.sets, ind.train=ind.train, covar.train = covar.train)
}
end_time <- Sys.time()
cat("OK.\n")
print(end_time - start_time)

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
Beta <- BIM %>% mutate(Beta = beta.variants, Gamma = gamma.best,
                       Scale=scaling.factors,
                       Z=Beta*Scale) %>%
  filter(Z!=0) %>%
  arrange(desc(abs(Z)))
cat("OK.\n")

## Save importance measures
out.file <- sprintf("%s/transfer_gwas/gwas/stats/lasso_transfer_res%d_%s_%s.txt",
                    oak, res, pop, pheno.name) 
Beta %>% write_delim(out.file, delim=" ")
cat(sprintf("Lasso importance measures written to:\n  %s\n", out.file))
