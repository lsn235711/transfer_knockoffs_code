#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

pop <- "African"
n.causal <- 100
specificity <- 1
traget.snr <- 0.15

pop <- args[1]
n.causal <- as.integer(args[2])
specificity <- as.integer(args[3])/100
target.snr <- as.integer(args[4])/100

suppressMessages(library(tidyverse))
suppressMessages(library(bigstatsr))
suppressMessages(library(bigsnpr))

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

ind.row <- which(Subjects$Group == pop)

## Load list of variants
bim.file <- sprintf("%s/data/knockoffs/ukb_gen_chr1_res5_%s.bim", oak, pop)
Variants <- read_delim(bim.file, delim="\t",
                       col_names=c("CHR", "SNP", "cM", "BP", "A1", "A2"), col_types=cols()) %>%
  filter(!endsWith(SNP, ".k")) %>%
  select(CHR, SNP, BP)

## Load list of causal variants
in.file <- sprintf("%s/data/causal/causal_n%s_s%s.txt", oak, n.causal, round(100*specificity))
Causal <- read_delim(in.file, delim=" ")

## Load data file
fbm.file <- sprintf("%s/fbm/ukb_gen_chr1_res5_%s", scratch, pop)
rds.file <- sprintf("%s.rds", fbm.file)
if(file.exists(rds.file)){
  cat(sprintf("Found FBM in %s.\n", rds.file))
} else {
  cat(sprintf("Could not find FBM in %s.\n", rds.file))
  cat(sprintf("Converting BED data to FBM... "))
  ## Load data file
  bed.file <- sprintf("%s/data/knockoffs/ukb_gen_chr1_res5_%s.bed", oak, pop)
  rds.file <- snp_readBed2(bed.file, backingfile=fbm.file, ind.row=ind.row)    
  cat("OK\n")
}
obj.bigSNP <- snp_attach(rds.file)

## Extract causal variants
idx.causal <- which(obj.bigSNP$map$marker.ID %in% Causal$SNP)
X <- obj.bigSNP$genotypes[,idx.causal]

## Define effects
gamma <- 10
set.seed(2021)
beta.signs <- (2*rbinom(n.causal, 1, 0.5)-1)
set.seed(digest::digest2int(pop))
beta.magn <- runif(n.causal, min=1/gamma, max=gamma)
beta <- beta.signs * beta.magn

## Compute vector of means
y.means = X %*% beta

## Choose noise level
noise.sd <- sqrt(var(y.means)/target.snr)

## Generate phenotypes
Y = do.call("cbind", lapply(1:10, function(seed) {
  set.seed(seed)
  y = y.means + noise.sd * rnorm(length(y.means))
}))

snr <- var(y.means) / (var(Y[,1])-var(y.means))
cat(sprintf("SNR = %.2f%%\n", snr*100))

phenotypes <- as_tibble(Y)
colnames(phenotypes) <- paste("Y", 1:ncol(Y), sep="")
phenotypes$FID <- Subjects$FID[ind.row]
phenotypes <- phenotypes %>% select(FID, everything())
  
## Store phenotypes
out.dir <- "/oak/stanford/groups/candes/matteo/transfer_knockoffs/data/phenotypes"
out.file <- sprintf("%s/phenotypes_n%s_snr%s_s%s_%s.txt", out.dir, n.causal, round(100*target.snr), round(100*specificity), pop)
phenotypes %>% write_delim(out.file, delim=" ")

cat(sprintf("Phenotypes written to: %s\n", out.file))
