library(tidyverse)

resolution <- 5

## Load partition information
in.file <- "/oak/stanford/groups/candes/popstruct/data/partitions/hap_chr1.txt"
Partitions <- read_delim(in.file, delim=" ")

Groups <- Partitions %>% mutate(Group = `res_5`) %>% select(SNP, Group)

pop.list <- c("African", "Asian", "Indian", "European", "British")
n.pop <- length(pop.list)
in.dir <- "/oak/stanford/groups/candes/matteo/transfer_knockoffs/data/knockoffs"
Variants.raw <- do.call("rbind", lapply(pop.list, function(pop) {
  in.file <- sprintf("%s/ukb_gen_chr1_res%s_%s.frq", in.dir, resolution, pop)
  df.1 <- read_table(in.file, col_types=cols()) %>% mutate(Population = pop)
  in.file <- sprintf("%s/ukb_gen_chr1_res%s_%s.bim", in.dir, resolution, pop)
  df.2 <- read_table2(in.file, col_types=cols(), col_names=c("CHR", "SNP", "X", "BP", "A1", "A2"))
  df.1 %>% left_join(select(df.2, SNP, BP)) %>%
    select(CHR, BP, SNP, everything())
}))

Variants <- Variants.raw %>%
  separate(SNP, into=c("SNP", "Knockoff"), "[.]") %>%
  mutate(Knockoff = ifelse(is.na(Knockoff), FALSE, TRUE)) %>%
  filter(!Knockoff) %>%
  select(-Knockoff, -NCHROBS, -A1, -A2) %>%
  arrange(CHR, BP, Population) %>%
  left_join(Groups, by="SNP") %>%
  select(CHR, Group, BP, SNP, everything())

Variants %>% group_by(Population) %>% summarise(N=n())

df <- Variants %>%
  arrange(CHR, Group, BP, SNP, desc(MAF)) %>%
  group_by(CHR, Group, BP, SNP) %>%
  mutate(idx = which.max(MAF), Population.max = Population[idx], MAF.max=MAF[idx],
         MAF.second = -sort(-MAF)[2],
         idx = which.min(MAF), Population.min = Population[idx], MAF.min=MAF[idx]) %>%
  select(-idx) %>%
  mutate(Hom = MAF.min/MAF.max, Het = MAF.max / MAF.second) %>%
  ungroup()

## Store summarised MAF information
MAF <- df %>%
  group_by(CHR, Group, BP, SNP, Population.max, MAF.max, Hom, Het) %>%
  summarise()
out.dir <- "/oak/stanford/groups/candes/matteo/transfer_knockoffs/data/qc"
out.file <- sprintf("%s/maf_info.txt", out.dir)
MAF %>% write_delim(out.file, delim=" ")

## Total number of causal variants
n.causal <- 100

## How many variants should be population-specific?
specificity <- 0.5

select_causal <- function(n.causal, specificity) {
  n.het <- round(n.causal * specificity)
  
  ## Make list of top homogeneous and common variants
  Variants.hom <- df %>%
    filter(MAF.min>0.1) %>%
    select(-Population.max, -MAF.min, -MAF.max, -MAF.second) %>%
    group_by(CHR, Group, BP, SNP, Hom, Het) %>%
    summarise(MAF=max(MAF)) %>%
    arrange(desc(Hom)) %>%
    head(2*n.causal)

  ## Sort variants from the most in decreasing order of specifificy
  Variants.het <- df %>%
    filter(MAF.max>0.1) %>%
    select(-MAF.min, -MAF.max, -MAF.second) %>%
    mutate(Population = Population.max) %>%
    group_by(CHR, Group, BP, SNP, Population, Hom, Het) %>%
    summarise(MAF=max(MAF)) %>%
    arrange(desc(Het)) %>%
    group_by(Population) %>%
    top_n(Het, n=ceiling(2*nrow(Variants.hom)/n.pop)) %>%
    arrange(Population, desc(Het))

  ## Pick the homogeneous variants, making sure they belong to distinct groups
  n.hom <- n.causal - n.het
  Causal.hom <- Variants.hom %>%
    head(1.5*n.hom) %>%
    mutate(Population = NA) %>%
    ungroup() %>%
    distinct(Group, .keep_all=TRUE) %>%
    head(n.hom)

  ## Pick the population-specific variants, making sure they belong to distinct groups
  n.het <- round(n.causal*specificity)
  Causal.het <- Variants.het %>%
    group_by(Population) %>%
    top_n(Het, n=ceiling(2*n.het/n.pop)) %>%
    ungroup() %>%
    distinct(Group, .keep_all=TRUE) %>%
    group_by(Population) %>%
    top_n(Het, n=ceiling(1.5*n.het/n.pop)) %>%
    arrange(Population, desc(Het)) %>%
    head(1.5*n.het) %>%
    ungroup()

  ## Remove heterogeneous variants whose groups are already found in the homogeneous list
  Causal.het <- Causal.het %>%
    filter(!Group %in% Causal.hom$Group) %>%
    group_by(Population) %>%
    top_n(Het, n=ceiling(n.het/n.pop)) %>%
    arrange(Population, desc(Het)) %>%
    head(n.het) %>%
    ungroup()

  ## Combine the list of causal variants, making sure they all belong to distinct groups
  Causal <- rbind(Causal.hom, Causal.het)
  
  ## Make sure we have no duplicates
  stopifnot(length(unique(Causal$Group))==n.causal)

  return(Causal)
}

for(n.causal in c(100)) {
  for(specificity in c(1.0, 0.75, 0.5, 0.25, 0)) {
    Causal <- select_causal(n.causal, specificity)

    ## Store list of causal variants
    out.dir <- "/oak/stanford/groups/candes/matteo/transfer_knockoffs/data/causal"
    out.file <- sprintf("%s/causal_n%s_s%s.txt", out.dir, n.causal, round(100*specificity))
    Causal %>% write_delim(out.file, delim=" ")

    cat(sprintf("List of %d causal variants saved in: %s\n", n.causal, out.file))
  }
}
