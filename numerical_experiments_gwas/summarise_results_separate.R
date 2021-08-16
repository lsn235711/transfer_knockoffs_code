library(tidyverse)
library(gridExtra)

oak <- "/oak/stanford/groups/candes/matteo/transfer_knockoffs"

pop.list <- c("African", "African-small", "Asian", "Asian-small", "Indian", "Indian-small", "European", "European-small", "British", "Everyone")
ncausal.list <- c(100)
snr.list <- c(5, 10)
spec.list <- c(0, 50, 100)
fold.list <- c(1:10)
transfer.list <- c("_transfer", "")

## Load MAF info
in.file <- sprintf("%s/data/qc/maf_info.txt", oak)
MAF <- read_delim(in.file, delim=" ") %>%
  select(-CHR, -Group, -BP, -Hom, -Het, -MAF.max)

Config <- expand.grid(Population=pop.list, N.causal=ncausal.list, SNR=snr.list, Specificity=spec.list, Fold=fold.list,
                      Transfer = transfer.list) %>%
  as_tibble()

Results <- do.call("rbind", lapply(1:nrow(Config), function(row) {
  n.causal <- Config$N.causal[row]
  snr <- Config$SNR[row]
  specificity <- Config$Specificity[row]  
  pop <- Config$Population[row]
  fold <- Config$Fold[row]
  transfer = Config$Transfer[row]
  ifile <- sprintf("%s/discoveries/discoveries%s_res5_n%s_snr%s_s%s_%s_fold%d.txt",
                   oak, transfer, n.causal, round(snr), round(specificity), pop, fold)
  if(!file.exists(ifile)) {
    return(tibble())
  }
  df <- read_delim(ifile, delim=" ", col_types=cols())
  ## Load list of causal variants
  in.file <- sprintf("%s/data/causal/causal_n%s_s%s.txt", oak, n.causal, round(specificity))
  Causal <- read_delim(in.file, delim=" ", col_types=cols()) %>%
    mutate(Causal = TRUE)
  df2 <- df %>%
    mutate(SNP = ifelse(is.na(SNP), SNP.lead, SNP)) %>%
    select(-SNP.lead) %>%
    full_join(Causal) %>%
#    mutate(Pop.c = Population) %>%
    select(Group, SNP, Causal, Discovered) %>%
    left_join(MAF) %>%
    mutate(Pop.m = Population.max) %>%
    select(-Population.max)
  res.fdp <- df2 %>%
    filter(Discovered) %>%
    group_by(Pop.m) %>%
    summarise(FDP = mean(!Causal)) %>%
    mutate(Pop = Pop.m) %>% select(-Pop.m)
  res.pow <- df2 %>%
    filter(Causal) %>%
    group_by(Pop.m) %>%
    summarise(Power = mean(Discovered)) %>%
    mutate(Pop = as.character(Pop.m)) %>% select(-Pop.m)    
  res <- Config[row,] %>% cbind(full_join(res.fdp, res.pow)) %>%
    mutate(Transfer = ifelse(Transfer=="", FALSE, TRUE))
  return(res)
}))
Results <- Results %>%
  mutate(Method = ifelse(Transfer, "Transfer knockoffs - Lasso", "Vanilla knockoffs")) %>%
  select(-Transfer) %>%
  as_tibble()

Config <- expand.grid(Population=pop.list, N.causal=ncausal.list, SNR=snr.list, Specificity=spec.list, Fold=fold.list) %>%
  as_tibble()

Results.linear <- do.call("rbind", lapply(1:nrow(Config), function(row) {
  n.causal <- Config$N.causal[row]
  specificity <- Config$Specificity[row]
  pop <- Config$Population[row]
  snr <- Config$SNR[row]
  fold <- Config$Fold[row]
  ifile <- sprintf("%s/discoveries_others/discoveries_transfer_2_res5_n%s_snr%s_s%s_%s_fold%d.txt",
                   oak, n.causal, round(snr), round(specificity), pop, fold)
  if(!file.exists(ifile)) {
    return(tibble())
  }
  df <- read_delim(ifile, delim=" ", col_types=cols())
  ## Load list of causal variants
  ##in.file <- sprintf("%s/data/causal/causal_n%s_s%s.txt", oak, n.causal, round(specificity))
  ##Causal <- read_delim(in.file, delim=" ", col_types=cols()) %>%
  ## mutate(Causal = TRUE)
  df2 <- df %>%
    mutate(SNP=as.character(SNP)) %>%
    mutate(SNP = ifelse(is.na(SNP), SNP.c, SNP)) %>%
    select(Theta, Group, SNP, Causal, Discovered) %>%
    left_join(MAF, by="SNP") %>%
    mutate(Pop.m = Population.max) %>%
    select(-Population.max)
  res.fdp <- df2 %>%
    filter(Discovered) %>%
    group_by(Theta, Pop.m) %>%
    summarise(FDP = mean(!Causal)) %>%
    mutate(Pop = Pop.m) %>% select(-Pop.m) %>%
    ungroup()
  res.pow <- df2 %>%
    group_by(Theta, Pop.m) %>%
    summarise(True=sum(Causal[Discovered]), Power = mean(Discovered[Causal])) %>%
    mutate(Pop = as.character(Pop.m)) %>% select(-Pop.m) %>%
    ungroup()
  res <- Config[row,] %>% cbind(full_join(res.fdp, res.pow, by=c("Theta", "Pop"))) %>%
    mutate(Transfer = Theta) %>%
    select(-Theta)
  return(res)
}))
Results.linear <- Results.linear %>% as_tibble()

Results.linear <- Results.linear %>%
  mutate(Method = "Transfer knockoffs - Linear", Theta = Transfer) %>%
  group_by(Population, N.causal, SNR, Specificity, Fold, Method, Theta) %>%
  mutate(True.total = sum(True)) %>%
  group_by(Population, N.causal, SNR, Specificity, Fold, Method, Pop) %>%
  summarise(idx = which.max(Power), Power=Power[idx], FDP=FDP[idx]) %>%
  select(-idx) %>%
  ungroup()
  
if(TRUE) {
  Config <- expand.grid(Population=pop.list, N.causal=ncausal.list, SNR=snr.list, Specificity=spec.list, Fold=fold.list) %>%
    as_tibble()
  Results.adaptive <- do.call("rbind", lapply(1:nrow(Config), function(row) {
    n.causal <- Config$N.causal[row]
    specificity <- Config$Specificity[row]
    pop <- Config$Population[row]
    snr <- Config$SNR[row]
    fold <- Config$Fold[row]
    ifile <- sprintf("%s/discoveries_others/discoveries_transfer_3_res5_n%s_snr%s_s%s_%s_fold%d.txt",
                     oak, n.causal, round(snr), round(specificity), pop, fold)
    if(!file.exists(ifile)) {
      return(tibble())
    }
    df <- read_delim(ifile, delim=" ", col_types=cols())
    df2 <- df %>%
      mutate(SNP=as.character(SNP)) %>%
      mutate(SNP = ifelse(is.na(SNP), SNP.c, SNP)) %>%
      select(Group, SNP, Causal, Discovered) %>%
      left_join(MAF, by="SNP") %>%
      mutate(Pop.m = Population.max) %>%
      select(-Population.max)
    res.fdp <- df2 %>%
      filter(Discovered) %>%
      group_by(Pop.m) %>%
      summarise(FDP = mean(!Causal)) %>%
      mutate(Pop = Pop.m) %>% select(-Pop.m) %>%
      ungroup()
    res.pow <- df2 %>%
      group_by(Pop.m) %>%
      summarise(Power = mean(Discovered[Causal])) %>%
      mutate(Pop = as.character(Pop.m)) %>% select(-Pop.m) %>%
      ungroup()
    res <- Config[row,] %>% cbind(full_join(res.fdp, res.pow, by=c("Pop")))
    return(res)
  }))
  Results.adaptive <- Results.adaptive %>%
    mutate(Method = "Transfer knockoffs - Adaptive")
} else {
  Results.adaptive <- tibble()
}
Results.adaptive <- Results.adaptive %>% as_tibble()

Samples <- do.call("rbind", lapply(pop.list, function(pop) {
  ## Load list of samples
  if(length(str_split(pop, "-small")[[1]])==2) {
    n <- 3284
  } else {
    ifile <- sprintf("%s/data/qc/samples_%s.txt", oak, pop)
    n <- nrow(read_delim(ifile, delim=" ", col_types=cols(), col_names=FALSE))
  }
  return(tibble(Population=pop, Samples=n))
}))

Summary <- rbind(Results, Results.linear, Results.adaptive) %>%
  mutate(FDP = ifelse(is.nan(FDP), 0, FDP)) %>%
  gather(FDP, Power, key="Key", value="Value") %>%
  group_by(Population, N.causal, SNR, Specificity, Pop, Method, Key) %>%
  summarise(Value.mean=mean(Value), Value.se=2*sd(Value)/sqrt(n())) %>%
  mutate(Value.mean = ifelse(is.na(Value.mean), 0, Value.mean),
         Value.se = ifelse(is.na(Value.se), 0, Value.se)) %>%
  left_join(Samples, by=c("Population"))
  
if(FALSE) {
  p2 <- Summary %>%
    ggplot(aes(x=Specificity, y=Value.mean, color=Transfer, linetype=Transfer, shape=Transfer)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=Value.mean-Value.se, ymax=Value.mean+Value.se)) +
    facet_grid(Key~Population) +
    theme_bw()
  p2
}

## Save results
out.file <- sprintf("%s/results/summary_separate.txt", oak)
Summary %>% write_delim(out.file, delim=" ")
 
