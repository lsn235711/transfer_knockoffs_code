library(tidyverse)
library(gridExtra)

oak <- "/oak/stanford/groups/candes/matteo/transfer_knockoffs"

pop.list <- c("African", "African-small", "Asian", "Asian-small", "Indian", "Indian-small", "European", "European-small", "British", "Everyone")
ncausal.list <- c(100)
snr.list <- c(5, 10)
spec.list <- c(0, 50, 100)
fold.list <- c(1:10)
transfer.list <- c("_transfer", "")

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
  res.fdp <- df %>%
    filter(Discovered) %>%  
    summarise(FDP = mean(!Causal)) %>%
    as.numeric()
  res.pow <- df %>%
    filter(Causal) %>%  
    summarise(Power = mean(Discovered)) %>%
    as.numeric()
  res <- Config[row,] %>% mutate(FDP=res.fdp, Power=res.pow) %>%
    mutate(Transfer = ifelse(Transfer=="", FALSE, TRUE))
  return(res)
}))
Results <- Results %>%
  mutate(Method = ifelse(Transfer, "Transfer knockoffs - Lasso", "Vanilla knockoffs")) %>%
  select(-Transfer)

Config <- expand.grid(Population=pop.list, N.causal=ncausal.list, SNR=snr.list, Specificity=spec.list, Fold=fold.list) %>%
  as_tibble()

Results.linear <- do.call("rbind", lapply(1:nrow(Config), function(row) {
  n.causal <- Config$N.causal[row]
  snr <- Config$SNR[row]
  specificity <- Config$Specificity[row]
  pop <- Config$Population[row]
  fold <- Config$Fold[row]
  ifile <- sprintf("%s/discoveries_others/discoveries_transfer_2_res5_n%s_snr%s_s%s_%s_fold%d.txt",
                   oak, n.causal, round(snr), round(specificity), pop, fold)
  if(!file.exists(ifile)) {
    return(tibble())
  }
  df <- read_delim(ifile, delim=" ", col_types=cols())
  res.fdp <- df %>%
    group_by(Theta) %>%
    summarise(FDP = mean(!Causal[Discovered])) %>%
    ungroup()
  res.pow <- df %>%
    group_by(Theta) %>%
    summarise(Power = mean(Discovered[Causal]))
  res <- Config[row,] %>% cbind(full_join(res.fdp, res.pow)) %>% as_tibble() %>%
    mutate(Transfer = Theta) %>%
    select(-Theta)
  return(res)
}))

## Combine results for method with hyper-parameters
Results.linear <- Results.linear %>%
  mutate(Method = "Transfer knockoffs - Linear", Theta = Transfer) %>%
  group_by(Population, N.causal, SNR, Specificity, Fold, Method) %>%
  summarise(idx = which.max(Power), Power=Power[idx], FDP=FDP[idx]) %>%
  select(-idx) %>%
  ungroup()

Config <- expand.grid(Population=pop.list, N.causal=ncausal.list, SNR=snr.list, Specificity=spec.list, Fold=fold.list) %>%
  as_tibble()

Results.adaptive <- do.call("rbind", lapply(1:nrow(Config), function(row) {
  n.causal <- Config$N.causal[row]
  snr <- Config$SNR[row]
  specificity <- Config$Specificity[row]
  pop <- Config$Population[row]
  fold <- Config$Fold[row]
  ifile <- sprintf("%s/discoveries_others/discoveries_transfer_3_res5_n%s_snr%s_s%s_%s_fold%d.txt",
                   oak, n.causal, round(snr), round(specificity), pop, fold)
  if(!file.exists(ifile)) {
    return(tibble())
  }
  df <- read_delim(ifile, delim=" ", col_types=cols())
  res.fdp <- df %>%
    summarise(FDP = mean(!Causal[Discovered]))
  res.pow <- df %>%
    summarise(Power = mean(Discovered[Causal]))
  res <- Config[row,] %>% mutate(FDP=res.fdp$FDP, Power=res.pow$Power) %>% as_tibble()
  return(res)
}))

Results.adaptive <- Results.adaptive %>%
  mutate(Method = "Transfer knockoffs - Adaptive")

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
  group_by(Population, N.causal, SNR, Specificity, Method, Key) %>%
  summarise(Value.mean=mean(Value), Value.se=2*sd(Value)/sqrt(n())) %>%
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
out.file <- sprintf("%s/results/summary.txt", oak)
Summary %>% write_delim(out.file, delim=" ")
