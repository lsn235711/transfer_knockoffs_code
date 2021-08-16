library(tidyverse)
library(gridExtra)

snr <- 5

ifile <- "../results/summary.txt"
Summary <- read_delim(ifile, delim=" ")
methods.values <- c("Vanilla knockoffs", "Transfer knockoffs - Linear", "Transfer knockoffs - Adaptive", "Transfer knockoffs - Lasso")
methods.labels <- c("Vanilla", "Transfer - linearly re-ordered (oracle)", "Transfer - adaptive (gam)", "Transfer - weighted-lasso")
color.scale <- c("#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")
shape.scale <- c(17,15,3,7)
linetype.scale <- c(1,1,1,1)
Summary <- Summary %>%
    mutate(Method = factor(Method, methods.values, methods.labels))
df.dummy <- tibble(Key="FDP", Value=0.1)
## Plot with equal sample sizes
p1 <- Summary %>%
    filter(Population!="Everyone", SNR==snr) %>%
    mutate(Full = ifelse(endsWith(Population, "-small"), FALSE, TRUE), Population = str_replace(Population, "-small", "")) %>%
    filter(!Full) %>%
    mutate(Population = sprintf("%s (n = %d)", Population, Samples)) %>%
    ggplot(aes(x=Specificity, y=Value.mean, color=Method, linetype=Method, shape=Method)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=pmax(0,Value.mean-Value.se), ymax=Value.mean+Value.se), width=5) +
    geom_hline(data=df.dummy, aes(yintercept=Value), linetype=2) +
    facet_grid(Key~Population) +
    scale_x_continuous(breaks=c(0,50,100)) +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_linetype_manual(values=linetype.scale) +
    xlab("Heterogeneity of causal variants (%)") +
    ylab("") +
    ylim(0,0.4) +
    guides(color = guide_legend(nrow = 2)) +
    theme_bw() +
    theme(legend.position="bottom", legend.key.size = grid::unit(2, "lines"))

ggsave(sprintf("../figures/transfer_snr%s_small.png", snr), p1, height=4, width=6, units="in")

## Plot with different sample sizes
p2 <- Summary %>%
    filter(SNR==snr) %>%
    filter(Population!="British", ! ((Population == "British")*(Method!="Vanilla"))) %>%
    mutate(Full = ifelse(endsWith(Population, "-small"), FALSE, TRUE), Population = str_replace(Population, "-small", "")) %>%
    mutate(Population = ifelse(Population=="Everyone", "Pooled", Population)) %>%
    filter(Full) %>%
    mutate(Population = sprintf("%s (n = %d)", Population, Samples)) %>%
    ggplot(aes(x=Specificity, y=Value.mean, color=Method, linetype=Method, shape=Method)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=pmax(0,Value.mean-Value.se), ymax=Value.mean+Value.se), width=5) +
    geom_hline(data=df.dummy, aes(yintercept=Value), linetype=2) +
    facet_grid(Key~Population) +
    scale_x_continuous(breaks=c(0,50,100)) +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_linetype_manual(values=linetype.scale) +
    xlab("Heterogeneity of causal variants (%)") +
    ylab("") +
    ylim(0,1) +
    guides(color = guide_legend(nrow = 2)) +
    theme_bw() +
    theme(legend.position="bottom", legend.box="vertical", legend.key.size = grid::unit(2, "lines"))

ggsave(sprintf("../figures/transfer_snr%s.png", snr), p2, height=4, width=8, units="in")


#########################
## Plots by population ##
#########################

ifile <- "../results/summary_separate.txt"
Summary <- read_delim(ifile, delim=" ")

methods.values <- c("Pooling", "Vanilla knockoffs", "Transfer knockoffs - Linear", "Transfer knockoffs - Adaptive", "Transfer knockoffs - Lasso")
methods.labels <- c("Vanilla on British population", "Transfer - Linear combination (oracle)", "Vanilla", "Transfer - Adaptive", "Transfer - Weighted lasso")
color.scale <- c("#377EB8", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")
shape.scale <- c(1,17,15,3,7)
linetype.scale <- c(2,1,1,1,1)

Summary <- Summary %>%
    mutate(Method = factor(Method, methods.values, methods.labels))

df.dummy <- tibble(Key="FDP", Value=0.1)

df.dummy <- tibble(Key="FDR (population-specific)", Value=0.1)
methods.values <- c("Vanilla on British population", "Vanilla", "Transfer - Linear combination (oracle)", "Transfer - Adaptive", "Transfer - Weighted lasso")
methods.labels <- c("Vanilla on British", "Vanilla",  "Transfer - linearly-reordered combination (oracle)", "Transfer - adaptive (gam)", "Transfer - weighted-lasso")

p1 <- Summary %>%
    filter(SNR==snr) %>%
    filter((Population==Pop)|(Population=="British")) %>%
    mutate(Method = ifelse(Population=="British", "Vanilla on British population", as.character(Method))) %>%
    mutate(Method = factor(Method, methods.values, methods.labels)) %>%
    mutate(Key = ifelse(Key=="FDP", "FDR (population-specific)", Key),
           Key = ifelse(Key=="Power", "Power (population-specific)", Key)) %>%
    ggplot(aes(x=Specificity, y=Value.mean, color=Method, linetype=Method, shape=Method)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=pmax(0,Value.mean-Value.se), ymax=Value.mean+Value.se), width=5) +
    geom_hline(data=df.dummy, aes(yintercept=Value), linetype=2) +
    facet_grid(Key~Pop) +
    scale_x_continuous(breaks=c(0,50,100)) +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_linetype_manual(values=linetype.scale) +
    xlab("Heterogeneity of causal variants (%)") +
    ylab("") +
    ylim(0,1) +
    guides(color = guide_legend(nrow = 2)) +
    theme_bw() +
    theme(legend.position="bottom", legend.key.size = grid::unit(2, "lines"))

ggsave(sprintf("../figures/transfer_specific_snr%s.png", snr), p1, height=5, width=8, units="in")

methods.values <- c("Pooling", "Vanilla", "Transfer - Linear combination (oracle)", "Transfer - Adaptive", "Transfer - Weighted lasso")
methods.labels <- c("Heuristic (pool)", "Vanilla", "Transfer - linearly-reordered combination (oracle)", "Transfer - adaptive (gam)", "Transfer - weighted-lasso")
color.scale <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")
shape.scale <- c(16,17,15,3,7)
linetype.scale <- c(1,1,1,1,1)


p2 <- Summary %>%
    filter(SNR==snr) %>%
    filter((Population==Pop)|(Population=="Everyone")) %>%
    mutate(Method = ifelse(Population=="Everyone", "Pooling", as.character(Method))) %>%
    mutate(Method = factor(Method, methods.values, methods.labels)) %>%
    mutate(Key = ifelse(Key=="FDP", "FDR (population-specific)", Key),
           Key = ifelse(Key=="Power", "Power (population-specific)", Key)) %>%
    ggplot(aes(x=Specificity, y=Value.mean, color=Method, linetype=Method, shape=Method)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=pmax(0,Value.mean-Value.se), ymax=Value.mean+Value.se), width=5) +
    geom_hline(data=df.dummy, aes(yintercept=Value), linetype=2) +
    facet_grid(Key~Pop) +
    scale_x_continuous(breaks=c(0,50,100)) +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_linetype_manual(values=linetype.scale) +
    xlab("Heterogeneity of causal variants (%)") +
    ylab("") +
    ylim(0,1) +
    guides(color = guide_legend(nrow = 2)) +
    theme_bw() +
    theme(legend.position="bottom", legend.key.size = grid::unit(2, "lines"))

ggsave(sprintf("../figures/transfer_specific_pooled_snr%s.png", snr), p2, height=5, width=8, units="in")
