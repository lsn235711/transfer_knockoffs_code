library(tidyverse)
library(gridExtra)

N = 10
FDPs = NULL
POWERs = NULL
overlaps = NULL
Methods = NULL 
for (arg in 1:300){
    filename = sprintf("/Users/lsn/Desktop/Can/project_transfer_learning/simulations/simulation_20210709_transfer_gam/sherlock/output_trans/output%i.Rdata", arg)
    if(file.exists(filename)){
        load(filename)
        overlap = 20*((arg-1)%%6)
        FDPs = c(FDPs, fdps)
        POWERs = c(POWERs, powers)
        overlaps = c(overlaps, rep(overlap,N*10))
        Methods = c(Methods, methods)
    }
}
results = data.frame(overlap = overlaps, power = POWERs, fdp = FDPs, method = Methods)

df.nominal <- tibble(key="FDR", value=0.1)

method.values <- c("pool", "vanilla", "adaptive", paste("linear_com", 0.1*(1:5), sep = "_"), "oracle", "gam", "weighted_lasso")
method.labels <- c("Heuristic (pool)", "Vanilla", "Adaptive (EM)", rep("Transfer - linearly re-ordered",5),"Transfer - linearly re-ordered (oracle)" , "Transfer - adaptive (gam)", "Transfer - weighted-lasso")
#method.labels <- c("Heuristic (pool)", "Vanilla knockoffs", "Adaptive knockoffs", paste("Linear combination statistics, theta = ", 0.1*(1:5), sep = ""),"Weighted Lasso Statistics")

## find the oracle method
# NN = 300
# power_matrix = cbind(POWERs[10*(0:(NN-1)) + 5],POWERs[10*(0:(NN-1)) + 6],POWERs[10*(0:(NN-1)) + 7], POWERs[10*(0:(NN-1)) + 8], POWERs[10*(0:(NN-1)) + 9])
# which_theta = apply(power_matrix, 1, which.max)
# oracle_index = 10*(0:(NN-1)) + which_theta + 4
# results_oracle = data.frame(overlap = overlaps[oracle_index], power = POWERs[oracle_index], fdp = FDPs[oracle_index], method = rep("oracle", NN))
# results = rbind(results, results_oracle)


linear_summary = results %>%
    filter(method == "linear_com_0.1"|
               method == "linear_com_0.2"|
               method == "linear_com_0.3"|
               method == "linear_com_0.4"|
               method == "linear_com_0.5") %>%
    mutate(Power = power, FDR = fdp) %>%
    group_by(overlap, method) %>%
    summarise(FDR = mean(FDR), Power = mean(Power)) %>%
    group_by(overlap) %>%
    summarise(idx = which.max(Power), Power=Power[idx], FDR = FDR[idx]) %>%
    select(-idx) %>%
    gather(Power, FDR, key="key", value="value") %>%
    mutate(method = "oracle")

others_summary = results %>%
    filter(method != "linear_com_0.1") %>%
    filter(method != "linear_com_0.2") %>%
    filter(method != "linear_com_0.3") %>%
    filter(method != "linear_com_0.4") %>%
    filter(method != "linear_com_0.5") %>%
    filter(method != "adaptive") %>%
    mutate(Power = power, FDR = fdp) %>%
    gather(Power, FDR, key="key", value="value") %>%
    group_by(overlap, method, key) %>%
    summarise(num=n(), value.se=2*sd(value)/sqrt(num), value=mean(value, na.rm = TRUE))

full_summary = full_join(linear_summary, others_summary)
p1 = full_summary %>%
    mutate(Method = factor(method, method.values, method.labels)) %>%
    #mutate(value = ifelse(is.na(value), 0, value)) %>%
    ggplot(aes(x=overlap, y=value, color=Method, shape=Method)) +
    #ggplot(aes(x=overlap, y=value, color=Method)) +
    geom_point(alpha=0.8) +
    geom_line(alpha=0.8) +
    #geom_errorbar(aes(ymin=value-value.se, ymax=value+value.se), width=0.2) +
    geom_hline(data=df.nominal, aes(yintercept=value), linetype=2) +
    scale_color_brewer(palette="Set1") +
    #scale_shape_manual(values=shape.scale) +
    facet_wrap(.~key, nrow=1, scales="free_y") +
    xlab("Overlap (in %)") +
    ylab(NULL) +
    ylim(0,1) +
    theme_bw() +
    theme(#strip.background = element_blank(),
          #strip.placement = "outside",
          legend.position = "right",
          plot.title = element_text(size = 10), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
          strip.text.x = element_text(size = 10), strip.text.y = element_text(size = 10), legend.title=element_text(size = 10)) +
    guides(col = guide_legend(nrow = 5)) #+
    #ggtitle("(a) Setting 1 (high invariance)")

p1

p1 %>% ggsave(file="experiment_transfer.pdf", width=7.8, height=2.7, units="in")


#### Compare different theta
#library(viridis)
library(RColorBrewer)
library(wesanderson)
method.values <- c("pool", "vanilla", "adaptive", paste("linear_com", 0.1*(1:5), sep = "_"), "weighted_lasso")
#method.labels <- c("Heuristic (pool)", "Vanilla knockoffs", "Adaptive knockoffs", rep("Linear combination statistics",5),"Weighted Lasso Statistics")
method.labels <- c("Heuristic (pool)", "Vanilla", "Transfer - adaptive (gam)", paste("Transfer - linearly reordered\n theta = ", 0.1*(1:5), sep = ""),"Transfer - weighted lasso")

#color_range = colorRampPalette(c("darkorchid4", "gold3")) 
color.scale <- c(brewer.pal(n=5,"Set1")[1:2], wes_palette("Moonrise1")[2],wes_palette("Moonrise2")[2], wes_palette("Cavalcanti1")[c(3,2)])


p2 <- results %>%
    filter(method != "adaptive") %>%
    filter(method != "gam") %>%
    filter(method != "weighted_lasso") %>%
    filter(method != "linear_com_0.5") %>%
    filter(method != "oracle") %>%
    mutate(Method = factor(method, method.values, method.labels)) %>%
    mutate(Power = power, FDR = fdp) %>%
    gather(Power, FDR, key="key", value="value") %>%
    group_by(overlap, Method, key) %>%
    summarise(num=n(), value.se=2*sd(value)/sqrt(num), value=mean(value, na.rm = TRUE)) %>%
    #mutate(value = ifelse(is.na(value), 0, value)) %>%
    ggplot(aes(x=overlap, y=value, color=Method, shape=Method)) +
    #ggplot(aes(x=overlap, y=value, color=Method)) +
    geom_point(alpha=0.7) +
    geom_line(alpha=0.7) +
    #geom_errorbar(aes(ymin=value-value.se, ymax=value+value.se), width=0.2) +
    geom_hline(data=df.nominal, aes(yintercept=value), linetype=2) +
    scale_color_manual(values=color.scale) +
    facet_wrap(.~key, nrow=1, scales="free_y") +
    xlab("Overlap (in %)") +
    ylab(NULL) +
    ylim(0,1) +
    theme_bw() +
    theme(#strip.background = element_blank(),
        #strip.placement = "outside",
        legend.position = "right",
        plot.title = element_text(size = 10), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        strip.text.x = element_text(size = 10), strip.text.y = element_text(size = 10), legend.title=element_text(size = 10)) +
    guides(col = guide_legend(nrow = 6)) #+
#ggtitle("(a) Setting 1 (high invariance)")

p2
p2 %>% ggsave(file="experiment_transfer_linear.pdf", width=6.9, height=2.7, units="in")

