
#####################
### Making Figure - 3


###########
### Fig-3D: 

library(ggplot2)
library(tidyr)
library(dplyr)

data_path = "/Users/cxqiu/GitHub/mouse_sprint/figures/fig_3"
save_path = "/Volumes/f0085ts/work/tapemouse/making_figures"

dat = read.csv(paste0(data_path, "/panel_F_temporal_resolution.csv"))

dat$rate[dat$architecture == "empirical_11x6" & dat$day == 0] = 13.33

p = ggplot(dat, aes(x = day, y = rate, color = architecture)) +
    geom_line(size = 2) +
    scale_x_continuous(breaks = c(seq(1.5, 13.5, 2))) +
    scale_y_continuous(breaks = seq(0, 16, 2), limits = c(0, 16)) +
    labs(x = "Time (embryonic day)", y = "Edits per cell per day", fill = NULL) +
    theme_classic() +
    scale_color_manual(values = c("sequential_1x66_flat" = "#4076bb", "non_sequential_1x66" = "#d14b5c",
                                  "empirical_11x6" = "black", "sequential_11x6_constant_rate" = "#108441")) +
    theme(legend.position = "none")

ggsave(paste0(save_path, "/Fig3/Fig3D.pdf"), p, height = 4, width = 5.5)





###########
### Fig-3E: 

library(ggplot2)
library(tidyr)
library(dplyr)

data_path = "/Users/cxqiu/GitHub/mouse_sprint/figures/fig_3"
save_path = "/Volumes/f0085ts/work/tapemouse/making_figures"

dat = read.csv(paste0(data_path, "/panel_D_editing_rate.csv"))

p <- ggplot(dat[!is.na(dat$rate),], aes(x = day_mid, y = rate, color = blastomere)) +
    geom_point() +
    geom_line() +
    scale_x_continuous(breaks = c(seq(1.5, 13.5, 2))) +
    scale_y_continuous(breaks = seq(0, 16, 2), limits = c(0, 16)) +
    scale_color_manual(values = c("A" = "#4783B5", "B" = "#F78C1E")) +
    theme_classic() +
    theme(legend.position = "none") +
    labs(x = "Time (days)", y = "Editing rate (edits/day)", fill = NULL)

ggsave(paste0(save_path, "/Fig3/Fig3E.pdf"), p, height = 4, width = 5)




###########
### Fig-3F: 

library(ggplot2)
library(tidyr)
library(dplyr)

data_path = "/Users/cxqiu/GitHub/mouse_sprint/figures/fig_3"
save_path = "/Volumes/f0085ts/work/tapemouse/making_figures"

dat = read.csv(paste0(data_path, "/panel_E_support_over_time.csv"))

# Reshape to long format so both pct columns can be plotted as separate lines
dat_long <- dat %>%
    pivot_longer(cols = c(pct_ge1, pct_ge2),
                 names_to = "threshold", values_to = "pct") %>%
    mutate(threshold = recode(threshold,
                              pct_ge1 = ">1",
                              pct_ge2 = ">2"))

p <- ggplot(dat_long[!is.na(dat_long$pct),], aes(x = day_mid, y = pct,
                          color = blastomere, linetype = threshold,
                          group = interaction(blastomere, threshold))) +
    geom_line() +
    geom_point() +
    scale_x_continuous(breaks = c(seq(1.5, 13.5, 2))) +
    scale_color_manual(values = c("A" = "#4783B5", "B" = "#F78C1E")) +
    scale_linetype_manual(values = c(">1" = "solid", ">2" = "dotted")) +
    theme_classic() +
    theme(legend.position = "none") +
    labs(x = "Time (days)", y = "% of internal nodes",
         color = "Blastomere", linetype = "Threshold")

ggsave(paste0(save_path, "/Fig3/Fig3F.pdf"), p, height = 4, width = 5)

