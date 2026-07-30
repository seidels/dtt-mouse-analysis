

#################################
### Making Supplementary Figure 5

##########################################################
### Supplementary Fig-5a:  active_typeheads_by_integration

library(ggplot2)
library(tidyr)
library(dplyr)

data_path = "/Users/cxqiu/GitHub/mouse_sprint/tape_pipeline/tables"
save_path = "/Volumes/f0085ts/work/tapemouse/making_figures"

dat = read.csv(paste0(data_path, "/e3.active_typeheads_by_integration.csv"), skip=1)

dat_long = pivot_longer(dat,
                        cols = -c(time_mid_days, n_nodes),
                        names_to = "integration",
                        values_to = "pct_active")
dat_long$time_mid_days = factor(dat_long$time_mid_days, levels = unique(dat$time_mid_days))

p = ggplot(dat_long, aes(x = time_mid_days, y = pct_active,
                     color = integration, group = integration)) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.8) +
    scale_color_manual(values = tapebc_color_plate) +
    theme_classic(base_size = 12) +
    theme(axis.text = element_text(color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Time (days)",
         y = "% of nodes with active typehead",
         color = "Integration") +
    ylim(0, 100)

ggsave(paste0(save_path, "/FigS4_active_typeheads_by_integration.pdf"), p, height =4, width = 8)



######################################################
### Supplementary Fig-5b:  e3.active_typeheads_by_time

library(ggplot2)
library(tidyr)
library(dplyr)
library(viridis)

data_path = "/Users/cxqiu/GitHub/mouse_sprint/tape_pipeline/tables"
save_path = "/Volumes/f0085ts/work/tapemouse/making_figures"

dat = read.csv(paste0(data_path, "/e3.active_typeheads_by_time.csv"), skip=2)
dat$time_mid_days = factor(dat$time_mid_days, levels = unique(dat$time_mid_days))

p = ggplot(dat, aes(x = time_mid_days, y = pct_active_mean,
                         color = time_mid_days, group = 1)) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.8) +
    scale_color_viridis(discrete=TRUE) +
    theme_classic(base_size = 12) +
    theme(axis.text = element_text(color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Time (days)",
         y = "% of nodes with active typehead",
         color = "Integration") +
    ylim(0, 100)

ggsave(paste0(save_path, "/FigS4_active_typeheads_by_time.pdf"), p, height =4, width = 8)



