
#####################W
### Making Figure - S3


############
### Fig-S3A: 

library(ggplot2)
library(tidyr)
library(dplyr)

data_path = "/Users/cxqiu/GitHub/mouse_sprint/figures/fig_S3"
save_path = "/Volumes/f0085ts/work/tapemouse/making_figures"

dat = read.csv(paste0(data_path, "/panel_B_ltt_vs_ceiling.csv"))

dat$log10_count = log10(dat$count)

p = ggplot(dat[dat$log10_count > 1,], aes(x = day, y = log10_count, color = series)) +
    geom_point() +
    geom_line() +
    scale_x_continuous(breaks = seq(1.5, 13.5, 2)) +
    scale_y_continuous(breaks = 1:7) +
    labs(x = "Time (embryonic day)", y = "Number of inferred lineages or cells in embryo", fill = NULL) +
    theme_classic() +
    scale_color_manual(values = c("unconstrained dating" = "#4076bb", "cell ceiling, used (≤ E7.5)" = "#d14b5c",
                                  "cell ceiling, unused (> E7.5)" = "#8a8881", "constrained dating" = "#108441")) +
    theme(legend.position = "none")

ggsave(paste0(save_path, "/FigS3/FigS3A.pdf"), p, height = 5, width = 6.5)




############
### Fig-S3B: 

library(ggplot2)
library(tidyr)
library(dplyr)
library(ggrastr)

data_path = "/Users/cxqiu/GitHub/mouse_sprint/figures/fig_S3"
save_path = "/Volumes/f0085ts/work/tapemouse/making_figures"

dat = read.csv(paste0(data_path, "/panel_C_date_shift.csv"))

p = ggplot(dat, aes(x = old, y = shift, color = blast)) +
    rasterise(geom_point(size = 0.1), dpi = 300) + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_x_continuous(breaks = seq(1.5, 13.5, 2)) +
    scale_y_continuous(breaks = seq(0, 4, 1)) +
    labs(x = "Unconstrained node date (embryonic day)", 
         y = "Date shift, new - old (days)", fill = NULL) +
    theme_classic() +
    scale_color_manual(values = c("B1" = "#4783B5", "B2" = "#F78C1E")) +
    theme(legend.position = "none")

ggsave(paste0(save_path, "/FigS3/FigS3B.pdf"), p, height = 5, width = 6.5)


############
### Fig-S3C: 

library(ggplot2)
library(tidyr)
library(dplyr)

data_path = "/Users/cxqiu/GitHub/mouse_sprint/figures/fig_S3"
save_path = "/Volumes/f0085ts/work/tapemouse/making_figures"

dat = read.csv(paste0(data_path, "/panel_A_support_histogram.csv"))

dat_frac <- dat %>%
    group_by(blastomere, support) %>%
    summarise(n = n(), .groups = "drop_last") %>%
    mutate(frac = n / sum(n)) %>%
    ungroup()

medians <- dat %>%
    group_by(blastomere) %>%
    summarise(med = median(support), .groups = "drop")

p = ggplot(dat_frac, aes(x = support, y = frac, fill = blastomere)) +
    geom_col(position = position_dodge(preserve = "single")) +
    geom_vline(data = medians,
               aes(xintercept = med, color = blastomere),
               linetype = "dashed", linewidth = 0.6, show.legend = FALSE) +
    scale_x_continuous(breaks = seq(0, max(dat_frac$support), by = 5)) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Support", y = "Fraction of nodes", fill = NULL) +
    theme_classic() +
    scale_fill_manual(values = c("B1" = "#4783B5", "B2" = "#F78C1E")) +
    theme(legend.position = "none")

ggsave(paste0(save_path, "/FigS3/FigS3C.pdf"), p, height = 5, width = 6.5)
