
#################################
### Making Supplementary Figure 1

#################################################################################
### Supplementary Fig-1c: TAPE-BC read-count rank-abundance, one facet per embryo

library(ggplot2)
library(dplyr)

data_path = "/Users/cxqiu/GitHub/mouse_sprint/tape_pipeline/tables"
save_path = "/Volumes/f0085ts/work/tapemouse/making_figures"

dat = read.csv(paste0(data_path, "/suppfig1a_barcode_rank_abundance.csv"))
dat = dat[dat$embryo %in% paste0("embryo ", c(2,3,6)),]
dat$log_reads = log(dat$reads)
dat$copies = factor(dat$copies, levels = names(table(dat$copies)))

dat_cutoff = read.csv(paste0(data_path, "/suppfig1a_cutoffs.csv"))
dat_cutoff = dat_cutoff[dat_cutoff$embryo %in% paste0("embryo ", c(2,3,6)),]
dat_cutoff$log_reads_at_freq_cutoff = log(dat_cutoff$reads_at_freq_cutoff)

p = dat %>%
    ggplot(aes(rank, log_reads, color = copies)) +
    geom_point(size = 0.6) +
    geom_hline(data = dat_cutoff, aes(yintercept = log_reads_at_freq_cutoff), linetype = "dashed") +
    facet_wrap(~ embryo, nrow = 1) +
    labs(x = "Barcode rank (ordered high to low)", y = "Log (Reads per barcode)", fill = NULL) +
    scale_y_continuous(breaks = seq(0, 20, by = 2)) +
    theme_classic()

ggsave(paste0(save_path, "/FigS1/FigS1_TAPE-BC_read-count_rank-abundance.pdf"), p, height = 3, width = 10)


##############################################
### Supplementary Fig-1e: rarefaction analysis

library(ggplot2)
library(dplyr)

data_path = "/Users/cxqiu/GitHub/mouse_sprint/tape_pipeline/tables"
save_path = "/Volumes/f0085ts/work/tapemouse/making_figures"

dat = read.csv(paste0(data_path, "/DTTz_3_S3.rarefaction_curve.csv"))

p = ggplot(data = filter(dat, region == "observed"), aes(x = cells_sampled, y = expected_lineages, color = tapebc)) +
    geom_point() +
    scale_color_manual(values = tapebc_color_plate) +
    labs(x = "# of genome equivalents sampled", y = "Expected distinct lineage genotypes", fill = NULL) +
    theme_classic()

ggsave(paste0(save_path, "/FigS1/FigS1_rarefaction_curve.pdf"), p, height = 3, width = 6)













