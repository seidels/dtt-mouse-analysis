
#####################
### Making Figure - 2

#######################################
### Fig-2c: saturation_sites_written

library(ggplot2)
library(tidyr)
library(dplyr)

data_path = "/Users/cxqiu/GitHub/mouse_sprint/tape_pipeline/tables"
save_path = "/Volumes/f0085ts/work/tapemouse/making_figures"

dat = read.csv(paste0(data_path, "/e3v5v6.saturation_sites_written.csv"), skip = 2)

p = ggplot(dat, aes(x = sites_written, y = pct_of_cells)) +
    geom_bar(stat="identity") +
    labs(x = "# of sites written", y = "% of cells", fill = NULL) +
    theme_classic()

ggsave(paste0(save_path, "/Fig2/Fig2_saturation_sites_written.pdf"), p, height =4, width = 6)


#######################################
### Fig-2d: lineage_rarefaction_curve

library(ggplot2)
library(tidyr)
library(dplyr)

data_path = "/Users/cxqiu/GitHub/mouse_sprint/tape_pipeline/tables"
save_path = "/Volumes/f0085ts/work/tapemouse/making_figures"

dat = read.csv(paste0(data_path, "/e3v5v6.lineage_rarefaction_curve.csv"))

p = ggplot(data = filter(dat, region == "observed"), aes(x = cells_sampled, y = expected_genotypes)) +
#    geom_line(linewidth = 1) +
    geom_point() +
    labs(x = "# of genome equivalents sampled", y = "Expected distinct lineage genotypes", fill = NULL) +
    theme_classic()

ggsave(paste0(save_path, "/Fig2/Fig2_lineage_rarefaction_curve.pdf"), p, height = 4, width = 6)



#####################################################
### Fig-2e: cell-type-compositions between two clades

library(ggplot2)
library(tidyr)
library(dplyr)

data_path = "/Users/cxqiu/GitHub/mouse_sprint/tape_pipeline/tables"
save_path = "/Volumes/f0085ts/work/tapemouse/making_figures"
work_path = "/Volumes/f0085ts/work/tapemouse/tree_analysis"

dat = read.table(paste0(data_path, "/v6/e3v5v6.routing_labels.tsv.gz"), sep='\t', header=T)
celltype = read.table(paste0(work_path, "/cell_metadata.v6.txt"), sep='\t', header=T)
major_trajectory_celltype = read.table(paste0(work_path, "/major_trajectory_celltype_table.txt"), sep='\t', header=T)

all_celltypes = unique(celltype$celltype[celltype$cell_id %in% dat$cell[dat$blastomere %in% c("B1","B2")]])
### 135 cell types

cell_num_1 = celltype %>% 
    filter(cell_id %in% dat$cell[dat$blastomere == "B1"]) %>% 
    group_by(celltype) %>% 
    tally() %>% 
    complete(celltype = all_celltypes, fill = list(n = 0)) %>%
    mutate(total_n = sum(n)) %>%
    mutate(log2_frac = log2(100*(n/total_n)+1)) %>%
    select(celltype, B1_log2_frac = log2_frac)

cell_num_2 = celltype %>% 
    filter(cell_id %in% dat$cell[dat$blastomere == "B2"]) %>% 
    group_by(celltype) %>% 
    tally() %>% 
    complete(celltype = all_celltypes, fill = list(n = 0)) %>%
    mutate(total_n = sum(n)) %>%
    mutate(log2_frac = log2(100*(n/total_n)+1)) %>%
    select(celltype, B2_log2_frac = log2_frac)

df = cell_num_1 %>% left_join(cell_num_2, by = "celltype") %>%
    left_join(major_trajectory_celltype, by = "celltype")

p = ggplot(df, aes(x = B1_log2_frac, y = B2_log2_frac, color = major_trajectory)) +
    geom_point(size = 3) +
    theme_classic(base_size = 12) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    labs(x = "Log2[Fraction (%) + 1], B1 blastomere", y = "Log2[Fraction (%) + 1], B2 blastomere") +
    scale_color_manual(values=major_trajectory_color_plate)

ggsave(paste0(save_path, "/Fig2/Fig2_celltype_frac_two_blastomere.pdf"), p, height = 5, width = 5)

fit = cor.test(df$B1_log2_frac, df$B2_log2_frac, method = "spearman")
print(fit$estimate) ### 0.9952278
print(fit$p.value)  ### 2.664643e-136


############################################
### jensen shannon divergence on proportions

data_path = "/Users/cxqiu/GitHub/mouse_sprint/tape_pipeline/tables"
save_path = "/Volumes/f0085ts/work/tapemouse/making_figures"
work_path = "/Volumes/f0085ts/work/tapemouse/tree_analysis"

dat = read.table(paste0(data_path, "/v6/e3v5v6.routing_labels.tsv.gz"), sep='\t', header=T)
celltype = read.table(paste0(work_path, "/cell_metadata.v6.txt"), sep='\t', header=T)
major_trajectory_celltype = read.table(paste0(work_path, "/major_trajectory_celltype_table.txt"), sep='\t', header=T)

all_celltypes = unique(celltype$celltype[celltype$cell_id %in% dat$cell[dat$blastomere %in% c("B1","B2")]])

cell_num_1 = celltype %>% 
    filter(cell_id %in% dat$cell[dat$blastomere == "B1"]) %>% 
    group_by(celltype) %>% 
    tally() %>% 
    complete(celltype = all_celltypes, fill = list(n = 0)) %>%
    mutate(total_n = sum(n)) %>%
    mutate(frac = n/total_n) %>%
    select(celltype, B1_frac = frac)

cell_num_2 = celltype %>% 
    filter(cell_id %in% dat$cell[dat$blastomere == "B2"]) %>% 
    group_by(celltype) %>% 
    tally() %>% 
    complete(celltype = all_celltypes, fill = list(n = 0)) %>%
    mutate(total_n = sum(n)) %>%
    mutate(frac = n/total_n) %>%
    select(celltype, B2_frac = frac)

df = cell_num_1 %>% left_join(cell_num_2, by = "celltype") %>%
    left_join(major_trajectory_celltype, by = "celltype")

p <- df$B1_frac
q <- df$B2_frac
m <- (p + q) / 2

kl <- function(x, y) sum(x * log2(x / y), na.rm = TRUE)
jsd <- 0.5 * kl(p, m) + 0.5 * kl(q, m)
cat("JSD:", jsd, "\n")
### 0.002752065

library(philentropy)

props <- rbind(df$B1_frac, df$B2_frac)
jsd <- JSD(props, unit = "log2")
cat("Jensen-Shannon divergence (log2):", jsd, "\n")







