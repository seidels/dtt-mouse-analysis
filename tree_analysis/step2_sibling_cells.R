
#######################################################
### Do sibling cells share the same cell type identity?

source("~/work/scripts/utils.R")
library(dplyr)
library(tidyr)
library(ape)
library(ggplot2)
library(ggrepel)

work_path <- "/net/shendure/vol2/projects/cxqiu/work/tapemouse"

# ---- Load data ----
cell_meta <- read.table(paste0(work_path, "/tree_analysis/cell_metadata.v6.txt"),
                        header = TRUE, sep = "\t")

# Major trajectory
# cell_meta <- cell_meta %>% select(cell_id, celltype = major_trajectory) %>% as.data.frame()

# Cell type
cell_meta <- cell_meta %>% select(cell_id, celltype) %>% as.data.frame()

tree <- read.tree(paste0(work_path, "/tree_analysis/tree_qc_pass.nwk"))

B1 <- read.table(paste0(work_path, "/tree_analysis/edit_rate/e3v5v6.B1_tape_consensus.tsv.gz"), header = TRUE)
B2 <- read.table(paste0(work_path, "/tree_analysis/edit_rate/e3v5v6.B2_tape_consensus.tsv.gz"), header = TRUE)

tree_tips_B1 <- tree$tip.label[tree$tip.label %in% B1$cell_id]
tree_tips_B2 <- tree$tip.label[tree$tip.label %in% B2$cell_id]


### 129 cell types overlapped between A and B
### 25 major trajectories overlapped between A and B
celltype_B1 = unique(cell_meta$celltype[cell_meta$cell_id %in% tree_tips_B1])
celltype_B2 = unique(cell_meta$celltype[cell_meta$cell_id %in% tree_tips_B2])
celltype_common = celltype_B1[celltype_B1 %in% celltype_B2]


# ---- Function 1: extract sibling pairs and annotate ----
analyze_clade <- function(tree, cell_meta, clade_tips, clade_name) {
  subtree <- keep.tip(tree, clade_tips)
  
  tip_parents <- data.frame(
    tip_id   = 1:Ntip(subtree),
    tip_name = subtree$tip.label,
    parent   = subtree$edge[match(1:Ntip(subtree), subtree$edge[, 2]), 1]
  )
  
  sibling_pairs <- tip_parents %>%
    group_by(parent) %>%
    filter(n() >= 2) %>%
    do(as.data.frame(t(combn(.$tip_name, 2)))) %>%
    ungroup() %>%
    rename(cell_A = V1, cell_B = V2)
  
  cat(clade_name, ":", Ntip(subtree), "tips,", nrow(sibling_pairs), "sibling pairs\n")
  
  sibling_celltypes <- sibling_pairs %>%
    left_join(cell_meta %>% select(cell_id, celltype), by = c("cell_A" = "cell_id")) %>%
    rename(celltype_A = celltype) %>%
    left_join(cell_meta %>% select(cell_id, celltype), by = c("cell_B" = "cell_id")) %>%
    rename(celltype_B = celltype)
  
  list(subtree = subtree,
       sibling_pairs = sibling_pairs,
       sibling_celltypes = sibling_celltypes)
}


compute_fold_change <- function(subtree, sibling_pairs, sibling_celltypes,
                                 cell_meta, n_perm = 100, seed = 1,
                                 min_pairs = 10) {
  set.seed(seed)
  all_tips <- subtree$tip.label
  
  tip_celltype <- cell_meta$celltype[match(all_tips, cell_meta$cell_id)]
  names(tip_celltype) <- all_tips
  
  # Observed same-type pair counts
  obs_counts <- sibling_celltypes %>%
    filter(celltype_A == celltype_B, !is.na(celltype_A)) %>%
    count(celltype_A, name = "n_obs") %>%
    rename(celltype = celltype_A)
  
  # Null: shuffle cell-type labels across tips, keep sibling pair structure
  null_counts_mat <- replicate(n_perm, {
    shuffled_ct <- sample(tip_celltype)                    # shuffle the labels
    names(shuffled_ct) <- names(tip_celltype)               # keep names aligned
    ct_a <- shuffled_ct[sibling_pairs$cell_A]
    ct_b <- shuffled_ct[sibling_pairs$cell_B]
    same <- ct_a[ct_a == ct_b & !is.na(ct_a) & !is.na(ct_b)]
    table(same)
  })
  
  all_celltypes <- unique(cell_meta$celltype)
  null_mean <- sapply(all_celltypes, function(ct) {
    vals <- sapply(null_counts_mat, function(x) ifelse(ct %in% names(x), x[ct], 0))
    mean(as.numeric(vals))
  })
  
  fc_df <- data.frame(celltype = all_celltypes, n_null_mean = null_mean) %>%
    left_join(obs_counts, by = "celltype") %>%
    mutate(
      n_obs = ifelse(is.na(n_obs), 0, n_obs),
      fold_change = (n_obs + 1) / (n_null_mean + 1),
      log2_fc = log2(fold_change)
    ) %>%
    arrange(desc(log2_fc))
  
  return(fc_df)
}


# ---- Run analysis for each blastomere ----
res_B1 <- analyze_clade(tree, cell_meta, tree_tips_B1, "B1")
res_B2 <- analyze_clade(tree, cell_meta, tree_tips_B2, "B2")
res_all <- analyze_clade(tree, cell_meta, c(tree_tips_B1, tree_tips_B2), "All")

fc_B1 <- compute_fold_change(res_B1$subtree, res_B1$sibling_pairs,
                              res_B1$sibling_celltypes, cell_meta,
                              n_perm = 100)
fc_B2 <- compute_fold_change(res_B2$subtree, res_B2$sibling_pairs,
                              res_B2$sibling_celltypes, cell_meta,
                              n_perm = 100)
fc_all <- compute_fold_change(res_all$subtree, res_all$sibling_pairs,
                              res_all$sibling_celltypes, cell_meta,
                              n_perm = 100)

cat("\nB1: top 10 cell types by fold change\n")
head(fc_B1, 10)
cat("\nB2: top 10 cell types by fold change\n")
head(fc_B2, 10)
cat("\nAll: top 10 cell types by fold change\n")
head(fc_all, 10)

saveRDS(fc_B1, paste0(work_path, "/tree_analysis/sibling_cells/fc_major_trajectory_B1.rds"))
saveRDS(fc_B2, paste0(work_path, "/tree_analysis/sibling_cells/fc_major_trajectory_B2.rds"))
saveRDS(fc_all, paste0(work_path, "/tree_analysis/sibling_cells/fc_major_trajectory_All.rds"))

saveRDS(fc_B1, paste0(work_path, "/tree_analysis/sibling_cells/fc_celltype_B1.rds"))
saveRDS(fc_B2, paste0(work_path, "/tree_analysis/sibling_cells/fc_celltype_B2.rds"))
saveRDS(fc_all, paste0(work_path, "/tree_analysis/sibling_cells/fc_celltype_All.rds"))



# ---- Plot obs vs. null across major trajectories ----
fc_all = readRDS(paste0(work_path, "/tree_analysis/sibling_cells/fc_major_trajectory_All.rds"))

fc_all$celltype = factor(fc_all$celltype, levels = as.vector(fc_all$celltype))

p_compare = ggplot(fc_all, aes(x = celltype, y = fold_change, fill = celltype)) +
  geom_bar(stat="identity") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black")) +
  theme(legend.position="none") +
  scale_fill_manual(values = major_trajectory_color_plate) +
  labs(x = NULL, y = "Fold change (sibling / random)")

ggsave("~/share/Fig5_sibling_major_trajectory.pdf", p_compare, width = 5, height = 4)



# ---- Compare B1 and B2: major trajectory ----
fc_B1 = readRDS(paste0(work_path, "/tree_analysis/sibling_cells/fc_major_trajectory_B1.rds"))
fc_B2 = readRDS(paste0(work_path, "/tree_analysis/sibling_cells/fc_major_trajectory_B2.rds"))

fc_compare <- fc_B1 %>%
  select(celltype, fc_B1 = fold_change) %>%
  inner_join(fc_B2 %>% select(celltype, fc_B2 = fold_change),
             by = "celltype")

sum(fc_compare$fc_B1 >= 2 & fc_compare$fc_B2 >= 2) ### 17

cor_test <- cor.test(fc_compare$fc_B1, fc_compare$fc_B2,
                     method = "spearman")
cat("\nSpearman correlation between B1 and B2:", round(cor_test$estimate, 3),
    ", p =", format.pval(cor_test$p.value, digits = 3), "\n")
### Spearman correlation between B1 and B2: 0.939 , p = 3.54e-12

labels_df <- fc_compare %>% filter(fc_B1 > 10 & fc_B2 > 10)

p_compare <- ggplot(fc_compare, aes(x = fc_B1, y = fc_B2)) +
  geom_point(aes(color = celltype), size = 3, alpha = 0.6) +
  geom_text_repel(data = labels_df,
                  aes(label = celltype),
                  size = 3, color = "black",
                  box.padding = 0.4, max.overlaps = Inf,
                  seed = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey40") +
  theme_classic(base_size = 12) +
  scale_color_manual(values = major_trajectory_color_plate) +
  theme(axis.text = element_text(color = "black"), legend.position="none") +
  labs(x = "FC (B1)",
       y = "FC (B2)")

ggsave("~/share/FigS6_sibling_identicial_major_trajectory.pdf", p_compare, width = 5, height = 5)

# 25 major trajectories; Spearman's correlation = 0.94, p < 1e-11


# ---- Compare B1 and B2: celltype ----
fc_B1 = readRDS(paste0(work_path, "/tree_analysis/sibling_cells/fc_celltype_B1.rds"))
fc_B2 = readRDS(paste0(work_path, "/tree_analysis/sibling_cells/fc_celltype_B2.rds"))

fc_compare <- fc_B1 %>%
  select(celltype, fc_B1 = fold_change) %>%
  inner_join(fc_B2 %>% select(celltype, fc_B2 = fold_change),
             by = "celltype")

fc_compare = fc_compare %>% filter(celltype %in% celltype_common)

sum(fc_compare$fc_B1 >= 2 & fc_compare$fc_B2 >= 2) ### 82

cor_test <- cor.test(fc_compare$fc_B1, fc_compare$fc_B2,
                     method = "spearman")
cat("\nSpearman correlation between B1 and B2:", round(cor_test$estimate, 3),
    ", p =", format.pval(cor_test$p.value, digits = 3), "\n")
### Spearman correlation between B1 and B2: 0.93 , p = 4.551754e-57

major_trajectory_celltype = read.table(paste0(work_path, "/tree_analysis/major_trajectory_celltype_table.txt"), header=T, sep="\t")

fc_compare <- fc_compare %>% left_join(major_trajectory_celltype, by = "celltype")

labels_df <- fc_compare %>% filter(fc_B1 > 40 & fc_B2 > 40)

p_compare <- ggplot(fc_compare, aes(x = fc_B1, y = fc_B2)) +
  geom_point(aes(color = major_trajectory), size = 3, alpha = 0.6) +
    geom_text_repel(data = labels_df,
                  aes(label = celltype),
                  size = 3, color = "black",
                  box.padding = 0.4, max.overlaps = Inf,
                  seed = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey40") +
  theme_classic(base_size = 12) +
  scale_color_manual(values = major_trajectory_color_plate) +
  theme(axis.text = element_text(color = "black"), legend.position="none") +
  labs(x = "FC (B1)",
       y = "FC (B2)")

ggsave("~/share/FigS6_sibling_identicial_celltype.pdf", p_compare, width = 5, height = 5)

# 129 cell types; Spearman's correlation = 0.93, p < 1e-56


fc_all = readRDS(paste0(work_path, "/tree_analysis/sibling_cells/fc_major_trajectory_All.rds"))
fc_all$blastomere = "Both"
print(sum(fc_all$n_obs)/sum(fc_all$n_null_mean)) ### 3.4 fold
fc_B1 = readRDS(paste0(work_path, "/tree_analysis/sibling_cells/fc_major_trajectory_B1.rds"))
fc_B1$blastomere = "Blastomere_A"
fc_B2 = readRDS(paste0(work_path, "/tree_analysis/sibling_cells/fc_major_trajectory_B2.rds"))
fc_B2$blastomere = "Blastomere_B"
res_out = rbind(fc_all, fc_B1, fc_B2)
write.table(res_out, "~/share/FigS6_sibling_identicial_major_trajectory.txt", row.names=F, col.names=T, sep="\t", quote=F)


fc_all = readRDS(paste0(work_path, "/tree_analysis/sibling_cells/fc_celltype_All.rds"))
fc_all$blastomere = "Both"
print(sum(fc_all$n_obs)/sum(fc_all$n_null_mean)) ### 8.9 fold
fc_B1 = readRDS(paste0(work_path, "/tree_analysis/sibling_cells/fc_celltype_B1.rds"))
fc_B1$blastomere = "Blastomere_A"
fc_B2 = readRDS(paste0(work_path, "/tree_analysis/sibling_cells/fc_celltype_B2.rds"))
fc_B2$blastomere = "Blastomere_B"
res_out = rbind(fc_all, fc_B1, fc_B2)
write.table(res_out, "~/share/FigS6_sibling_identicial_celltype.txt", row.names=F, col.names=T, sep="\t", quote=F)



