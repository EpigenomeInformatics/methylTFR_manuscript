## Getting the most important TFs for each celltype
set.seed(13)

library(methylTFR)
library(ggplot2)
library(dplyr)
plot_dir <- "/icbb/projects/nitschre/methylTFR/figures"


# Import raw deviations
mtfr <- readRDS("/icbb/projects/nitschre/methylTFR/r_objects/jaspar2020_distal_deviations_.RDS")

# Get Raw deviations
mtfr_raw <- deviations(mtfr)

# Remove TEMRA
mtfr_raw <- mtfr_raw[, -which(colnames(mtfr_raw) == "51_Hf03_BlTR_Ct_WGBS_S_1.MCSv3.20170714.GRCh38.cpg.filtered.CG.bed")]

# Row-wise z-scoring
zscores <- methylTFR:::computeRowZScore(mtfr_raw)

# Row means of each celltype
zscores <- as.data.frame(zscores)

zscores$CM_mean <- rowMeans(zscores[, c(
  "51_Hf03_BlCM_Ct_WGBS_S_1.MCSv3.20170714.GRCh38.cpg.filtered.CG.bed",
  "51_Hf04_BlCM_Ct_WGBS_S_1.MCSv3.20170714.GRCh38.cpg.filtered.CG.bed"
)], na.rm = TRUE)

zscores$EM_mean <- rowMeans(zscores[, c(
  "51_Hf03_BlEM_Ct_WGBS_S_1.MCSv3.20170714.GRCh38.cpg.filtered.CG.bed",
  "51_Hf04_BlEM_Ct_WGBS_S_1.MCSv3.20170714.GRCh38.cpg.filtered.CG.bed"
)], na.rm = TRUE)

zscores$TN_mean <- rowMeans(zscores[, c(
  "51_Hf03_BlTN_Ct_WGBS_S_1.MCSv3.20170714.GRCh38.cpg.filtered.CG.bed",
  "51_Hf04_BlTN_Ct_WGBS_S_1.MCSv3.20170714.GRCh38.cpg.filtered.CG.bed"
)], na.rm = TRUE)

# Remove celltypes
zscores <- zscores[, -c(1:6)]

# Differences TN - TEM and TN - TCM 
zscores$TN_TEM_diff <- zscores$TN_mean - zscores$EM_mean
zscores$TN_TCM_diff <- zscores$TN_mean - zscores$CM_mean

# Add pvalues from diff test
diff_cm <- readRDS("//icbb/projects/nitschre/methylTFR/r_objects/diff_cm.RDS")
diff_em <- readRDS("//icbb/projects/nitschre/methylTFR/r_objects/diff_em.RDS")

zscores$diff_cm <- diff_cm$p_value_adjusted
zscores$diff_em <- diff_em$p_value_adjusted

# Find sign TFs for each group
zscores <- zscores %>%
  mutate(group = case_when(
    abs(zscores$TN_TEM_diff) >= 0.5 & diff_em < 0.05 ~ "signif_TEM",
    abs(zscores$TN_TCM_diff) >= 0.5 & diff_cm < 0.05 ~ "signif_TCM",
    TRUE ~ "other"
  ))

# Correlation coefficient
cor <- cor(zscores$TN_TCM_diff, zscores$TN_TEM_diff, method="spearman")

# Scatterplot
fn_scatterplot <- file.path(plot_dir, "scatterplot.pdf")
pdf(fn_scatterplot)
ggplot(zscores, aes(x=TN_TCM_diff, y=TN_TEM_diff, color=group)) +
  geom_point(alpha=0.4, size=1) +
  geom_text(
    data = subset(zscores, group == "signif_TEM" | group == "signif_TCM"),
    aes(label = rownames(subset(zscores, group == "signif_TEM" | group == "signif_TCM"))),
    vjust = -0.5, size = 2) +
  theme_minimal()+
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("signif_TEM" = "red", "signif_TCM" = "blue", "other" = "grey")) +
  annotate("text", x = -1, y = 2, label = paste0("Cor. Coef.:",round(cor, 1), size = 6))
dev.off()



### Using samples that were merged before running methyltfr
# Import raw deviations
mtfr <- readRDS("/icbb/projects/nitschre/methylTFR/r_objects/jaspar2020_distal_deviations_merged.RDS")

# Get Raw deviations
mtfr_raw <- deviations(mtfr)

# Remove TEMRA
mtfr_raw <- mtfr_raw[, -which(colnames(mtfr_raw) == "TEMRA")]

# Row-wise z-scoring
zscores <- methylTFR:::computeRowZScore(mtfr_raw)

zscores <- as.data.frame(zscores)

# Differences TN - TEM and TN - TCM 
zscores$TN_TEM_diff <- zscores$TN - zscores$TEM
zscores$TN_TCM_diff <- zscores$TN - zscores$TCM

# Add pvalues from diff test
diff_cm <- readRDS("//icbb/projects/nitschre/methylTFR/r_objects/diff_cm.RDS")
diff_em <- readRDS("//icbb/projects/nitschre/methylTFR/r_objects/diff_em.RDS")

zscores$diff_cm <- diff_cm$p_value_adjusted
zscores$diff_em <- diff_em$p_value_adjusted

# Find sign TFs for each group
zscores <- zscores %>%
  mutate(group = case_when(
    abs(zscores$TN_TEM_diff) >= 0.5 & diff_em < 0.05 ~ "signif_TEM",
    abs(zscores$TN_TCM_diff) >= 0.5 & diff_cm < 0.05 ~ "signif_TCM",
    TRUE ~ "other"
  ))

# Correlation coefficient
cor <- cor(zscores$TN_TCM_diff, zscores$TN_TEM_diff, method="spearman")

# Scatterplot
fn_scatterplot <- fn_pca_ind <- file.path(plot_dir, "scatterplot_merged.pdf")

labels <- subset(zscores, group %in% c("signif_TEM", "signif_TCM") 
                  & diff_em < 0.05 | diff_cm < 0.05)

pdf(fn_scatterplot)
ggplot(zscores, aes(x=TN_TCM_diff, y=TN_TEM_diff, color=group)) +
  geom_point(alpha=0.6, size=1) +
  geom_text_repel(
    data = labels,
    aes(label = rownames(labels)),
    max.overlaps = 180,
    size = 2) +
  theme_minimal()+
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("signif_TEM" = "red", "signif_TCM" = "blue", "other" = "grey")) +
  annotate("text", x = -1, y = 2, label = paste0("Cor. Coef.:",round(cor, 1), size = 6))
dev.off()
# 92 missing

#
fn_scatterplot <- fn_pca_ind <- file.path(plot_dir, "scatterplot_merged.pdf")

labels <- subset(zscores, group %in% c("signif_TEM", "signif_TCM") 
                  & diff_em < 0.05 | diff_cm < 0.05)

pdf(fn_scatterplot)
ggplot(zscores, aes(x=TN_TCM_diff, y=TN_TEM_diff, color=group)) +
  geom_point(alpha=0.6, size=1) +
  geom_text(
    data = labels,
    aes(label = rownames(labels)),
    size = 2) +
  theme_minimal()+
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("signif_TEM" = "red", "signif_TCM" = "blue", "other" = "grey")) +
  annotate("text", x = -1, y = 2, label = paste0("Cor. Coef.:",round(cor, 1), size = 6))
  dev.off()