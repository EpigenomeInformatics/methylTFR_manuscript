## Getting the most important TFs for each celltype
set.seed(13)

library(methylTFR)
library(ggplot2)
library(ggrepel)
library(dplyr)
plot_dir <- "/icbb/projects/nitschre/methylTFR/figures"


# Import raw deviations
mtfr <- readRDS("/icbb/projects/nitschre/methylTFR/r_objects/jaspar2020_distal_deviations.RDS")

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
    abs(zscores$TN_TEM_diff) >= 2 & diff_em < 0.035 ~ "signif",
    abs(zscores$TN_TCM_diff) >= 1.5 & diff_cm < 0.035 ~ "signif",
    TRUE ~ "other"
  ))

# Correlation coefficient
cor <- cor(zscores$TN_TCM_diff, zscores$TN_TEM_diff, method="spearman")

# Scatterplot
labels <- subset(zscores, group == "signif")

fn_scatterplot <- file.path(plot_dir, "scatterplot_003.pdf")
pdf(fn_scatterplot)
ggplot(zscores, aes(x=TN_TCM_diff, y=TN_TEM_diff, color=group)) +
  geom_point(alpha=0.6, size=2) +
  geom_text_repel(
    data = labels,
    aes(label = rownames(labels)),
    max.overlaps = 100,
    size = 2,
    segment.color = NA) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("signif" = "red", "other" = "grey")) +
  annotate("text", x = -1, y = 2, label = paste0("Cor. Coef.:",round(cor, 1), size = 6))
dev.off()

#### With legend

p <-ggplot(zscores, aes(x=TN_TCM_diff, y=TN_TEM_diff, color=group)) +
  geom_point(alpha=0.4, size=1) +
  geom_text_repel(
    data = labels,
    aes(label = rownames(labels)),
    max.overlaps = 200,
    size = 2.5,
    segment.color = NA) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("signif" = "red", "other" = "grey")) +
  annotate("text", x = -1, y = 2, label = paste0("Cor. Coef.:",round(cor, 1), size = 6))

# All TFs are plotted at the right of the plot
df_tf <- subset(zscores, group == "signif" & TN_TEM_diff > 0)
df_labels <- data.frame(
  tf = rownames(df_tf),
  y = rev(seq(1, 200, length.out =length(rownames(df_tf))))
)

p_legend <- ggplot(df_labels, aes(x=0, y=y, label = tf)) +
                geom_text(hjust=0, size= 2)+
                theme_void() +
                labs(x = NULL, y = NULL) +
                ylim(0,200)



fn_scatterplot <- file.path(plot_dir, "scatterplot_legend.pdf")
pdf(fn_scatterplot, height = 8)
p + p_legend + plot_layout(widths = c(3,1))
dev.off()

