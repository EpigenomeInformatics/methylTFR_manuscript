#!/usr/bin/env Rscript

#####################################################################
# mofa_integration.R
# Created by IBG on 25-10-2025
# Script to perform TF activity integration using MOFA2
#####################################################################
suppressPackageStartupMessages({
library(MOFA2)
library(dplyr)
library(tidyr)
library(viridis)
library(tibble)
library(patchwork)
library(ggrepel)
library(scales)})

set.seed(12)

# Define custom colors for cell types
cell_type_colors <- c(
  "B-cell" = "#AE017E",
  "Monocyte" = "#CC4C02",
  "NK-cell" = "#A65628",
  "Th-Mem" = "#41B6C4",
  "Tc-Mem" = "#4292C6",
  "Tc-Naive" = "#888FB5",
  "Th-Naive" = "#C7E9B4"
)


# Paths
plot_dir <- "/icbb/projects/igunduz/methylTFR/figures/MOFA/"
r_objects_dir <- "/icbb/projects/nitschre/methylTFR/r_objects/"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(r_objects_dir, recursive = TRUE, showWarnings = FALSE)

# Load pseudobulk TF activity matrices 
chromVar <- readRDS(gzfile("/icbb/projects/nitschre/methylTFR/scripts/other/cvar_zscores_adj_psuedobulk.R"))
mtfr     <- readRDS(gzfile("/icbb/projects/nitschre/methylTFR/scripts/other/mtfr_zscores_adj_psuedobulk.R"))

# Quick sanity checks
stopifnot(is.matrix(chromVar) | is.data.frame(chromVar))
stopifnot(is.matrix(mtfr) | is.data.frame(mtfr))
if(!all(colnames(chromVar) == colnames(mtfr))){
  warning("Column names (samples) do not match exactly between chromVar and mtfr.")
}

# Create MOFA object and metadata
data_list <- list(chromVar = as.matrix(chromVar), mtfr = as.matrix(mtfr))
MOFAobject <- create_mofa(data_list)

metadata <- data.frame(
  sample = colnames(mtfr),
  celltype = sub("_.*", "", colnames(mtfr)),
  stringsAsFactors = FALSE
)
samples_metadata(MOFAobject) <- metadata


# Prepare and train MOFA (save HDF5)
MOFAobject <- prepare_mofa(MOFAobject)
outfile <- file.path(r_objects_dir, "model.hdf5")

# If model file exists, load trained object; otherwise run training
if(file.exists(outfile)){
  MOFAobject.trained <- load_model(outfile)
} else {
  MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk = TRUE)
}

# Save trained object in rds 
saveRDS(MOFAobject.trained, file.path(plot_dir, "MOFAobject_trained.rds"))

# Extract factors and compute ANOVA R2 for celltype
factors_long <- get_factors(MOFAobject.trained, as.data.frame = TRUE)
factors_long$celltype <- metadata$celltype

r2_df <- data.frame(Factor = character(), R2 = numeric(), stringsAsFactors = FALSE)

for(fa in unique(factors_long$factor)){
  df <- subset(factors_long, factor == fa)
  mod <- aov(value ~ celltype, data = df)
  ss_tab <- summary(mod)[[1]]
  ss_factor <- ss_tab["celltype", "Sum Sq"]
  ss_resid  <- ss_tab["Residuals", "Sum Sq"]
  r2_val <- ss_factor / (ss_factor + ss_resid)
  r2_df <- rbind(r2_df, data.frame(Factor = fa, R2 = r2_val, stringsAsFactors = FALSE))
}
r2_df <- r2_df %>% arrange(desc(R2))
write.csv(r2_df, file.path(plot_dir, "r2Values.csv"), row.names = FALSE)

# Choose top factors (by R2) for plotting
top_n <- 5
top_factors <- head(r2_df$Factor, top_n)

# Plot 1: Factor scatter (example F1 vs F2) colored by celltype
factors_wide <- factors_long %>%
  pivot_wider(names_from = factor, values_from = value)

# default pair: top two factors by R2 (if <2 fallback to 1,2)
f_x <- ifelse(length(top_factors) >= 1, top_factors[1], "1")
f_y <- ifelse(length(top_factors) >= 2, top_factors[2], "2")

if(!(f_x %in% colnames(factors_wide)) | !(f_y %in% colnames(factors_wide))){
  # fallback to numeric factor names if needed
  available_factors <- setdiff(colnames(factors_wide), "celltype")
  f_x <- available_factors[1]; f_y <- available_factors[2]
}

p_factors_scatter <- ggplot(factors_wide, aes_string(x = f_x, y = f_y, color = "celltype")) +
  geom_point(alpha = 0.9, size = 2) +
  stat_ellipse(aes(group = celltype), linetype = 2, alpha = 0.5, size = 0.4) +
  scale_color_manual(values = cell_type_colors) +  # Use custom colors
  labs(x = paste0("Factor ", f_x), y = paste0("Factor ", f_y),
       title = paste0("MOFA: Factor ", f_x, " vs ", f_y, " (colored by celltype)")) +
  theme_classic(base_size = 13) +
  theme(legend.position = "right")

ggsave(filename = file.path(plot_dir, "factor_scatter_Fx_Fy.pdf"), plot = p_factors_scatter,
       width = 7, height = 5)

# Variance explained per view × factor (from MOFA)
ve <- get_variance_explained(MOFAobject.trained)
# Convert to long table (robust for MOFA2 output)
r2_mat <- as.data.frame(ve$r2_per_factor[[1]])
r2_mat$view <- rownames(r2_mat)
r2_long <- r2_mat %>% pivot_longer(cols = -view, names_to = "factor", values_to = "r2")

# Keep only top factors (from R2 ranking)
top_factors <- head(r2_df$Factor, 5)
r2_long <- r2_long %>% filter(view %in% top_factors)
r2_long$factor <- factor(r2_long$view, levels = top_factors)

# Plot
pvar <- ggplot(r2_long, aes(x = factor, y = r2, fill = view)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_viridis_d(option = "C") +
  labs(x = "Factor", y = "Variance explained (R²)", fill = "View",
       title = "Variance explained per factor by view (top factors)") +
  theme_classic(base_size = 12)

ggsave(filename = file.path(plot_dir, "variance_explained_topFactors.pdf"), plot = pvar,
       width = 8, height = 4)
# Get weights and prepare modality contribution plots
weights_df <- get_weights(MOFAobject.trained, as.data.frame = TRUE)
weights_df <- weights_df %>% mutate(value_abs = abs(value), value_signed = value)

# Aggregate absolute weights per view × factor and compute fraction per factor
agg <- weights_df %>%
  filter(factor %in% top_factors) %>%
  group_by(factor, view) %>%
  summarise(sum_abs = sum(value_abs, na.rm = TRUE), .groups = "drop") %>%
  group_by(factor) %>%
  mutate(frac = sum_abs / sum(sum_abs)) %>%
  ungroup()

agg$factor <- factor(agg$factor, levels = top_factors)

# Plot stacked contribution (fraction)
p_modality_frac <- ggplot(agg, aes(x = factor, y = frac, fill = view)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = ifelse(frac > 0.03, paste0(round(frac*100,1), "%"), "")),
            position = position_stack(vjust = 0.5), size = 3) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_viridis_d(option = "C") +
  labs(x = "Factor", y = "Relative contribution (|weights| fraction)", fill = "View",
       title = "Modality contribution per factor (weights-based)") +
  theme_classic(base_size = 12)

ggsave(filename = file.path(plot_dir, "modality_contribution_fraction_topFactors.pdf"), plot = p_modality_frac,
       width = 8, height = 4)

# Also plot absolute sums (dodged), with optional log-scale if large dynamic range
p_modality_abs <- ggplot(agg, aes(x = factor, y = sum_abs, fill = view)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_text(aes(label = round(sum_abs, 1)), position = position_dodge(width = 0.8), vjust = -0.5, size = 3) +
  scale_y_continuous(trans = "log10", labels = scales::comma_format()) +
  scale_fill_viridis_d(option = "C") +
  labs(x = "Factor", y = "Sum of |weights| (log10 scale)", fill = "View",
       title = "Absolute modality contribution per factor (sum |weights|)") +
  theme_classic(base_size = 12)

ggsave(filename = file.path(plot_dir, "modality_contribution_abs_topFactors_log.pdf"), plot = p_modality_abs,
       width = 9, height = 4)

# Top TF dotplot (signed + magnitude)
dot_topN <- 6
top_by_factor <- weights_df %>%
  group_by(factor) %>%
  arrange(desc(value_abs)) %>%
  slice_head(n = dot_topN) %>%
  ungroup() %>%
  mutate(feature = factor(feature, levels = unique(feature)))

p_dot <- ggplot(top_by_factor, aes(x = factor, y = feature, size = value_abs, color = value_signed)) +
  geom_point(alpha = 0.9) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(min(top_by_factor$value_signed), max(top_by_factor$value_signed))) +
  scale_size_continuous(range = c(2, 6)) +
  labs(x = "Factor", y = "TF (feature)", color = "Signed weight", size = "|weight|",
       title = paste0("Top ", dot_topN, " TFs per factor (signed weights)")) +
  theme_classic(base_size = 11) +
  theme(axis.text.y = element_text(size = 8))

ggsave(filename = file.path(plot_dir, "topTFs_dotplot_signed.pdf"), plot = p_dot, width = 10, height = 6)

# Combine panels into a single figure for main text (example)
# Use p_factors_scatter, p_modality_frac, p_variance_expl, p_dot
combined_fig <- (p_factors_scatter + p_modality_frac) / (p_variance_expl + p_dot) + plot_annotation(tag_levels = "A")
ggsave(filename = file.path(plot_dir, "figure_combined_panels.pdf"), plot = combined_fig, width = 12, height = 10)

# Save intermediate tables for reproducibility
write.csv(agg, file.path(r_objects_dir, "modality_contribution_by_factor.csv"), row.names = FALSE)
write.csv(weights_df, file.path(r_objects_dir, "MOFA_weights_long.csv"), row.names = FALSE)
write.csv(r2_df, file.path(r_objects_dir, "factor_celltype_R2.csv"), row.names = FALSE)

message("All done. Figures saved to: ", plot_dir, "  |  R objects / tables saved to: ", r_objects_dir)

# Extract factors
factors_long <- get_factors(MOFAobject.trained, as.data.frame = TRUE)
factors_long$celltype <- metadata$celltype  # ensure metadata aligned

# Filter to top factors (optional)
top_factors <- c("Factor1","Factor2","Factor3","Factor4","Factor7") 
factors_long <- factors_long %>% filter(factor %in% top_factors)

# Plot: strip / dot plot per factor
p <- ggplot(factors_long, aes(x = factor, y = value, color = celltype)) +
  geom_jitter(width = 0.2, height = 0, size = 2, alpha = 0.8) +
  scale_color_manual(values = cell_type_colors) +  # Use custom colors
  labs(x = "Factor", y = "Factor value", color = "Cell type",
       title = "MOFA factor values per sample by cell type") +
  theme_classic(base_size = 13) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())
ggsave(filename = file.path(plot_dir, "factors_stripplot_by_celltype.pdf"), plot = p,
       width = 8, height = 5)

# Total variance explained per view
r2_per_factor_df <- get_variance_explained(MOFAobject.trained)$r2_per_factor$group1
r2_per_factor_df <- as.data.frame(r2_per_factor_df)

# Add Factor as a column
r2_per_factor_df <- r2_per_factor_df %>%
  rownames_to_column(var = "Factor")

# Convert to long format and use fractions
r2_long <- r2_per_factor_df %>%
  pivot_longer(cols = -Factor, names_to = "View", values_to = "R2_percent") %>%
  mutate(R2_fraction = R2_percent / 100)

# Filter for the top factors (e.g., Factor 1-7, where contribution drops off)
top_n_factors <- 7
top_factors_levels <- paste0("Factor", 1:top_n_factors)

r2_long_top <- r2_long %>%
  filter(Factor %in% top_factors_levels) %>%
  mutate(Factor = factor(Factor, levels = top_factors_levels))

# Define custom colors for the views
view_colors <- c("green", "blue")

# Plot 2: Variance Explained per Factor by View
pvar_per_factor <- ggplot(r2_long_top, aes(x = Factor, y = R2_fraction, fill = View)) +
  geom_bar(stat = "identity", position = "dodge") +  # Use dodge for side-by-side bars
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = view_colors) +  # Use custom green and blue colors
  labs(x = "Factor", y = "Variance Explained (R²)", fill = "View",
       title = paste0("Variance Explained per Factor by View (Top ", top_n_factors, " Factors)")) +
  theme_classic(base_size = 14)

ggsave(filename = file.path(plot_dir, "variance_explained_per_factor_topFactors.pdf"), plot = pvar_per_factor,
       width = 9, height = 5)

# Plot a correlation for single top TF between chromVar and mtfr (e.g., BATF)
chromVar_BATF <- chromVar["BATF",]
mtfr_BATF <- mtfr["BATF",]
data_scatter <- data.frame(
  chromVar = chromVar_BATF,
  mtfr = mtfr_BATF,
  celltype = metadata$celltype
)

p_data_scatter <- ggplot(data_scatter, aes(x = chromVar, y = mtfr, color = celltype)) +
  geom_point() +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  scale_color_manual(values = cell_type_colors) +
  theme_classic(base_size = 13) +
  labs(title = "TF Activity Correlation (chromVar vs mtfr) for BATF")
ggsave(filename = file.path(plot_dir,"BATF_chromVar_vs_mtfr_scatter.pdf"),
 plot = p_data_scatter, width = 6, height = 5)


