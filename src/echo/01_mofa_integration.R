#!/usr/bin/env Rscript

#####################################################################
# 01_mofa_integration.R
# Created by IBG on 25-10-2025
# Script to perform TF activity integration using MOFA2
#####################################################################

suppressPackageStartupMessages({
  library(MOFA2)
  library(reticulate)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(viridis)
  library(tibble)
  library(patchwork)
  library(ggrepel)
  library(scales)
})
set.seed(12)
use_python(Sys.which("python"), required = TRUE)

# View display names and colours, matching the manuscript figure. The
# object names the views chromVar and mtfr; the legend spells them out.
view_labels <- c("chromVar" = "chromVAR", "mtfr" = "methylTFR")
view_colors <- c("chromVAR" = "#8FA8D4", "methylTFR" = "#E4675C")

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


# Paths. The model and the intermediate tables are written under the
# analysis directory, which is writable, rather than into a shared
# project folder.
analysis_dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/echo"
plot_dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/figures/echo/mofa_integration"
r_objects_dir <- file.path(analysis_dir, "mofa_230826")
table_dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/tables"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(r_objects_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# Load pseudobulk TF activity matrices
chromVar <- readRDS(gzfile("/icbb/projects/nitschre/methylTFR/scripts/other/cvar_zscores_adj_psuedobulk.R"))
mtfr <- readRDS(gzfile("/icbb/projects/nitschre/methylTFR/scripts/other/mtfr_zscores_adj_psuedobulk.R"))

# Quick sanity checks
stopifnot(is.matrix(chromVar) | is.data.frame(chromVar))
stopifnot(is.matrix(mtfr) | is.data.frame(mtfr))
if (!all(colnames(chromVar) == colnames(mtfr))) {
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
if (file.exists(outfile)) {
  MOFAobject.trained <- load_model(outfile)
} else {
  MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk = FALSE)
}

# Save trained object in rds
saveRDS(MOFAobject.trained, file.path(plot_dir, "MOFAobject_trained.rds"))

# Extract factors and compute ANOVA R2 for celltype.
# get_factors returns one row per sample AND factor, so the cell type is
# joined on the sample rather than assigned as a per-sample vector.
factors_long <- get_factors(MOFAobject.trained, as.data.frame = TRUE)
factors_long$sample <- as.character(factors_long$sample)
factors_long <- left_join(factors_long, metadata, by = "sample")
if (any(is.na(factors_long$celltype))) {
  stop("Some factor rows could not be matched to a cell type")
}

r2_df <- data.frame(Factor = character(), R2 = numeric(), stringsAsFactors = FALSE)

for (fa in unique(factors_long$factor)) {
  df <- subset(factors_long, factor == fa)
  mod <- aov(value ~ celltype, data = df)
  ss_tab <- summary(mod)[[1]]
  ss_factor <- ss_tab["celltype", "Sum Sq"]
  ss_resid <- ss_tab["Residuals", "Sum Sq"]
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

if (!(f_x %in% colnames(factors_wide)) | !(f_y %in% colnames(factors_wide))) {
  # fallback to numeric factor names if needed
  available_factors <- setdiff(colnames(factors_wide), "celltype")
  f_x <- available_factors[1]
  f_y <- available_factors[2]
}

p_factors_scatter <- ggplot(factors_wide, aes(x = .data[[f_x]], y = .data[[f_y]], color = celltype)) +
  geom_point(alpha = 0.9, size = 2) +
  stat_ellipse(aes(group = celltype), linetype = 2, alpha = 0.5, linewidth = 0.4) +
  scale_color_manual(values = cell_type_colors) + # Use custom colors
  labs(
    x = paste0("Factor ", f_x), y = paste0("Factor ", f_y),
    title = paste0("MOFA: Factor ", f_x, " vs ", f_y, " (colored by celltype)")
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "right")

ggsave(
  filename = file.path(plot_dir, "factor_scatter_Fx_Fy.pdf"), plot = p_factors_scatter,
  width = 7, height = 5
)

# Variance explained per view × factor (from MOFA)
ve <- get_variance_explained(MOFAobject.trained)
# r2_per_factor is factors in the rows and views in the columns, which is
# also how the second variance panel further down reads it
r2_mat <- as.data.frame(ve$r2_per_factor[[1]])
r2_mat$factor <- rownames(r2_mat)
r2_long <- r2_mat %>% pivot_longer(cols = -factor, names_to = "view", values_to = "r2")

top_factors <- head(r2_df$Factor, 5)
r2_long <- r2_long %>% filter(factor %in% top_factors)
if (nrow(r2_long) == 0) {
  stop(
    "No variance explained rows left. Factors in the model: ",
    paste(utils::head(rownames(r2_mat), 5), collapse = ", "),
    "; requested: ", paste(top_factors, collapse = ", ")
  )
}
r2_long$factor <- factor(r2_long$factor, levels = top_factors)

# Modality contribution per factor, with the variance the cell type
# explains drawn over it. One panel rather than two, as in the Blueprint
# figure, so the bars and the R2 line share an x axis.
weights_df <- get_weights(MOFAobject.trained, as.data.frame = TRUE)
weights_df <- weights_df %>% mutate(value_abs = abs(value), value_signed = value)

agg <- weights_df %>%
  filter(factor %in% top_factors) %>%
  group_by(factor, view) %>%
  summarise(sum_abs = sum(value_abs, na.rm = TRUE), .groups = "drop") %>%
  group_by(factor) %>%
  mutate(frac = sum_abs / sum(sum_abs)) %>%
  ungroup()

unmapped_views <- setdiff(unique(agg$view), names(view_labels))
if (length(unmapped_views) > 0) {
  stop("No display name for view(s): ", paste(unmapped_views, collapse = ", "))
}
agg$view <- factor(unname(view_labels[as.character(agg$view)]), levels = names(view_colors))
agg$factor <- factor(agg$factor, levels = top_factors)

r2_line <- r2_df[r2_df$Factor %in% top_factors, ]
r2_line$Factor <- factor(r2_line$Factor, levels = top_factors)

p_modality_frac <- ggplot(agg, aes(x = factor, y = frac)) +
  geom_bar(aes(fill = view), stat = "identity", width = 0.7) +
  geom_text(aes(group = view, label = ifelse(frac > 0.03, paste0(round(frac * 100, 1), "%"), "")),
    position = position_stack(vjust = 0.5), size = 3
  ) +
  geom_line(data = r2_line, aes(x = Factor, y = R2, group = 1, colour = "Variance Explained (R2)")) +
  geom_point(data = r2_line, aes(x = Factor, y = R2, colour = "Variance Explained (R2)")) +
  scale_colour_manual(values = c("Variance Explained (R2)" = "purple"), name = NULL) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    sec.axis = sec_axis(~., name = expression(Variance ~ Explained ~ (R^2)))
  ) +
  scale_fill_manual(values = view_colors) +
  labs(x = "Factor", y = "Relative contribution (|weights| fraction)", fill = "View") +
  theme_classic(base_size = 12) +
  theme(axis.title.y.right = element_text(colour = "purple"))

ggsave(
  filename = file.path(plot_dir, "modality_contribution_fraction_R2.pdf"),
  plot = p_modality_frac, width = 8, height = 4.5
)

# The aggregate behind the panel, so 04 can redraw it without retraining
write.csv(agg, file.path(table_dir, "echo_mofa_modality_contribution.csv"), row.names = FALSE)
write.csv(r2_df, file.path(table_dir, "echo_mofa_factor_celltype_R2.csv"), row.names = FALSE)

# Top weighted features per factor, kept for the tables rather than a dotplot
dot_topN <- 6
top_by_factor <- weights_df %>%
  group_by(factor) %>%
  arrange(desc(value_abs)) %>%
  slice_head(n = dot_topN) %>%
  ungroup()

# The ECHO top TFs, written out so the footprint scripts can use this
# dataset's own list instead of the Blueprint one
write.csv(
  top_by_factor[, c("feature", "factor", "view", "value", "value_abs")],
  file.path(table_dir, "mofa_echo_topTFs.csv"),
  row.names = FALSE
)

# The training diagnostics on one page. The manuscript figure is
# assembled in 04, which adds the SELEX panel and the paired heatmaps.
combined_fig <- (p_factors_scatter / p_modality_frac) + plot_annotation(tag_levels = "A")
ggsave(filename = file.path(plot_dir, "figure_combined_panels.pdf"), plot = combined_fig, width = 10, height = 9)

# Save intermediate tables for reproducibility
write.csv(agg, file.path(r_objects_dir, "modality_contribution_by_factor.csv"), row.names = FALSE)
write.csv(weights_df, file.path(r_objects_dir, "MOFA_weights_long.csv"), row.names = FALSE)
write.csv(r2_df, file.path(r_objects_dir, "factor_celltype_R2.csv"), row.names = FALSE)

message("All done. Figures saved to: ", plot_dir, "  |  R objects / tables saved to: ", r_objects_dir)

# Extract factors, cell type joined on the sample as above
factors_long <- get_factors(MOFAobject.trained, as.data.frame = TRUE)
factors_long$sample <- as.character(factors_long$sample)
factors_long <- left_join(factors_long, metadata, by = "sample")

# Filter to top factors (optional)
top_factors <- c("Factor1", "Factor2", "Factor3", "Factor4", "Factor7")
factors_long <- factors_long %>% filter(factor %in% top_factors)

# Plot: strip / dot plot per factor
p <- ggplot(factors_long, aes(x = factor, y = value, color = celltype)) +
  geom_jitter(width = 0.2, height = 0, size = 2, alpha = 0.8) +
  scale_color_manual(values = cell_type_colors) + # Use custom colors
  labs(
    x = "Factor", y = "Factor value", color = "Cell type",
    title = "MOFA factor values per sample by cell type"
  ) +
  theme_classic(base_size = 13) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
ggsave(
  filename = file.path(plot_dir, "factors_stripplot_by_celltype.pdf"), plot = p,
  width = 8, height = 5
)

#####################################################################
