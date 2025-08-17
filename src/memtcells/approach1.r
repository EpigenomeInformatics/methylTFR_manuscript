set.seed(12)
library(glmnet)
library(ggplot2)


# Load Diff TFs between TN and TEM
diff <- readRDS("/icbb/projects/nitschre/methylTFR/r_objects/diff_em.RDS")

# pvalue < 0.05
diff <- subset(diff, p_value_adjusted < 0.05)

# Find TFs that are associated with certain gene
genes = c("FOXP1","BCL6","E2F2","RUNX3")

results <- lapply(genes, function(gene){
    source("/icbb/projects/nitschre/methylTFR/scripts/other/06_rna_process.R")
    TFs <- subset(dorothea_with_tracking, target == gene)$tf
    id <- subset(dorothea_with_tracking,target==gene)$tracking_id[1]

    # Keep TFs that are also in Diff 
    diffTFs <- intersect(TFs, diff$motifs)

    # Expressionsmatrix: keep only for certain gene and transpose it
    expr_matrix_log <- t(subset(expr_matrix_log, tracking_id == id))
    colnames(expr_matrix_log) <- expr_matrix_log[1,]
    expr_matrix_log <- expr_matrix_log[-1, , drop = FALSE]

    expr_matrix_log <- as.data.frame(expr_matrix_log)
    expr_matrix_log$tracking_id <- rownames(expr_matrix_log)


    # Transpose deviation matrix
    deviation_scores <- t(subset(deviation_scores, rownames(deviation_scores) %in% diffTFs))
    deviation_scores <- as.data.frame(deviation_scores)
    deviation_scores$tracking_id <- rownames(deviation_scores)

    # Join both tables
    matrix <- full_join(expr_matrix_log, deviation_scores, by="tracking_id")

    # Model
    # Step 1: Setup
    gene_id <- colnames(matrix)[1]  # target gene
    tf_cols <- setdiff(colnames(matrix), c(gene_id, "tracking_id"))

    # Step 2: Prepare data frame with all TFs
    df <- matrix[, c(gene_id, tf_cols)]

    # Step 3: Remove columns with all NA (some TFs may have only NAs)
    df <- df[, !sapply(df, function(x) all(is.na(x)))]

    # Step 4: Build formula with main effects + pairwise TF:TF interactions
    all_tfs <- setdiff(colnames(df), gene_id)
    interaction_terms <- combn(all_tfs, 2, function(x) paste0(x[1], ":", x[2]))
    full_formula <- as.formula(
    paste(gene_id, "~", paste(c(all_tfs, interaction_terms), collapse = " + "))
    )

    # Step 5: Build model matrix and response vector
    options(na.action='na.pass') # keep NA rows
    X <- model.matrix(full_formula, data = df)
    y <- as.numeric(df[[gene_id]])

    # Step 6: Drop rows with NA in y or X
    valid_rows <- complete.cases(X, y)
    X <- X[valid_rows, ]
    y <- y[valid_rows]

    # Step 7: Fit ridge model
    fit <- cv.glmnet(X, y, alpha = 0)

    # Step 8: Extract non-zero coefficients
    coef_mat <- coef(fit, s = "lambda.min")
    selected <- as.matrix(coef_mat)
    selected_nonzero <- selected[selected[, 1] != 0, , drop = FALSE]

    # Step 9: Print non-zero effects
    print(selected_nonzero)

    # Step 10: Plot selected effects (excluding intercept)
    df_plot <- data.frame(
        term = rownames(selected_nonzero)[-1],
        coef = selected_nonzero[-1, 1]
    )

    p <- ggplot(df_plot, aes(x = reorder(term, coef), y = coef)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        labs(title = paste("TF and TF:TF interactions predictive of", gene_id),
            x = "TF or Interaction",
            y = "Coefficient") +
    theme_minimal()

    ggsave(paste0("/icbb/projects/nitschre/methylTFR/figures/memoryTcells/Model_coefficients_",gene, ".pdf"), plot=p)
    return(selected_nonzero)
})