set.seed(12)

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

    ### Linear model
    res <- lm(as.formula(paste(id,"~",colnames(matrix)[3],
        "+",colnames(matrix)[4], 
        "+",colnames(matrix)[5],
        "+","tracking_id")),
        data=matrix)

    return(res)
})