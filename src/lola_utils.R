lolaVolcanoPlot <- function(lolaRes,outputDir,comparison,region,
database = "TF_motif_clusters",signifCol = "qValue"){

df <- lolaRes$region[[comparison]][[region]] %>%
  dplyr::filter(collection == database) %>%
  dplyr::filter(userSet %in% c("rankCut_1000_hyper","rankCut_1000_hypo")) %>%
  dplyr::mutate(condition = ifelse(userSet == "rankCut_1000_hyper","gain","loss"))
  
  if (signifCol == "qValue") {
    signifCol <- "qValueLog"
    df$qValueLog <- -log10(df$qValue)
    df$differential <- ifelse(df$qValueLog > 1.3, "Differential", "Non-differential")
  }
oddsRatioCol <- "log2OR"
df$log2OR <- log2(df$oddsRatio)
# Get top 10 motifs for gain
top_motifs_gain <- df[df$condition == "gain", ]
top_motifs_gain <- top_motifs_gain[top_motifs_gain$differential == "Differential", ]
top_motifs_gain <- top_motifs_gain[order(-top_motifs_gain[[signifCol]]), ]
top_motifs_gain <- head(top_motifs_gain, n = 10)
# Select the description and conditions column
top_motifs_gain <- dplyr::select(top_motifs_gain, description, condition,log2OR,qValueLog)


# Get top 10 motifs for loss
top_motifs_loss <- df[df$condition == "loss", ]
top_motifs_loss <- top_motifs_loss[top_motifs_loss$differential == "Differential", ]
top_motifs_loss <- top_motifs_loss[order(-top_motifs_loss[[signifCol]]), ]
top_motifs_loss <- head(top_motifs_loss, n = 10)
top_motifs_loss <- dplyr::select(top_motifs_loss, description, condition,log2OR,qValueLog)
top_motifs_loss$log2OR <- -top_motifs_loss$log2OR

top_motifs <- rbind(top_motifs_gain, top_motifs_loss)
pp <- plotVolcano(df, top_motifs, oddsRatioCol, signifCol)

pdf(file.path(outputDir, paste0("lolaVolcanoPlot_",region,"_",comparison,".pdf")), width = 8, height = 6)
print(pp)
dev.off()

}


plotVolcano <- function(df, top_motifs, oddsRatioCol, signifCol) {
  # Ensure the condition column is a factor with the desired order
  df$condition <- factor(df$condition, levels = c("loss", "gain"))
  top_motifs$condition <- factor(top_motifs$condition, levels = c("loss", "gain"))

  # Split the data for loss and gain to apply different x-axis limits
  df_loss <- df %>%
    filter(condition == "loss") %>%
    mutate(log2OR = pmin(log2OR, 0)) # Cap at 0
  df_gain <- df %>%
    filter(condition == "gain") %>%
    mutate(log2OR = pmax(log2OR, 0)) # Start at 0
  df_combined <- bind_rows(df_loss, df_gain)

  top_motifs_loss <- top_motifs %>%
    filter(condition == "loss") %>%
    mutate(log2OR = pmin(log2OR, 0))
  top_motifs_gain <- top_motifs %>%
    filter(condition == "gain") %>%
    mutate(log2OR = pmax(log2OR, 0))
  top_motifs_combined <- bind_rows(top_motifs_loss, top_motifs_gain)

  pp <- ggplot(df_combined) +
    aes_string(oddsRatioCol, signifCol) +
    geom_point() +
    theme_classic() +
    geom_point(data = df_combined %>% filter(differential == "Non-differential"), color = "black") +
    ggrepel::geom_text_repel(aes(label = description),
      data = top_motifs_combined,
      size = 5,
      color = ifelse(top_motifs_combined$condition == "gain", "red", "blue"),
      min.segment.length = 0,
      seed = 42,
      box.padding = 1,
      max.overlaps = Inf,
      segment.size = 0.5
    ) +
    ylab("-log10(Q-Value)") +
    xlab("log2(Odds-Ratio)") +
    theme(
      axis.text = element_text(size = 14),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
      strip.text = element_text(size = 12, face = "bold"),
      strip.background = element_blank()
    ) +
    guides(alpha = "none", size = guide_legend(order = 1)) +
    facet_wrap(~condition, scales = "free_x") # Allow free x-axis scaling per facet

  return(pp)
}
lolaVolcanoPlotC19 <- function(cell, lolaDb, outputDir, reg, pValCut = 2, database = "TF_motif_clusters",
                               signifCol = "qValue", cnames = c("loss", "gain"), n = 3, colorpanel = c()) {
  if (endsWith(outputDir, ".rds")) {
    bed.files <- readRDS(outputDir)
    df <- bed.files$region
    bed <- names(df)
    df <- df[[reg]]$distal %>% #Change when promoter or other regions

      dplyr::filter(collection == database) %>%
      dplyr::filter(userSet %in% c("rankCut_1000_hyper", "rankCut_1000_hypo")) %>%
      dplyr::mutate(condition = ifelse(userSet == "rankCut_1000_hyper", "gain", "loss"))
    df <- muRtools:::lolaPrepareDataFrameForPlot(lolaDb, df,
      scoreCol = "log2OR", signifCol = signifCol,
      orderCol = "maxRnk", includedCollections = database,
      pvalCut = pValCut, maxTerms = Inf, perUserSet = FALSE,
      groupByCollection = TRUE, orderDecreasing = NULL
    )
    df$log2OR <- ifelse(df$condition == "loss", -df$log2OR, df$log2OR)
    df$name <- df$description
  } else {
    bed.files <- list.files(paste0(outputDir, cell, "/reports/differential_data"), pattern = "lolaRes", full.names = TRUE)
    bed.files <- grep(bed.files, pattern = region, value = TRUE)
    bed.files_l <- grep(pattern = paste0("_cutL2fcPadj05loss.rds"), bed.files, value = TRUE)[n]
    bed.files_g <- grep(pattern = paste0("_cutL2fcPadj05gain.rds"), bed.files, value = TRUE)[n]
    bed.files <- c(bed.files_l, bed.files_g)
    names(bed.files) <- cnames
    data <- lapply(names(bed.files), function(bed) {
      a <- readRDS(bed.files[[bed]]) %>%
        dplyr::filter(collection == database)
      a <- muRtools:::lolaPrepareDataFrameForPlot(lolaDb, a,
        scoreCol = "log2OR", signifCol = signifCol,
        orderCol = "maxRnk", includedCollections = database,
        pvalCut = pValCut, maxTerms = Inf, perUserSet = FALSE,
        groupByCollection = TRUE, orderDecreasing = NULL
      ) %>%
        dplyr::filter(log2OR >= 0)

      if (!is.null(a)) {
        a <- a %>%
          dplyr::mutate(condition = rep(bed, nrow(a)), log2OR = log2(oddsRatio))
        if (bed == "loss") {
          a <- a %>%
            dplyr::mutate(log2OR = -log2OR)
        }
      }
      a
    })
    df <- rlist::list.rbind(data)
    df$name <- gsub("\\s*\\[.*\\]", "", df$name)
  }
  if (length(bed.files) == 0) {
    stop("No files found with the specified pattern.")
  }
  if (signifCol == "qValue") {
    signifCol <- "qValueLog"
    df$qValueLog <- -log10(df$qValue)
    df$differential <- ifelse(df$qValueLog > 1.3, "Differential", "Non-differential")
  }


  top_motifs_gain <- df %>%
    dplyr::filter(differential == "Differential", condition == "gain") %>%
    arrange(desc(qValueLog)) %>%
    dplyr::filter(qValueLog > 2) %>%
    top_n(10, qValueLog)

  top_motifs_loss <- df %>%
    dplyr::filter(differential == "Differential", condition == "loss") %>%
    arrange(desc(qValueLog)) %>%
    dplyr::filter(qValueLog > 2) %>%
    top_n(10, qValueLog)

  top_motifs <- bind_rows(top_motifs_gain, top_motifs_loss)

  oddsRatioCol <- "log2OR"
  if (!is.element(oddsRatioCol, colnames(df))) {
    logger.error("Invalid LOLA result. Could not find a valid column containing odds ratios.")
  }
  pp <- plotVolcano(df, top_motifs, oddsRatioCol, signifCol)

  return(list(plot = pp, df = df))
}