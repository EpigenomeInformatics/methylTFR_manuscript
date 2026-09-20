#!/usr/bin/env Rscript

# 1. Install Bioconductor packages without updating existing dependencies
cat("Installing Bioconductor packages...\n")
BiocManager::install(
  c("AnnotationHubData", "MOFA2", "viper", "dorothea"),
  update = FALSE,
  ask = FALSE
)

# 2. Install ArchR from GitHub without upgrading any dependencies
cat("Installing ArchR...\n")
devtools::install_github(
  "GreenleafLab/ArchR",
  ref = "master",
  repos = BiocManager::repositories(),
  upgrade = "never"
)

# 3. Install methylTFR and its Hg38 annotation package from GitHub
cat("Installing methylTFR and annotations...\n")
devtools::install_github(
  c(
    "EpigenomeInformatics/methylTFR"
  ),
  repos = BiocManager::repositories(),
  upgrade = "never"
)

cat("All manual installations complete!\n")