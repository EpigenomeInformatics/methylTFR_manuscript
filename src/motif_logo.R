# Load necessary libraries
suppressPackageStartupMessages(library(JASPAR2020))
suppressPackageStartupMessages(library(TFBSTools))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(motifmatchr))
suppressPackageStartupMessages(library(ChrAccR))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(BSgenome))
suppressPackageStartupMessages(library(muLogR))

# Define the directory to save motif logos
motif_dir <- "/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/figures/motif_logos/"
if (!dir.exists(motif_dir)) {
  dir.create(motif_dir)
}
source("/scratch/icbb/igunduz/methylTFR_manuscript/github/methylTFR_manuscript/src/utils.R")

# Function to extract and save motif logos for specific transcription factors
save_motif_logos <- function(tfNames, motifDb = "jaspar2020", motif_dir) {
  
  # Fetch the motif data from the specified database
  motifObj <- prepareMotifmatchr("hg38", motifDb)$motifs
  
  for (tf in tfNames) {
    # Use a regex pattern to match the exact transcription factor name
    full_mn <- grep(tf, names(motifObj), value = TRUE, ignore.case = TRUE)
    
    if (length(full_mn) > 0) {
      # Open a PDF device
      pdf_file <- paste0(motif_dir, tf, "_motif_logo.pdf")
      pdf(pdf_file, width = 6, height = 2)
      
      for (mn in full_mn) {
        pwm <- motifObj[[mn]]  # Extract the PWM for the current match

        # Generate the motif logo using hmSeqLogo and save it
        grid.newpage()
        hmSeqLogo(pwm)
        
        # Add a title with the original motif name
        grid.text(mn, y = unit(1, "npc") - unit(1, "lines"), just = "center", gp = gpar(fontsize = 10))
      }

      # Close the PDF device
      dev.off()
      
      cat("Saved motif logos for", tf, "as", pdf_file, "\n")
    } else {
      cat("Motif for transcription factor", tf, "not found in the database.\n")
    }
  }
}

# Example usage
tfNames <-  c(
  "TFAP2B", "PAX5", "PAX9", "PAX1", "NHLH2", "ASCL1", "NHLH1", "BHLHE22", "FERD3L", "PAX6",
  "EMX1", "PAX4", "EN1", "LHX1", "ELF4", "ELF2", "ETV5", "ETV6", "ELF5", "SPIC",
  "SPIB", "SPI1", "EHF", "ELF3", "IKZF1", "TFAP2C", "TFAP2B", "TFAP2A", "TFAP2E", "TFAP2C",
  "TGIF1", "CREB3L4", "PBX3", "TAL1::TCF3", "MYOG", "ATOH1", "MYF5", "BHLHA15", "ZBTB18", "EBF3",
  "EBF1", "TFAP4", "NEUROD1", "VSX2", "EGR4", "DPRX", "SOX8", "CUX1", "CUX2", "TEAD3"
)
save_motif_logos(tfNames, motifDb = "jaspar2020", motif_dir = motif_dir)

