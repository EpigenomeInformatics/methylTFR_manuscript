# Computational quantification of transcription factor activity from DNA methylation

Code repository for the manuscript *"methylTFR: Computational quantification of
transcription factor activity from DNA methylation"*.

methylTFR estimates transcription factor activity from bisulfite sequencing by
comparing the observed methylation around motif occurrences against the
methylation expected from GC content alone, and reporting the bias-corrected
deviation per motif and sample.

---

## Software

| Package | Purpose |
|---|---|
| `methylTFR` | analysis package |
| `methylTFRAnnotationHg38` | motif binding sites, GC frequency tables and the genome-wide GC distribution for hg38 |
| `methylTFRAnnotationMm10` | the same resources for mm10 |
| `methylTFRAnnotationBuilder` | builds annotation resources for a new genome or motif set |

### Installation

```bash
conda env create -f mtfr.yaml
conda activate mtfr
```

Two packages are installed separately because they conflict with the conda
solver:

```r
BiocManager::install("RnBeads")
BiocManager::install("AnnotationHubData")   # only needed to build annotations
```

```r
devtools::install_github("EpigenomeInformatics/methylTFR")
devtools::install_github("EpigenomeInformatics/methylTFRAnnotationHg38")
```

Annotation resources are retrieved from AnnotationHub on first use and cached
locally. The available motif sets are `altius`, `cisbpv2`, `jaspar2020` and
`jaspar2020_distal`.

### Motif sets

Sets ending in `_distal` provide GC frequency tables computed over distal
regulatory regions only and share unrestricted binding sites with their base
set. Pass the base set name to `getTFbindsites()` and the `_distal` name to
`getGCfreq()`:

```r
tf_bindsites <- getTFbindsites(motifSet = "jaspar2020")
gcfreqs      <- getGCfreq(motifSet = "jaspar2020_distal")
gc_dist      <- getGenomeGC()
```

The distal restriction uses Ensembl Regulatory Build v104, filtered to hg38 and
to distal elements, and is read from
`methylTFRAnnotationHg38_old/inst/extdata/distal_regions.RDS`. This is the
region set the `jaspar2020_distal` GC frequency tables were built against, and
the analyses require the same one.

---

## Repository layout

```
src/
  utils.R          shared plotting helpers
  lola_utils.R     LOLA enrichment plots
  motif_logo.R     motif logo helper
  blueprint/       BLUEPRINT bulk WGBS
  memtcells/       CD4+ T memory cell subtypes
  echo/            ECHO single cell WGBS with matched scATAC-seq
tables/            sample annotation and analysis tables
figures/           figure panels
mtfr.yaml          conda environment
```

Every script opens with a settings block holding paths, motif set, thresholds
and palettes. Dated tags such as `RnBeads_230826` and `mTFR_devs_230826`
identify a processing run and appear in the output directory names.

---

## BLUEPRINT

`src/blueprint/` — 147 healthy samples, run in order.

| Script | Description |
|---|---|
| `01_run_rneads_blueprint.R` | RnBeads set filtered to the healthy samples |
| `02_run_mTFR_blueprint.R` | methylTFR deviations for both motif sets |
| `03_pca_stability_blueprint.R` | PCA of the deviations, tiling1kb and distal methylation, and a random forest stability analysis on the resulting principal components |
| `04_mFoots_allcelltypes.R` | observed minus expected methylation footprints across five cell type groups |
| `05_EO_mFoots_allcelltypes.R` | expected and observed footprint curves |
| `06_differential_heatmap.R` | differential deviations across cell types, Z-score heatmap and cohort composition |
| `07_rna_mofa_prep.R` | RNA counts matrix and the WGBS to RNA sample map |
| `08_bp_mofa.R` | MOFA2 integration of deviations and expression, factor and modality figures, paired heatmaps |

Cancer samples are identified through the `DISEASE` column of the sample
annotation, which leaves 147 of 157 samples. Cell type groups are recorded in
`cellTypeGroup` and mapped to display names within each script.

## CD4+ T memory cells

`src/memtcells/` — two donors, subtypes TN, TCM and TEM, run in order.

| Script | Description |
|---|---|
| `01_run_rnbeads_memTcells.R` | RnBeads import, preprocessing, differential methylation and LOLA motif enrichment |
| `02_run_mTFR_memTcells.R` | methylTFR deviations |
| `03_pca_memTcells.R` | PCA of the deviations by subtype and donor |
| `04_mFoots_memTcells.R` | observed minus expected footprints per subtype |
| `05_EO_mFoots_memTcells.R` | expected and observed footprint curves |
| `06_differential_TFs.R` | differential deviations against naive, paired activity and expression heatmap, LOLA enrichment and its comparison with methylTFR, differential methylation panels |

TEMRA is a single unreplicated sample and is excluded from the group
comparisons.

With two donors per group the differential test is limited: for TEM against TN
145 motifs reach an adjusted p-value below 0.05, while for TCM against TN the
smallest attainable adjusted p-value is 0.095. `06_differential_TFs.R` reports
the smallest raw and adjusted p-value of each contrast, and exposes the
significance column, an effect size floor and a paired test as settings.

Differential regions in the methylation panels are defined by the RnBeads
combined rank, chosen per comparison with `auto.select.rank.cut()`. This matches
the `rankCut` user sets of the LOLA enrichment, so both describe the same
regions.

## ECHO

`src/echo/` — single cell WGBS with matched scATAC-seq.

| Script | Description |
|---|---|
| `01_mofa_integration.R` | MOFA2 on the chromVAR and methylTFR pseudobulk matrices |
| `02_echo_cell_pseudobulks.R` | cell type pseudobulks from the methylation data |
| `03_mfoot_meth.R` | observed minus expected methylation footprints |
| `03_mfoot_atac.R` | chromatin accessibility footprints for the same factors |
| `04_mtfr_scECHO.R` | methylTFR on the single cell allc files of the samples overlapping the ATAC data |
| `05_plot.R` | paired methylTFR and chromVAR heatmap, correlation by methyl-SELEX call, per factor activity scatters |
| `debug_sc_mtfr.R` | diagnostic that runs the pipeline serially on a single allc file |

Both footprint scripts use subtraction, so the methylation and accessibility
panels are on the same scale. The transcription factors shown are taken from the
ECHO MOFA factor table.

Single cell coverage is sparse relative to the distal regions, and a cell that
does not cover every GC bin inside them cannot be scored. `04_mtfr_scECHO.R`
therefore runs a coverage check first, writes the per cell result to
`cell_gcbin_qc.tsv` and excludes the cells that fall short.

---

## Bias-corrected deviations

The deviations are provided as `methylTFRdeviations` objects in RDS files.

- [The BLUEPRINT project](#)
- [CD4+ T memory cells](#)

```r
library(methylTFR)
dev <- readRDS("jaspar2020_distal_deviations.RDS")
deviations(dev)        # bias-corrected deviations
deviationZScores(dev)  # row-wise Z-scores
```

---

## Reproducing an analysis

1. Create the environment and install the packages listed above.
2. Set `analysis.dir`, the run tags and `distal.file` in the settings block of
   the scripts of interest.
3. Run the scripts of that dataset in numerical order. Steps whose output
   already exists are skipped, so a rerun only fills in what is missing.

Figures are written alongside the analysis outputs under `analysis.dir`. Tables
reported in the manuscript are written to `tables/`.
