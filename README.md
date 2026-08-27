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
| [`methylTFR`](https://github.com/EpigenomeInformatics/methylTFR) | analysis package |
| [`methylTFRAnnotationHg38`](https://github.com/EpigenomeInformatics/methylTFRAnnotationHg38) | motif binding sites, GC frequency tables and the genome-wide GC distribution for hg38 |
| [`methylTFRAnnotationMm10`](https://github.com/EpigenomeInformatics/methylTFRAnnotationMm10) | the same resources for mm10 |
| [`methylTFRAnnotationBuilder`](https://github.com/EpigenomeInformatics/methylTFRAnnotationBuilder) | builds annotation resources for a new genome or motif set |
| [`viper`](https://bioconductor.org/packages/viper) | transcription factor activity from expression, used by the Blueprint VIPER integration |
| [`dorothea`](https://bioconductor.org/packages/dorothea) | the regulon VIPER is run against |

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
BiocManager::install(c("viper", "dorothea"))  # only needed for 07b and 09
```

```r
devtools::install_github("EpigenomeInformatics/methylTFR")
devtools::install_github("EpigenomeInformatics/methylTFRAnnotationHg38")

# Only needed outside hg38
devtools::install_github("EpigenomeInformatics/methylTFRAnnotationMm10")
devtools::install_github("EpigenomeInformatics/methylTFRAnnotationBuilder")
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
to distal elements.

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
data/              bias-corrected deviations, one folder per dataset
tables/            sample annotation and analysis tables
figures/           figure panels, one folder per dataset and step
mtfr.yaml          conda environment
```

Every script opens with a settings block holding paths, motif set, thresholds
and palettes. Dated tags such as `RnBeads_230826` and `mTFR_devs_230826`
identify a processing run and appear in the output directory names.

Heatmaps, PCA plots, scatters and the MOFA panels are written into `figures/`
and travel with the repository. Footprints are not: there is one PDF per motif
per run, which is more than belongs in git, so they stay under `analysis.dir`
alongside the methylation data they were computed from.

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
| `07b_viper_prep.R` | VIPER transcription factor activity from those counts, against a DoRothEA regulon |
| `08_bp_mofa.R` | MOFA2 integration of deviations and expression, factor and modality figures, paired heatmaps |
| `09_bp_viper_mofa.R` | MOFA2 integration of deviations and VIPER activity, and the paired methylTFR / VIPER heatmap of the differential motifs |

`08` and `09` are two readings of the same samples: `08` pairs a motif with the
expression of its transcription factor, `09` with the activity VIPER infers for
it. `09` needs `06` for the differential motifs and `07b` for the activities.

A motif name is not a gene symbol. JASPAR writes heterodimers as `FOS::JUNB` and
variants as `JUN(var.2)`, so `09` splits on `::`, strips the variant suffix and
matches a motif when any of its subunits has an activity. Which subunit each row
was matched on is recorded in `tables/bp_motif_viper_pairs.csv`.


## CD4+ T memory cells

`src/memtcells/` — two donors, subtypes TN, TCM and TEM, run in order.

| Script | Description |
|---|---|
| `01_run_rnbeads_memTcells.R` | RnBeads import, preprocessing, differential methylation and LOLA motif enrichment |
| `02_run_mTFR_memTcells.R` | methylTFR deviations |
| `03_pca_memTcells.R` | PCA of the deviations by subtype and donor |
| `04_mFoots_memTcells.R` | observed minus expected footprints per subtype |
| `05_EO_mFoots_memTcells.R` | expected and observed footprint curves |
| `06_differential_TFs.R` | differential deviations against naive, paired activity and expression heatmap, expression against deviation scatters, LOLA enrichment and its comparison with methylTFR, differential methylation panels |


## ECHO

`src/echo/` — single cell WGBS with matched scATAC-seq.

| Script | Description |
|---|---|
| `00a_read_pseudobulk_devs.R` | reads the per cell type sample pseudobulk deviations and checks that they share one motif set |
| `00b_merge_pseudobulk_devs.R` | merges them into a single `methylTFRdeviations` object carrying the sample annotation |
| `01_mofa_integration.R` | MOFA2 on the chromVAR and methylTFR pseudobulk matrices |
| `02_echo_cell_pseudobulks.R` | cell type pseudobulks from the methylation data |
| `03_mfoot_meth.R` | observed minus expected methylation footprints |
| `03_mfoot_atac.R` | chromatin accessibility footprints for the same factors |
| `04_plot.R` | modality contribution panel, paired methylTFR and chromVAR heatmap, correlation by methyl-SELEX call, per factor activity scatters, and the assembled figure |

The pseudobulk deviations are one object per cell type, with donor samples in
the columns, scored against `jaspar2020_distal`. Cell type is the token before
the first underscore of the sample name, and the rest of the annotation travels
in the `colData` of the objects themselves.

Both footprint scripts use subtraction, so the methylation and accessibility
panels are on the same scale, and both carry the mean activity of each cell type
in the legend.

---

## Bias-corrected deviations

`data/` holds the deviations as `methylTFRdeviations` objects, one RDS per
dataset and motif set, under `data/<dataset>/<motifSet>_deviations.RDS`.

| Dataset | Motif set | File |
|---|---|---|
| CD4+ T memory cells | `jaspar2020` | [`jaspar2020_deviations.RDS`](data/memtcells/jaspar2020_deviations.RDS) |
| CD4+ T memory cells | `jaspar2020_distal` | [`jaspar2020_distal_deviations.RDS`](data/memtcells/jaspar2020_distal_deviations.RDS) |
| BLUEPRINT | `jaspar2020` | [`jaspar2020_deviations.RDS`](data/blueprint/jaspar2020_deviations.RDS) |
| BLUEPRINT | `jaspar2020_distal` | [`jaspar2020_distal_deviations.RDS`](data/blueprint/jaspar2020_distal_deviations.RDS) |


### Downloading

Each link above opens the file on GitHub, where **Download raw file** saves it.
The whole set comes with the repository:

```bash
git clone https://github.com/EpigenomeInformatics/methylTFR_manuscript.git
```

A single object can also be read straight into R without cloning:

```r
library(methylTFR)

base <- "https://raw.githubusercontent.com/EpigenomeInformatics/methylTFR_manuscript/HEAD/data"
dev <- readRDS(gzcon(url(file.path(base, "memtcells", "jaspar2020_distal_deviations.RDS"), "rb")))
```

### Using them

```r
dev <- readRDS("data/memtcells/jaspar2020_distal_deviations.RDS")

deviations(dev)        # bias-corrected deviations, motifs x samples
deviationZScores(dev)  # row-wise Z-scores
colData(dev)           # sample annotation
rownames(dev)          # motif identifiers
```

The objects extend `SummarizedExperiment`, so the usual accessors and
subsetting apply. Loading them is enough to reproduce every downstream figure
of a dataset, which is to say every script after `02`.

---

## Reproducing an analysis

1. Create the environment and install the packages listed above.
2. Set `analysis.dir`, `github.dir`, the run tags and `distal.file` in the
   settings block of the scripts of interest.
3. Run the scripts of that dataset in numerical order. Steps whose output
   already exists are skipped, so a rerun only fills in what is missing.

To go straight to the analysis, point `dev.tag` at a directory holding the
files from `data/` and start at `03`.

Figures land in `figures/<dataset>/`, apart from the footprints, which stay
under `analysis.dir`. Tables reported in the manuscript are written to
`tables/`.