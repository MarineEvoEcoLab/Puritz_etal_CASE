# Puritz et al. — CASE: Coastal Acidification and Sewage Effluent

Data and reproducible analysis code for a multiple-stressor selection experiment on larval eastern oyster (*Crassostrea virginica*). The CASE study asks how coastal acidification (CA) and sewage effluent (SE) — two stressors that co-occur in urbanized estuaries — affect larval survival and drive selection on standing genetic variation, both alone and in combination.

> **Status:** Manuscript in preparation. This repository is under active development and results should be considered preliminary until the paper is published. A companion RNAseq paper using the same experimental design is also in progress.

## Background

Coastal marine organisms rarely face one stressor at a time. In urbanized estuaries, eutrophication from sewage effluent drives coupled low-oxygen and acidified conditions that climate change is expected to intensify. *C. virginica* larvae are especially vulnerable during their short pelagic stage, yet the combined genomic response to acidification and sewage effluent had not been characterized.

Wild-broodstock larvae were exposed for 24 hours to one of four treatments in a 2 × 2 factorial design crossing CA and SE:

| Treatment | Acidification | Sewage effluent |
|-----------|---------------|-----------------|
| **CON**   | ~400 µatm pCO₂ (control) | none |
| **CA**    | ~2,800 µatm pCO₂ | none |
| **SE**    | control | 5% v/v treated effluent |
| **CASE**  | ~2,800 µatm pCO₂ | 5% v/v treated effluent |

Pooled expressed exome capture sequencing (EecSeq) estimated allele frequencies before and after exposure across three independent spawning blocks (Blocks 10, 11, 12; Block 6 was phenotyped but not sequenced). Loci under selection were identified with Cochran–Mantel–Haenszel (CMH) tests across replicate blocks, with control-treatment outliers removed to minimize experimental artifacts and maternal effects.

### Key findings

- **Synergistic mortality** under combined CASE conditions: CA alone caused no significant mortality, SE alone caused significant mortality, and CASE caused the highest.
- **Genomic differentiation** across all three stressor treatments relative to controls.
- **Functional enrichment** of CASE selection candidates in protein folding, ubiquitin–proteasome catabolism, and cilium assembly.

## Repository structure

```
Puritz_etal_CASE/
├── analysis/                       # R Markdown reproducible analysis
│   ├── Final_reproducible_analysis.Rmd   # main pipeline (survival → outliers → genomic synergy → tiers)
│   ├── Final_reproducible_analysis.md     # knitted output
│   └── sorted.ref3.0.exon*.bed            # exon annotation (capture target regions)
├── scripts/                        # preprocessing & pop-gen utilities (see below)
├── data/                           # phenotype & water-chemistry inputs
│   ├── CASE_FINAL_Mortality.txt
│   ├── CASE_LARVAL_SIZE.txt
│   ├── CASE_LNDA.txt
│   └── Water_Chemistry*.xlsx
├── GO/                             # gene-ontology inputs and enrichment outputs
├── figures/                        # generated figures (main/ and supp/)
├── results/                        # generated tables and enrichment outputs
├── CASE_methods_draft.md           # materials & methods working draft
└── CASE_environment.yaml           # conda environment for the bioinformatics pipeline
```

## Pipeline overview

The full path from raw reads to figures has two stages.

**1. Bioinformatic preprocessing** (bash/Python/Perl, `scripts/`). Read assembly and mapping with `dDocent` (`dDocent_ngs_3.1`), variant calling and filtering (`CASE_PVCF_filter2.sh`), conversion of the filtered VCF to PoPoolation2 synchronized format (`VCFtoPopPool.py`), coverage augmentation (`add_cov_sync`, `add_cov_sync_IS`), allele-frequency polarization (`polarize_freqs`), pool subsampling (`sub_sample.py`, `subsample-synchronized.pl`), and CMH/frequency-difference estimation (`snp-frequency-diff.pl`). GO theme clustering uses `go_cluster_themes.py`. These steps require external tools and large intermediates and are documented as `eval=FALSE` commands in the analysis notebook.

**2. Statistical analysis and figures** (R, `analysis/`). `Final_reproducible_analysis.Rmd` is the single entry point. It loads the allele-frequency and CMH objects, applies FDR control (`qvalue`), assigns loci to confidence tiers, and builds the main figures and supplements. The document is organized around the manuscript figures (survival → outliers & repeatability → genomic synergy → tier architecture).

### Locus confidence tiers

Selection candidates are grouped into three tiers of increasing stringency:

- **Core** — replicated across multiple blocks; highest confidence.
- **Aggregate** — significant per-block plus a stringent across-block threshold.
- **Block-Specific** — top per-block outliers not captured by the other tiers.

Per-block significance thresholds are corrected for known differences among blocks (e.g., drift inflation from fewer broodstock, lower coverage/power); multi-block replication in the Core tier controls for drift without correction. See the analysis notebook for exact thresholds.

## Reproducing the analysis

**Bioinformatics environment** (conda/mamba):

```bash
conda env create -f CASE_environment.yaml
conda activate CASE
```

**R analysis.** Render the main notebook from the repository root:

```r
rmarkdown::render("analysis/Final_reproducible_analysis.Rmd")
```

Required R packages: `ggplot2`, `tidyr`, `plyr`, `dplyr`, `qvalue`, `stringr`, `pcadapt`, `poolSeq`, `ACER`, `scales`, `data.table`, `patchwork`, `ggridges`, `ggrain`, `grid`, `ggman`, `VennDiagram`, `topGO`, `kableExtra`, `glmmTMB`, and `igraph`. `png` is also required (called via `png::readPNG` when assembling figure panels). Figures are written to `figures/main/` and `figures/supp/`; tables to `results/tables/`.

> Raw sequencing reads and large intermediate files are not tracked in this repository. Add a data-availability note (e.g., NCBI SRA/BioProject accession) here once deposited.

## Data and permits

Wild broodstock were collected under Massachusetts Division of Marine Fisheries scientific collection permit #173000.

## Citation

A citation will be added here once the manuscript is published.

## Contact

Jonathan B. Puritz — Department of Biological Sciences, University of Rhode Island — jpuritz@uri.edu

## Funding

NSF and NOAA/Sea Grant. *[Add specific award numbers before publication.]*
