# Puritz et al. — CASE: Coastal Acidification and Sewage Effluent

Data and reproducible analysis code for a multiple-stressor selection experiment on larval eastern oyster (*Crassostrea virginica*). The CASE study asks how coastal acidification (CA) and sewage effluent (SE) — two stressors that co-occur in urbanized estuaries — affect larval survival and drive selection on standing genetic variation, both alone and in combination.

> **Status:** Manuscript in preparation. This repository is under active development and results should be considered preliminary until the paper is published. A companion RNAseq paper using the same experimental design is also in progress.

## Authors

- **Jonathan B. Puritz** — Department of Biological Sciences, University of Rhode Island
- **Johanna Harvey** — Department of Natural Resources Science, University of Rhode Island
- **Katie Lotterhos** — Marine Science Center, Northeastern University

## Background

Coastal marine organisms rarely face one stressor at a time. In urbanized estuaries, eutrophication couples coastal acidification (CA) to sewage effluent (SE): nutrient loading stimulates microbial CO₂ production, and even treated effluent runs lower in pH and higher in CO₂ than ambient seawater. CA and SE have been characterized separately but never together in an early life stage.

Wild-broodstock larvae were exposed for 24 hours to one of four treatments in a 2 × 2 factorial design crossing CA and SE:

| Treatment | Acidification | Sewage effluent |
|-----------|---------------|-----------------|
| **CON**   | ~400 µatm pCO₂ (control) | none |
| **CA**    | ~2,800 µatm pCO₂ | none |
| **SE**    | control | 5% v/v treated effluent |
| **CASE**  | ~2,800 µatm pCO₂ | 5% v/v treated effluent |

Larvae came from four independent spawns (SP1–SP4, wild broodstock from Ipswich and Barnstable, MA, 2016–2017), each an independent genetic background. SP1 was phenotyped for survival only; SP2, SP3, and SP4 were also sequenced. Pooled expressed exome capture sequencing (EecSeq) estimated allele frequencies before and after exposure in each sequenced spawn, and loci under selection were identified with an adapted Cochran–Mantel–Haenszel (CMH) test (ACER) within and across spawns, with control-significant loci removed to minimize experimental artifacts and family/maternal effects.

### Key findings

- **Mortality was additive, not synergistic.** Survival declined CON > CA > SE > CASE. CA alone never differed from control; CASE was statistically indistinguishable from SE alone (Tukey *p* ≈ 1.0). A factorial CA × SE GLMM found no significant interaction (*p* = 0.88) — effluent, not acidification, drove the mortality effect.
- **The genomic response was compositionally synergistic even though its magnitude stayed additive.** CASE recruited far more exclusive loci and genes than CA and SE combined (1,941 CASE-exclusive loci vs. 1,136 for CA + SE combined), while the aggregate size of allele-frequency change tracked the additive sum of the single-stressor effects (orthogonal-regression slope ≈ 0.99).
- **Outlier reproducibility across spawns rose with stressor complexity even as enrichment over baseline fell.** Loci replicating in ≥2 of 3 spawns increased from 222 (CA) to 647 (SE) to 1,197 (CASE), while fold-enrichment over a control-matched null fell from 403× to 235× to 116×.
- **Functional enrichment.** All three treatments were enriched for proteostasis (ubiquitin–proteasome pathways). CASE alone added ribosomal/translational-repression, mitochondrial calcium import, calcium-dependent proteolysis, ER-associated degradation, and a steroid-metabolizing (CYP17-like) module consistent with effluent-linked endocrine disruption.
- **Selection acted across a continuum of effect sizes.** Common, replicated (Core) loci showed small, consistent frequency shifts; rare (Private) loci carried larger, background-dependent shifts, with diversity (expected heterozygosity) rising at selected loci in a roughly 5–6-fold gradient from Core to Private.

## Repository structure

```
Puritz_etal_CASE/
├── analysis/                               # R Markdown reproducible analysis
│   ├── Final_reproducible_analysis.Rmd         # main pipeline (survival → outliers → genomic synergy → tiers)
│   ├── sorted.ref3.0.exon*.bed                 # exon annotation (capture target regions)
│   ├── GO/                                     # GO enrichment inputs/outputs, per treatment and tier
│   ├── figures/                                # generated figures (main/ and supp/)
│   └── results/                                # generated tables (additivity, interaction, He-null tests)
├── scripts/                                # preprocessing & pop-gen utilities (see below)
├── data/                                    # phenotype & water-chemistry inputs
│   ├── CASE_FINAL_Mortality.txt
│   ├── CASE_LARVAL_SIZE.txt
│   ├── CASE_LNDA.txt
│   ├── Water_Chemistry.csv                     # flat export of Water_Chemistry.xlsx (analysis input)
│   ├── Water_Chemistry*.xlsx                   # archival workbooks (CO2SYS output, summary)
│   └── Water_chemistry_raw/                    # raw VINDTA run files
├── CASE_environment.yaml                   # conda environment for the bioinformatics pipeline
└── random_draw_environment.yaml            # environment for T0 pseudo-replicate resampling
```


## Pipeline overview

The full path from raw reads to figures has two stages.

**1. Bioinformatic preprocessing** (bash/Python/Perl, `scripts/`). Read assembly and mapping with `dDocent` (`dDocent_ngs_3.1`), variant calling and filtering (`CASE_PVCF_filter2.sh`), conversion of the filtered VCF to PoPoolation2 synchronized format (`VCFtoPopPool.py`), coverage augmentation (`add_cov_sync`, `add_cov_sync_IS`), allele-frequency polarization (`polarize_freqs`), and pool subsampling (`sub_sample.py`, `subsample-synchronized.pl`). Per spawn, the four initial-timepoint (T0) libraries are pooled and resampled (`sum.sh`, `sub_sample.py`) into depth-matched pseudo-replicates matching the number of post-exposure replicates, so pre/post pools are directly comparable. CMH/frequency-difference estimation uses `snp-frequency-diff.pl`. These steps require external tools and large intermediates and are documented as `eval=FALSE` commands in the analysis notebook.

**2. Statistical analysis and figures** (R, `analysis/`). `Final_reproducible_analysis.Rmd` is the main entry point, using canonical spawn labels from `spawn_labels.R`. It loads the allele-frequency and CMH objects, applies FDR control (`qvalue`, π₀ = 1), assigns loci to confidence tiers, and builds the main figures and supplements. Four supplemental R Markdown documents cover the survival-additivity test, the He/matched-null diversity test, and additional GO-background and additivity checks. The pipeline is organized around the manuscript's results order (survival → outliers & repeatability → genomic synergy → tier architecture → diversity).

### Locus confidence tiers

Reproducible outliers are assigned hierarchically to one of three tiers, each locus appearing once:

- **Core** — q < 0.10 in all three sequenced spawns, or q < 0.05 in at least two; common alleles with small, consistent shifts reflecting parallel selection on standing variation. 5,312 loci total (662 CA, 1,661 SE, 2,989 CASE).
- **Convergent** — resolved only in the pooled across-spawn CMH (per-spawn q < 0.01 in ≥1 spawn and across-spawn q < 0.0001); individually weak, polygenic signal. 1,040 loci total (62 CA, 166 SE, 812 CASE).
- **Private** — top 1% of per-treatment q-values among control-insignificant loci, using spawn-specific thresholds; rarer alleles with larger, background-dependent shifts. 1,121 loci total (233 CA, 249 SE, 639 CASE).

Effective population size and per-spawn q thresholds are set per spawn to account for differences in breeder number and coverage (e.g., Ne = 250 for SP2 vs. 1,000 for SP3/SP4). See the analysis notebook and manuscript methods for exact values. *(Tiers were renamed from an earlier Aggregate/Block-Specific scheme; the names and thresholds above reflect the 2026-08-25 manuscript draft.)*

## Reproducing the analysis

**Bioinformatics environment** (conda/mamba):

```bash
conda env create -f CASE_environment.yaml
conda activate CASE
```

**R analysis.** Analyses in the manuscript were run under R 3.6.0. Render the main notebook from the repository root:

```r
rmarkdown::render("analysis/Final_reproducible_analysis.Rmd")
```

Required R packages: `ggplot2`, `tidyr`, `plyr`, `dplyr`, `qvalue`, `stringr`, `pcadapt`, `poolSeq`, `ACER`, `scales`, `data.table`, `patchwork`, `ggridges`, `ggrain`, `grid`, `ggman`, `VennDiagram`, `topGO`, `kableExtra`, `glmmTMB`, and `igraph`. `png` is also required (called via `png::readPNG` when assembling figure panels). Figures are written to `analysis/figures/main/` and `analysis/figures/supp/`; tables to `analysis/results/`.

> Raw sequencing reads are not yet deposited — the manuscript's data-availability statement currently has a placeholder for the NCBI SRA/BioProject accession. Large intermediate files are not tracked in this repository.

## Data and permits

Wild broodstock were collected under Massachusetts Division of Marine Fisheries scientific collection permit #173000.

## Citation

Working title (target journal: PLOS Genetics): *"Synergistic selection across the genome of eastern oyster larvae under coastal acidification and sewage effluent: a CASE study."* A full citation will be added once the manuscript is published.

## Contact

Jonathan B. Puritz — Department of Biological Sciences, University of Rhode Island — jpuritz@uri.edu

## Funding

NSF and NOAA/Sea Grant. *[Add specific award numbers before publication.]*
