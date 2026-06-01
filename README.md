# 🧬 illumina-16s-pipeline

> Reproducible workflow for Illumina 16S rRNA amplicon analysis — from raw paired-end reads to taxonomic classification, diversity analysis, and functional prediction.

![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20WSL2-blue)
![Conda](https://img.shields.io/badge/conda-required-green)
![R](https://img.shields.io/badge/R-%3E%3D4.1-276DC3?logo=r)
![Status](https://img.shields.io/badge/status-in%20development-yellow)

---

## Overview

This pipeline processes Illumina paired-end 16S amplicon sequencing data (e.g. V3–V4 region) using DADA2 for denoising and ASV calling, followed by taxonomic classification with MiDAS and SILVA databases, diversity analysis, and functional prediction with PICRUSt2.

```
Raw FASTQ (R1 + R2)
    │
    ├── [Step 1] Quality Assessment        → DADA2: plotQualityProfile
    │
    ├── [Step 2] Filter & Trim             → DADA2: filterAndTrim
    │                                        (maxEE = 3/3, truncLen = 277/217 bp)
    │
    ├── [Step 3] Error Learning & Denoising → DADA2: learnErrors, dada
    │
    ├── [Step 4] Merge Paired Reads        → DADA2: mergePairs
    │
    ├── [Step 5] Chimera Removal           → DADA2: removeBimeraDenovo
    │
    ├── [Step 6] Taxonomic Classification  → MiDAS 5.3 + SILVA 138.2 (IDTAXA)
    │
    ├── [Step 7] Diversity Analysis        → R: rarefaction, alpha diversity,
    │                                        taxa count, barplots, heatmaps
    │
    └── [Step 8] Functional Prediction     → PICRUSt2 + R: barplots by KO group
```

---

## Repository Structure

```
illumina-16s-pipeline/
├── R/
│   ├── 01_dada2_taxonomy.R          # DADA2 processing + MiDAS/SILVA taxonomy
│   ├── 02_diversity_abundance.R     # Rarefaction, alpha diversity, barplots, heatmaps
│   ├── 03_picrust2_preparation.R    # Export FASTA + BIOM for PICRUSt2
│   └── 04_picrust2_plots.R          # PICRUSt2 functional abundance barplots
├── data/                            # Place raw .fastq.gz files here (gitignored)
├── results/                         # Output files (gitignored)
├── README.md
├── .gitignore
└── LICENSE
```

---

## Requirements

- Ubuntu 20.04+ or WSL2
- R ≥ 4.1
- [Conda](https://docs.conda.io/en/latest/miniconda.html) (for PICRUSt2)

### R dependencies

```r
install.packages(c("tidyverse", "ggplot2", "reshape2", "readr", "readxl",
                   "writexl", "RColorBrewer", "vegan", "pheatmap",
                   "patchwork", "ggsci", "ggpubr", "scico", "scales",
                   "viridis", "devtools", "biomformat"))

if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("dada2", "phyloseq", "DECIPHER", "Biostrings",
                        "BiocStyle"))

# QsRutils (alpha diversity)
devtools::install_github("jfq3/QsRutils")
```

### PICRUSt2 (Conda)

```bash
conda create -n picrust2 -c bioconda -c conda-forge picrust2
```

### Databases

Download and update paths in `01_dada2_taxonomy.R` (`PATH_MIDAS` and `PATH_SILVA`):

| Database | Version | Download |
|----------|---------|----------|
| MiDAS    | 5.3     | [midasfieldguide.org](https://www.midasfieldguide.org/guide/downloads) |
| SILVA    | 138.2   | [SILVA SSU r138.2](https://www.arb-silva.de/download/archive/qiime) |

---

## Usage

### Step 1 — DADA2 processing and taxonomic classification

Place raw paired-end FASTQ files in `data/` (pattern: `*_R1_001.fastq.gz` / `*_R2_001.fastq.gz`).

Update the database paths at the top of the script:
```r
PATH_MIDAS <- "path/to/DADA2_taxonomy_MiDAS_5.3.fa"
PATH_SILVA  <- "path/to/SILVA_SSU_r138_2_2024.RData"
```

```r
source("R/01_dada2_taxonomy.R")
```

**Outputs:**
```
inter_files/out.rds
inter_files/errF.rds / errR.rds
inter_files/dadaFs.rds / dadaRs.rds
inter_files/mergers.rds / lengthsDist.rds
inter_files/seqtabNoC.rds
inter_files/track_df.Rda
inter_files/taxid_midas.rds / ps_midas.rds
inter_files/taxid_silva.rds / ps_silva.rds
```

### Step 2 — Diversity analysis and relative abundance

```r
source("R/02_diversity_abundance.R")
```

**Outputs:**
```
tabelas/ASV_table_bruto.xlsx
tabelas/indices_alpha_diversity.xlsx
figuras/rarecurve_total.tiff
figuras/taxa_numbers.tiff
figuras/barra_<rank>.tiff
figuras/heat_<rank>.tiff
```

### Step 3 — Export data for PICRUSt2

```r
source("R/03_picrust2_preparation.R")
```

Then run PICRUSt2 in the terminal (commands are documented at the bottom of the script):

```bash
conda activate picrust2
picrust2_pipeline.py -s ASVseqs.fasta -i ASVcounts.biom -o output_picrust -p 3 --stratified
conda deactivate
```

**Outputs:**
```
ASVseqs.fasta
ASVcounts.biom
output_picrust/KO_metagenome_out/pred_metagenome_contrib.tsv.gz
```

### Step 4 — Functional prediction plots

```r
source("R/04_picrust2_plots.R")
```

**Outputs:**
```
tabelas/functions_table.xlsx
figuras/func_norm_family_5.tiff
```

---

## Filtering Parameters

| Parameter    | Value        | Step                        |
|--------------|--------------|-----------------------------|
| Max N        | 0            | filterAndTrim               |
| Max EE (F/R) | 3 / 3        | filterAndTrim               |
| Trim left    | 14 bp (F)    | filterAndTrim (primer removal) |
| Truncate     | 277 / 217 bp | filterAndTrim               |
| Min overlap  | 12 bp        | mergePairs                  |
| Max NSTI     | 2.0          | PICRUSt2 metagenome_pipeline |

> Parameters are set for V3–V4 16S amplicons (~460 bp merged). Adjust `truncLen` and `trimLeft` based on your quality profiles and primer lengths.

---

## Functional Groups (PICRUSt2)

| Group   | Function                          | Key KOs |
|---------|-----------------------------------|---------|
| AOB     | Ammonia-oxidizing bacteria        | K10944, K10945, K10946, K10535, K05601, K15864 |
| NOB     | Nitrite-oxidizing bacteria        | K00370, K00371 |
| ANAMMOX | Anaerobic ammonium oxidation      | K20932, K20933, K20934, K20935 |
| GAO     | Glycogen-accumulating organisms   | K20812, K00975, K00688, K02438 |
| PAO     | Polyphosphate-accumulating org.   | K00937, K22468 |
| DNB     | Denitrifying bacteria             | K00372, K00360, K00367, K00370, K00371, K00373, K00374, K10534 |

---

## Tools and References

| Tool     | Version | Reference |
|----------|---------|-----------|
| DADA2    | ≥1.28   | [Callahan et al., 2016](https://doi.org/10.1038/nmeth.3869) |
| phyloseq | ≥1.38   | [McMurdie & Holmes, 2013](https://doi.org/10.1371/journal.ppat.1003531) |
| DECIPHER | ≥2.20   | [Wright, 2016](https://doi.org/10.1016/j.jmb.2016.04.020) |
| MiDAS   | 5.3     | [Karst et al., 2021](https://doi.org/10.1038/s41587-021-01095-3) |
| SILVA    | 138.2   | [Quast et al., 2013](https://doi.org/10.1093/nar/gks1219) |
| PICRUSt2 | ≥2.5   | [Douglas et al., 2020](https://doi.org/10.1038/s41587-020-0548-6) |
| vegan    | ≥2.6    | [Oksanen et al., 2022](https://CRAN.R-project.org/package=vegan) |

---

## Notes

- All scripts use relative paths from the project root. Run them from the `illumina-16s-pipeline/` directory.
- Database paths (`PATH_MIDAS`, `PATH_SILVA`) must be updated in `01_dada2_taxonomy.R` before running.
- Filtering parameters (`truncLen`, `trimLeft`, `maxEE`) should be adjusted based on your quality profiles and target amplicon region.
- PICRUSt2 bash commands are embedded as comments in `03_picrust2_preparation.R` — run them in the terminal after the R script.
- Raw data and intermediate files are gitignored. Use `.gitkeep` to preserve empty folders.
- This pipeline was developed and tested on WSL2 (Ubuntu 22.04).

---

## Author

**Gabriel Ibovi**
Bioinformatics | Metagenomics | Illumina Amplicon Sequencing
🔗 [github.com/gabrielibovi-bioinfo](https://github.com/gabrielibovi-bioinfo)

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.
