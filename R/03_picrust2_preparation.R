################################################################################
# Script: 03_picrust2_preparation.R
# Description: Exports ASV sequences (FASTA) and abundance table (BIOM) from
#              phyloseq for use in PICRUSt2 functional prediction.
#
# Steps performed:
#   1. Load phyloseq object and sample metadata
#   2. Export ASV representative sequences as FASTA (ASVseqs.fasta)
#   3. Export ASV count table as BIOM (ASVcounts.biom)
#
# NOTE: PICRUSt2 must be run from the terminal after this script.
#       See the ### BASH ### section at the bottom for the exact commands.
#
# Input files:
#   - inter_files/ps_silva.rds
#   - metadados.txt
#
# Output files:
#   - ASVseqs.fasta
#   - ASVcounts.biom
#   - output_picrust/KO_metagenome_out/pred_metagenome_contrib.tsv.gz
################################################################################

library(tidyverse)
library(phyloseq)
library(readr)
library(biomformat)

options(scipen = 999)
set.seed(123)

# ==============================================================================
# SECTION 1: Load phyloseq object and sample metadata
# ==============================================================================

ps    <- readRDS("inter_files/ps_silva.rds")
samdf <- read.table("metadados.txt", header = TRUE, encoding = "UTF-8", sep = "\t")

order_samples <- paste0(samdf$sampleName)
seqtab        <- prune_samples(sample_names(ps) %in% order_samples, ps)
rownames(samdf) <- samdf$sampleName
sample_data(seqtab) <- as.data.frame(samdf)

# ==============================================================================
# SECTION 2: Export ASV representative sequences as FASTA
# ==============================================================================

fasta_file <- "ASVseqs.fasta"

if (file.exists(fasta_file)) file.remove(fasta_file)

seqs_fasta <- seqtab %>%
  refseq() %>%
  data.frame() %>%
  rownames_to_column() %>%
  dplyr::rename(ASV = rowname, FASTA = ".")

for (i in 1:nrow(seqs_fasta)) {
  cat(paste0(">", seqs_fasta$ASV[i]), file = fasta_file, sep = "\n", append = TRUE)
  cat(seqs_fasta$FASTA[i],            file = fasta_file, sep = "\n", append = TRUE)
}

cat("FASTA file written:", fasta_file, "\n")

# ==============================================================================
# SECTION 3: Export ASV count table as BIOM
# ==============================================================================

asv_table <- seqtab %>%
  otu_table(taxa_are_rows = FALSE) %>%
  data.frame()

biomformat::make_biom(t(asv_table)) %>%
  biomformat::write_biom(biom_file = "ASVcounts.biom")

cat("BIOM file written: ASVcounts.biom\n")

# ==============================================================================
# SECTION 4 ### BASH ### — Run PICRUSt2 in the terminal
# ==============================================================================
# After this script, activate the PICRUSt2 conda environment and run:
#
#   conda activate picrust2
#
#   # Option A: Full pipeline (recommended first run)
#   picrust2_pipeline.py \
#     -s ASVseqs.fasta \
#     -i ASVcounts.biom \
#     -o output_picrust \
#     -p 3 --stratified
#
#   # Option B: Step-by-step pipeline
#   place_seqs.py \
#     --study_fasta ASVseqs.fasta \
#     --ref_dir ~/anaconda3/envs/picrust2/lib/python3.9/site-packages/picrust2/default_files/prokaryotic/pro_ref \
#     --out_tree output_picrust/out.tre \
#     --processes 6 \
#     --intermediate output_picrust/intermediate/place_seqs \
#     --min_align 0.8 \
#     --chunk_size 1000 \
#     --placement_tool epa-ng \
#     --verbose
#
#   hsp.py \
#     --tree output_picrust/out.tre \
#     --output output_picrust/marker_predicted_and_nsti.tsv.gz \
#     --processes 6 \
#     --in_trait 16S
#
#   hsp.py \
#     --tree output_picrust/out.tre \
#     --output output_picrust/KO_predicted.tsv.gz \
#     --processes 6 \
#     --in_trait KO
#
#   metagenome_pipeline.py \
#     --input ASVcounts.biom \
#     --function output_picrust/KO_predicted.tsv.gz \
#     --marker output_picrust/marker_predicted_and_nsti.tsv.gz \
#     --min_reads 1 \
#     --min_samples 1 \
#     --out_dir output_picrust/KO_metagenome_out \
#     --max_nsti 2.0 \
#     --strat_out
#
#   conda deactivate
#
# Results are visualized in: 04_picrust2_plots.R
