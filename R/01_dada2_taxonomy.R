################################################################################
# Script: 01_dada2_taxonomy.R
# Description: DADA2 amplicon processing pipeline for Illumina 16S rRNA data.
#              Performs quality filtering, error learning, denoising, chimera
#              removal, and taxonomic classification with MiDAS and SILVA.
#
# Steps performed:
#   1. Define directories and list raw FASTQ files
#   2. Quality assessment before filtering (QualityProfile plots)
#   3. Filter and trim reads
#   4. Learn error rates and run DADA2 denoising
#   5. Merge paired-end reads and build sequence table
#   6. Remove chimeras
#   7. Track reads through the pipeline
#   8. Assign taxonomy with MiDAS database
#   9. Assign taxonomy with SILVA database (IDTAXA method)
#  10. Build phyloseq objects and export
#
# Input files:
#   - data/*.fastq.gz           (raw paired-end FASTQ, R1/R2)
#   - metadados.txt             (sample metadata)
#   - <path>/DADA2_taxonomy.fa  (MiDAS 5.3 database)
#   - <path>/SILVA_SSU_r138_2_2024.RData (SILVA 138.2 database)
#
# Output files:
#   - inter_files/out.rds
#   - inter_files/errF.rds / errR.rds
#   - inter_files/dadaFs.rds / dadaRs.rds
#   - inter_files/mergers.rds / lengthsDist.rds
#   - inter_files/seqtabNoC.rds
#   - inter_files/track_df.Rda
#   - inter_files/taxid_midas.rds / ps_midas.rds
#   - inter_files/taxid_silva.rds / ps_silva.rds
################################################################################

library(tidyverse)
library(knitr)
library(BiocStyle)
library(ggplot2)
library(gridExtra)
library(phyloseq)
library(DECIPHER)
library(dada2)
library(stats)
library(parallel)
library(Biostrings)

options(scipen = 999)
set.seed(100)

# ==============================================================================
# SECTION 0: Database paths — update before running
# ==============================================================================

# Update these paths to match your local database locations
PATH_MIDAS <- "path/to/DADA2_taxonomy_MiDAS_5.3.fa"
PATH_SILVA  <- "path/to/SILVA_SSU_r138_2_2024.RData"

# ==============================================================================
# SECTION 1: Define directories and list raw FASTQ files
# ==============================================================================

fastq_folder       <- "data"
filter_folder      <- "filtered_reads"
intermediate_folder <- "inter_files"

filt_path         <- file.path(filter_folder)
intermediate_path <- file.path(intermediate_folder)

if (!file_test("-d", filt_path))         dir.create(filt_path)
if (!file_test("-d", intermediate_path)) dir.create(intermediate_path)

# List forward (R1) and reverse (R2) FASTQ files
fnFs <- sort(list.files(fastq_folder, pattern = "..._R1_001.fastq.gz"))
fnRs <- sort(list.files(fastq_folder, pattern = "..._R2_001.fastq.gz"))

sampleNames <- sapply(strsplit(fnFs, "-"), `[`, 1)

fnFs <- file.path(fastq_folder, fnFs)
fnRs <- file.path(fastq_folder, fnRs)

# ==============================================================================
# SECTION 2: Quality assessment — before filtering
# ==============================================================================

# Plot quality profiles for the first 5 samples (forward and reverse)
# Recommended figure size: 700 x 500 px
plotQualityProfile(fnFs[1:5])
plotQualityProfile(fnRs[1:5])

# ==============================================================================
# SECTION 3: Filter and trim reads
# ==============================================================================

# NOTE: truncLen values (277, 217) are set for V3-V4 amplicons (~460 bp merged).
# Adjust trimLeft and truncLen based on your primer lengths and quality profiles.
# Rule: forward + reverse - amplicon length >= 10 bp (minimum overlap for merging)

filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))
names(filtFs) <- sampleNames
names(filtRs) <- sampleNames

print("Filtering and trimming reads...")

out <- filterAndTrim(
  fnFs, filtFs, fnRs, filtRs,
  maxN      = 0,           # Remove reads containing Ns
  maxEE     = c(3, 3),     # Maximum expected errors per read (F, R)
  truncQ    = 2,           # Truncate at first base with quality < 2
  trimLeft  = c(14, 0),    # Trim primers from left: 14 bp (F), 0 bp (R)
  truncLen  = c(277, 217), # Truncate reads to fixed length (F, R)
  rm.phix   = TRUE,        # Remove PhiX spike-in contamination
  compress  = TRUE,
  verbose   = TRUE
)

saveRDS(out, paste0(intermediate_folder, "/out.rds"))

# Quality assessment — after filtering
plotQualityProfile(filtFs[1:5])
plotQualityProfile(filtRs[1:5])

# ==============================================================================
# SECTION 4: Learn error rates
# ==============================================================================

errF <- learnErrors(filtFs, multithread = TRUE)
saveRDS(errF, paste0(intermediate_folder, "/errF.rds"))

errR <- learnErrors(filtRs, multithread = TRUE)
saveRDS(errR, paste0(intermediate_folder, "/errR.rds"))

# ==============================================================================
# SECTION 5: Denoise reads with DADA2
# ==============================================================================

dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
saveRDS(dadaFs, paste0(intermediate_folder, "/dadaFs.rds"))

dadaRs <- dada(filtRs, err = errR, multithread = TRUE)
saveRDS(dadaRs, paste0(intermediate_folder, "/dadaRs.rds"))

gc()

# ==============================================================================
# SECTION 6: Merge paired-end reads
# ==============================================================================

mergers    <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, minOverlap = 12)
seqtabAll  <- makeSequenceTable(mergers)
rownames(seqtabAll) <- sampleNames
lengthsDist <- table(nchar(getSequences(seqtabAll)))

saveRDS(mergers,     paste0(intermediate_folder, "/mergers.rds"))
saveRDS(lengthsDist, paste0(intermediate_folder, "/lengthsDist.rds"))

gc()

# ==============================================================================
# SECTION 7: Remove chimeras
# ==============================================================================

seqtabNoC <- removeBimeraDenovo(seqtabAll,
                                method      = "consensus",
                                verbose     = TRUE,
                                multithread = 8)

saveRDS(seqtabNoC, paste0(intermediate_folder, "/seqtabNoC.rds"))

# ==============================================================================
# SECTION 8: Track reads through the pipeline
# ==============================================================================

getN <- function(x) sum(getUniques(x))

track <- cbind(
  out,
  sapply(dadaFs,  getN),
  sapply(dadaRs,  getN),
  sapply(mergers, getN),
  rowSums(seqtabNoC)
)

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sampleNames
track_df <- as.data.frame(track)
track_df$nonchim_percent <- track_df$nonchim * 100 / track_df$input

print(track_df)
save(track_df, file = paste0(intermediate_folder, "/track_df.Rda"))

gc()

# ==============================================================================
# SECTION 9: Taxonomic assignment — MiDAS 5.3
# ==============================================================================

taxid_midas <- assignTaxonomy(seqtabNoC, PATH_MIDAS, multithread = TRUE)
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
colnames(taxid_midas) <- ranks

saveRDS(taxid_midas, paste0(intermediate_folder, "/taxid_midas.rds"))

# Build phyloseq object
ps_midas <- phyloseq(
  otu_table(seqtabNoC, taxa_are_rows = FALSE),
  tax_table(taxid_midas)
)

dna_midas <- Biostrings::DNAStringSet(taxa_names(ps_midas))
names(dna_midas) <- taxa_names(ps_midas)
ps_midas <- merge_phyloseq(ps_midas, dna_midas)
taxa_names(ps_midas) <- paste0("ASV", seq(ntaxa(ps_midas)))

saveRDS(ps_midas, paste0(intermediate_folder, "/ps_midas.rds"))
print(paste0("Done! Phyloseq object saved to >> ", intermediate_folder, "/ps_midas.rds <<"))

# ==============================================================================
# SECTION 10: Taxonomic assignment — SILVA 138.2 (IDTAXA method)
# ==============================================================================

dna <- DNAStringSet(getSequences(seqtabNoC))
load(PATH_SILVA)

ids <- IdTaxa(dna, trainingSet, strand = "both", processors = 8, verbose = TRUE)

taxid_silva <- t(sapply(ids, function(x) {
  m    <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid_silva) <- ranks
rownames(taxid_silva) <- getSequences(seqtabNoC)

saveRDS(taxid_silva, paste0(intermediate_folder, "/taxid_silva.rds"))

# Build phyloseq object
ps_silva <- phyloseq(
  otu_table(seqtabNoC, taxa_are_rows = FALSE),
  tax_table(taxid_silva)
)

dna_silva <- Biostrings::DNAStringSet(taxa_names(ps_silva))
names(dna_silva) <- taxa_names(ps_silva)
ps_silva <- merge_phyloseq(ps_silva, dna_silva)
taxa_names(ps_silva) <- paste0("ASV", seq(ntaxa(ps_silva)))

saveRDS(ps_silva, paste0(intermediate_folder, "/ps_silva.rds"))
print(paste0("Done! Phyloseq object saved to >> ", intermediate_folder, "/ps_silva.rds <<"))
