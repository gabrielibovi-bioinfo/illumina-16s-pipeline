################################################################################
# Script: 04_picrust2_plots.R
# Description: Functional abundance barplot from PICRUSt2 output.
#              Visualizes predicted metabolic functions (KO level) by
#              functional group (AOB, NOB, ANAMMOX, GAO, PAO, DNB)
#              for the top 5 most abundant families.
#
# Steps performed:
#   1. Load phyloseq object and sample metadata
#   2. Define functional KO groups
#   3. Load PICRUSt2 metagenome contribution table
#   4. Merge with taxonomy and metadata; calculate normalized abundance
#   5. Filter top 5 families by total abundance per functional group
#   6. Generate barplot faceted by functional group
#
# Input files:
#   - inter_files/ps_silva.rds
#   - metadados.txt
#   - output_picrust/KO_metagenome_out/pred_metagenome_contrib.tsv.gz
#
# Output files:
#   - tabelas/functions_table.xlsx
#   - figuras/func_norm_family_5.tiff
################################################################################

library(tidyverse)
library(phyloseq)
library(readr)
library(ggsci)
library(writexl)
library(ggpubr)
library(scales)
library(dplyr)
library(ggplot2)
library(castor, quietly = TRUE)
library(viridis)
library(RColorBrewer)

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
# SECTION 2: Define functional KO groups
# ==============================================================================

AOB    <- c("K10944", "K10945", "K10946", "K10535", "K05601", "K15864")
# ammonia monooxygenase (AMO), hydroxylamine oxidoreductase (HAO)

NOB    <- c("K00370", "K00371")
# nitrite oxidoreductase (NOR/NXR)

ANAMMOX <- c("K20932", "K20933", "K20934", "K20935")
# hydrazine hydrolase, hydrazine dehydrogenase

GAO    <- c("K20812", "K00975", "K00688", "K02438")
# glycogen synthase (glgA), glucose-1-phosphate adenylyltransferase (glgC),
# glycogen phosphorylase (glgP), glycogen debranching enzyme (glgX)

PAO    <- c("K00937", "K22468")
# polyphosphate kinase (ppk, ppk2)

DNB    <- c("K00372", "K00360", "K00367", "K00370", "K00371", "K00373", "K00374", "K10534")
# assimilatory nitrate reductase (nasA, nasB) and related subunits

ko_ids  <- c(AOB, NOB, ANAMMOX, GAO, PAO, DNB)
groups  <- c(rep("AOB", length(AOB)), rep("NOB", length(NOB)),
             rep("ANAMMOX", length(ANAMMOX)), rep("GAO", length(GAO)),
             rep("PAO", length(PAO)), rep("DNB", length(DNB)))

KOs <- data.frame(group = groups, ko = ko_ids)

# ==============================================================================
# SECTION 3: Remove chloroplast and mitochondria ASVs
# ==============================================================================

exclude_asvs <- c(
  as.vector(na.omit(rownames(tax_table(seqtab)[tax_table(seqtab)[, "order"]  == "Chloroplast", ]))),
  as.vector(na.omit(rownames(tax_table(seqtab)[tax_table(seqtab)[, "family"] == "Mitochondria", ])))
)

taxons      <- data.frame(tax_table(seqtab))
taxons$ASV  <- rownames(taxons)

sample_info <- data.frame(sample_data(seqtab))

# ==============================================================================
# SECTION 4: Load PICRUSt2 output and build functional abundance table
# ==============================================================================

ko_table <- read.table(
  "output_picrust/KO_metagenome_out/pred_metagenome_contrib.tsv.gz",
  fill = TRUE, header = TRUE
)

func_table <- ko_table %>%
  merge(KOs, by.x = "function.", by.y = "ko") %>%
  filter(!taxon %in% exclude_asvs) %>%
  merge(taxons,      by.x = "taxon",  by.y = "ASV") %>%
  merge(sample_info, by.x = "sample", by.y = "sampleName") %>%
  group_by(sample, group, domain, phylum, class, order, family, genus) %>%
  summarize(abundance = sum(norm_taxon_function_contrib), .groups = "drop") %>%
  select(
    Domain   = domain,
    Phylum   = phylum,
    Class    = class,
    Order    = order,
    Family   = family,
    Genus    = genus,
    Sample   = sample,
    Group    = group,
    Abundance = abundance
  )

# write_xlsx(func_table, "tabelas/functions_table.xlsx")

# ==============================================================================
# SECTION 5: Filter top 5 families by total abundance
# ==============================================================================

top5_families <- func_table %>%
  filter(!is.na(Family) & !is.na(Abundance)) %>%
  group_by(Family) %>%
  summarise(Total = sum(Abundance), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  slice_head(n = 5) %>%
  pull(Family)

func_table_top5 <- func_table %>%
  filter(Family %in% top5_families) %>%
  group_by(Sample, Family, Group) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  arrange(desc(Abundance))

func_table_top5$Group <- factor(
  func_table_top5$Group,
  levels = c("AOB", "NOB", "DNB", "PAO", "GAO", "ANAMMOX")
)

# ==============================================================================
# SECTION 6: Barplot faceted by functional group
# ==============================================================================

color_palette <- colorRampPalette(brewer.pal(6, "PuOr"))(6)

func_barplot <- ggplot(func_table_top5,
                       aes(x = Sample, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", width = 0.9, position = position_dodge()) +
  facet_grid(cols = vars(Group), switch = "y", space = "free") +
  theme_classic(base_size = 14) +
  theme(
    text              = element_text(family = "serif", color = "black", face = "plain"),
    axis.title.x      = element_blank(),
    axis.text.x       = element_text(size = 14),
    axis.text.y       = element_blank(),
    axis.ticks.y      = element_blank(),
    panel.background  = element_rect(fill = "gray98"),
    panel.grid.major.y = element_line(colour = "white")
  ) +
  labs(
    y    = "Functional normalized abundance (adjusted count)",
    fill = "Family contribution"
  ) +
  scale_fill_manual(values = rev(color_palette)) +
  scale_x_discrete(
    labels = c("RF1" = "Inoculum", "R1" = "R1", "R2" = "R2", "R3" = "R3"),
    limits = c("RF1", "R1", "R2", "R3")
  )

print(func_barplot)
# ggsave("figuras/func_norm_family_5.tiff",
#        plot = func_barplot, width = 2900, height = 1400, units = "px", dpi = 310)
