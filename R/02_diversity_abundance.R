################################################################################
# Script: 02_diversity_abundance.R
# Description: Diversity analysis and relative abundance visualization from
#              DADA2/phyloseq output (Illumina 16S rRNA data).
#
# Steps performed:
#   1. Load phyloseq object and sample metadata
#   2. Export ASV abundance table (relative proportions)
#   3. Rarefaction curves
#   4. Alpha diversity indices (Chao1, Shannon, Simpson, etc.)
#   5. Unique taxa count plot by taxonomic rank
#   6. Relative abundance barplot by taxonomic rank
#   7. Heatmap of relative abundance by taxonomic rank
#
# Input files:
#   - inter_files/ps_silva.rds
#   - metadados.txt
#   - tables/ASV_table.xlsx  (from step 2 export)
#
# Output files:
#   - tables/ASV_table.xlsx
#   - tables/indices_alpha_diversity.xlsx
#   - figures/rarecurve.tiff
#   - figures/taxa_numbers.tiff
#   - figures/barplot_<rank>.tiff
#   - figures/heatmap_<rank>.tiff
################################################################################

library(phyloseq)
library(dada2)
library(dplyr)
library(tidyverse)
library(reshape2)
library(readr)
library(readxl)
library(writexl)
library(ggplot2)
library(RColorBrewer)
library(vegan)
library(QsRutils)
library(pheatmap)
library(patchwork)
library(ggsci)
library(ggpubr)
library(scico)

options(scipen = 999)
set.seed(123)

# ==============================================================================
# SECTION 1: Load phyloseq object and sample metadata
# ==============================================================================

ps     <- readRDS("inter_files/ps_silva.rds")
samdf  <- read.table("metadados.txt", header = TRUE, encoding = "UTF-8", sep = "\t")
order_samples <- paste0(samdf$sampleName)

seqtab <- prune_samples(sample_names(ps) %in% order_samples, ps)
rownames(samdf) <- samdf$sampleName
sample_data(seqtab) <- as.data.frame(samdf)

# Optional: subset by domain
# seqtab <- subset_taxa(seqtab, domain == "Bacteria")
# seqtab <- subset_taxa(seqtab, domain == "Archaea")

# ==============================================================================
# SECTION 2: Export ASV abundance table (relative proportions)
# ==============================================================================

asv_proportions_tab <- apply(otu_table(seqtab), 1, function(x) x / sum(x) * 100)
ASVs <- as.data.frame(asv_proportions_tab)
taxa <- as.data.frame(tax_table(seqtab))

ASVs$ASV <- rownames(ASVs)
tidy      <- ASVs %>% melt()
taxa$ASV  <- rownames(taxa)

merged <- merge(taxa, tidy, by = "ASV")
colnames(merged) <- c("ASV", "Domain", "Phylum", "Class", "Order",
                      "Family", "Genus", "Species", "Sample",
                      "Relative abundance (%)")

# write_xlsx(merged, "tables/ASV_table_bruto.xlsx")

# ==============================================================================
# SECTION 3: Rarefaction curves
# ==============================================================================

# tiff("figures/rarecurve.tiff", width = 3400, height = 2000, res = 300)
sample_names_vec <- sample_data(seqtab)[[2]]
par(cex.axis = 1.6, cex.lab = 1.6)

rarecurve_data <- rarecurve(
  otu_table(seqtab) %>% data.frame(),
  step  = 500, cex = 1,
  xlab  = "Sample size", ylab = "Species",
  label = FALSE, bty = "L", family = "serif"
)

for (i in seq_along(rarecurve_data)) {
  x_values <- attr(rarecurve_data[[i]], "Subsample")
  y_values <- rarecurve_data[[i]]
  label_x  <- x_values[length(x_values)]
  label_y  <- y_values[length(y_values)]

  name_length <- nchar(sample_names_vec[i])
  offset      <- name_length * 1800

  rect(
    xleft = label_x - offset, ybottom = label_y - 25,
    xright = label_x - 0,    ytop    = label_y + 25,
    border = "black", col = "white"
  )
  text(
    x = label_x, y = label_y,
    labels = sample_names_vec[i],
    pos = 2, cex = 1.1, col = "black", family = "serif"
  )
}
# dev.off()

# ==============================================================================
# SECTION 4: Alpha diversity indices
# ==============================================================================

a <- goods(otu_table(seqtab)) %>% rownames_to_column(var = "Samples")
b <- estimate_richness(seqtab,
       measures = c("Observed", "Chao1", "ACE", "Shannon",
                    "Simpson", "InvSimpson", "Fisher")) %>%
     select(-c("se.chao1", "se.ACE"))

indices       <- bind_cols(a, b)
indices_table <- indices[, c("Samples", "Observed", "Chao1", "ACE",
                             "Shannon", "Simpson", "InvSimpson",
                             "Fisher", "no.seqs")]

print(indices_table)

# write_xlsx(indices_table, "tables/indices_alpha_diversity.xlsx")

# Diversidade beta
ordination <- ordinate(seqtab, method="NMDS", distance="bray")
plot_ordination(seqtab, ordination, color="sampleName")

ordination <- ordinate(seqtab, method="PCoA", distance="bray")
plot_beta <- plot_ordination(seqtab, ordination, color="sampleName")
plot_beta

# ==============================================================================
# SECTION 5: Unique taxa count by taxonomic rank and sample
# ==============================================================================

merged_plot <- read_excel("tables/ASV_table.xlsx")

tax_levels <- c("Phylum", "Class", "Order", "Family")

filtered_data <- merged_plot %>%
  filter(`Relative abundance (%)` > 0)

taxa_counts_df <- filtered_data %>%
  pivot_longer(cols = all_of(tax_levels), names_to = "Level", values_to = "Taxon") %>%
  group_by(Level, Sample) %>%
  summarise(Taxa_Count = n_distinct(Taxon, na.rm = TRUE), .groups = "drop") %>%
  mutate(Sample_Label = case_when(
    Sample == "RF1" ~ "Inoculum",
    Sample == "R1"  ~ "R1",
    Sample == "R2"  ~ "R2",
    Sample == "R3"  ~ "R3",
    TRUE            ~ Sample
  ))

taxa_counts_df$Level <- factor(taxa_counts_df$Level,
                               levels = c("Phylum", "Class", "Order", "Family"))

taxa_counts_df$Sample_Label <- factor(taxa_counts_df$Sample_Label,
                                      levels = c("Inoculum", "R1", "R2", "R3"))

taxa_numbers <- ggplot(taxa_counts_df,
                       aes(x = Level, y = Taxa_Count, color = Sample_Label)) +
  geom_point(position = position_dodge(width = 0.6), size = 5) +
  scale_y_continuous(limits = c(0, max(taxa_counts_df$Taxa_Count) + 5)) +
  theme_bw(base_family = "serif", base_size = 13) +
  labs(x = "", y = "Number of taxa") +
  theme(
    panel.grid      = element_blank(),
    legend.position = c(0.15, 0.80),
    legend.title    = element_blank(),
    axis.text.x     = element_text(face = "bold", size = 13),
    axis.text.y     = element_text(size = 13)
  ) +
  scale_x_discrete(limits = c("Phylum", "Class", "Order", "Family")) +
  scale_color_manual(values = brewer.pal(n = 4, "PuOr"))

print(taxa_numbers)
# ggsave("figures/taxa_numbers.tiff",plot = taxa_numbers, height = 1200, width = 1600, units = "px", dpi = 310)

# ==============================================================================
# SECTION 6: Relative abundance barplot by taxonomic rank
# ==============================================================================

# --- Settings ---
hide_unclassified    <- FALSE
add_others_category  <- TRUE
min_abundance        <- 5   # minimum % to display a taxon individually

# Taxonomic rank (change to: "phylum", "class", "order", "family", "genus" and "species")
taxrank <- "phylum"

# Agglomerate and build proportions table
rank_seqtab     <- tax_glom(seqtab, taxrank = taxrank)
rank_counts_tab <- otu_table(rank_seqtab)
rank_tax_vec    <- as.vector(tax_table(rank_seqtab)[, taxrank])
colnames(rank_counts_tab) <- rank_tax_vec

rank_counts_tab2 <- otu_table(tax_glom(seqtab, taxrank = taxrank, NArm = FALSE))
rank_tax_vec2    <- as.vector(tax_table(tax_glom(seqtab, taxrank = taxrank, NArm = FALSE))[, taxrank])
colnames(rank_counts_tab2) <- rank_tax_vec2

unclassified_counts <- rowSums(rank_counts_tab2) - rowSums(rank_counts_tab)
rank_with_unclassified <- cbind(rank_counts_tab, "Unclassified" = unclassified_counts)

if (hide_unclassified) {
  rank_with_unclassified <- rank_with_unclassified[
    , !colnames(rank_with_unclassified) == "Unclassified"
  ]
}

asv_filt <- rank_with_unclassified[
  , colSums(rank_with_unclassified) * 100 / sum(rank_with_unclassified) >= min_abundance
]

if (add_others_category && min_abundance > 0) {
  asv_filt <- cbind(
    asv_filt,
    Others = rowSums(rank_with_unclassified[
      , !(colnames(rank_with_unclassified) %in% colnames(asv_filt))
    ])
  )
}

rank_taxa_proportions <- apply(asv_filt, 1, function(x) x / sum(x) * 100)
rank_taxa_proportions_non0 <- rank_taxa_proportions[rowSums(rank_taxa_proportions) > 0, ]

tidy_proportions <- reshape2::melt(
  rank_taxa_proportions_non0,
  value.name = "abundance",
  varnames   = c("taxa", "sampleName")
) %>%
  merge(data.frame(sample_data(seqtab)), by = "sampleName")

plot_labels <- sort(as.character(unique(tidy_proportions$taxa)))
if ("Others" %in% plot_labels)
  plot_labels <- c(plot_labels[plot_labels != "Others"], "Others")
if ("Unclassified" %in% plot_labels)
  plot_labels <- c(plot_labels[plot_labels != "Unclassified"], "Unclassified")

tidy_proportions$taxa <- factor(tidy_proportions$taxa, levels = plot_labels)
tidy_proportions <- tidy_proportions[rowSums(asv_filt) > 0, ]

tidy_proportions$sampleName2 <- factor(
  tidy_proportions$sampleName2,
  levels = c("Inoculum", "R1", "R2", "R3")
)

barplot_fig <- ggplot(tidy_proportions) +
  geom_col(aes(x = sampleName2, y = abundance, fill = taxa), position = "stack") +
  facet_grid(scales = "free", space = "fixed") +
  theme_bw(base_size = 13) +
  labs(fill = str_to_title(taxrank), y = "Abundance (%)") +
  theme(
    panel.grid      = element_line(colour = "white"),
    panel.spacing   = unit(0.4, "lines"),
    aspect.ratio    = 1.3,
    text            = element_text(family = "serif", color = "black", face = "plain"),
    strip.placement = "outside", strip.background = element_blank(),
    strip.text      = element_text(size = 13, face = "bold"),
    axis.title.x    = element_blank(),
    axis.text.x     = element_text(size = 13, angle = 0, hjust = 0.5, color = "black"),
    axis.text.y     = element_text(size = 13, hjust = 1),
    legend.text     = element_text(size = 12),
    legend.title    = element_text(size = 12),
    legend.key.size = unit(0.8, "lines")
  ) +
  scale_fill_manual(values = brewer.pal(n = 9, "PuOr")) +
  scale_y_continuous(
    expand = expansion(mult = c(0.006, 0.006)),
    labels = function(x) paste0(x, "%")
  ) +
  scale_x_discrete(expand = expansion(mult = c(0.16, 0.16)))

print(barplot_fig)

# ggsave("figures/barplot_phylum.tiff", plot = barplot_fig, height = 1300, width = 2700, units = "px", dpi = 300)

# ==============================================================================
# SECTION 7: Heatmap of relative abundance by taxonomic rank
# ==============================================================================

tidy_no_unclass <- tidy_proportions[
  !tidy_proportions$taxa %in% c("Unclassified", "Others"), ]

tidy_no_unclass$taxa <- factor(tidy_no_unclass$taxa, levels = rev(plot_labels))

heatmap_mat <- tidy_no_unclass %>%
  select(sampleName2, taxa, abundance) %>%
  pivot_wider(names_from = sampleName2, values_from = abundance, values_fill = 0) %>%
  column_to_rownames("taxa") %>%
  as.matrix()

purple_pastel <- colorRampPalette(c("#FDF4D9", "#E48751"))(100)

heatmap_fig <- pheatmap(
  heatmap_mat,
  color                     = purple_pastel,
  scale                     = "none",
  cluster_rows              = FALSE,
  cluster_cols              = FALSE,
  clustering_method         = "ward.D2",
  clustering_distance_rows  = "euclidean",
  clustering_distance_cols  = "euclidean",
  legend                    = TRUE,
  border_color              = "gray40",
  fontsize_row              = 12,
  fontsize_col              = 10,
  angle_col                 = 0
)

# ggsave("figures/heatmap_phylum.tiff", plot = heatmap_fig, height = 1300, width = 1490, units = "px", dpi = 350)
