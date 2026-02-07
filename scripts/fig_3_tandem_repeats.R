library(tidyverse)
library(GenomicRanges)
library(patchwork)

# Paths
FAI_FILE <- "analyses/genome_features/sequence_sizes/nxAuaRhod1.1.primary.fa.gz.fai"
MEMBERS_FILE <- "analyses/genome_features/repeats/TRF/trf_family_members.tsv"
GRS_FILE <- "analyses/genome_features/elim_coords/nxAuaRhod1_1.GRS.bed"
GENE_GFF <- "analyses/genes/GTFs/nxAuaRhod1_1.gff.gz"

# RNA source files
RNAMMER_GFF <- "analyses/genome_features/repeats/rnammer/nxAuaRhod1_1.rnammer.gff3.gz"
TRNASCAN_GFF <- "analyses/genome_features/repeats/tRNAscan/nxAuaRhod1_1.trnas.gff.gz"
INFERNAL_GFF <- "analyses/genome_features/repeats/infernal/nxAuaRhod1_1.infernal.gff.gz"

OUTPUT_PDF <- "report/figures/fig_3_tandem_repeats.pdf"
WINDOW_SIZE <- 100000

# 1. Load genome info
message("Loading genome info...")
seq_sizes <- read_tsv(FAI_FILE, col_names = c("Sequence", "size", "L1", "L2", "L3"), col_types = "ciiii") %>%
  select(Sequence, size) %>%
  filter(!grepl("MT|unloc|scaffold", Sequence))

gnm_gr <- Seqinfo(seqnames = seq_sizes$Sequence, seqlengths = seq_sizes$size, isCircular = rep(FALSE, nrow(seq_sizes)))
genome_windows <- tileGenome(gnm_gr, tilewidth = WINDOW_SIZE, cut.last.tile.in.chrom = TRUE)

# GRS regions
grs_bed <- read_tsv(GRS_FILE, col_names = c("Chrom", "Start", "End")) %>%
  filter(Chrom %in% seq_sizes$Sequence) %>%
  mutate(GRS_ID = paste0(Chrom, "_", row_number()))
grs_gr <- makeGRangesFromDataFrame(grs_bed, seqnames.field = "Chrom", keep.extra.columns = TRUE)

# 2. Load and Categorize Data
message("Loading gene regions...")
genes <- read_tsv(GENE_GFF, comment = "#", col_names = c("Chrom", "source", "type", "Start", "End", "score", "strand", "phase", "attributes")) %>%
  filter(type == "gene", Chrom %in% seq_sizes$Sequence) %>%
  mutate(chr = sub("SUPER_", "Chr ", Chrom))

message("Loading and categorizing structural RNAs...")
rrna_raw <- read_tsv(RNAMMER_GFF, comment = "#", col_names = c("Chrom", "src", "type", "Start", "End", "score", "strand", "phase", "attr")) %>%
  mutate(type = "rRNA")
trna_raw <- read_tsv(TRNASCAN_GFF, comment = "#", col_names = c("Chrom", "src", "type", "Start", "End", "score", "strand", "phase", "attr")) %>%
  mutate(type = "tRNA")
irna_raw <- read_tsv(INFERNAL_GFF, comment = "#", col_names = c("Chrom", "src", "type", "Start", "End", "score", "strand", "phase", "attr")) %>%
  mutate(type = case_when(
    grepl("rRNA", type, ignore.case = T) ~ "rRNA",
    grepl("tRNA", type, ignore.case = T) ~ "tRNA",
    grepl("U[12456]|snRNA|spliceosomal", type, ignore.case = T) ~ "snRNA",
    TRUE ~ "Other RNA"
  ))

# Resolve overlaps based on priority: rRNA > tRNA > snRNA > Other RNA
message("Resolving RNA overlaps...")
make_rna_gr <- function(df, t) {
  df_f <- df %>% filter(type == t, Chrom %in% seq_sizes$Sequence)
  if (nrow(df_f) == 0) {
    return(GRanges())
  }
  makeGRangesFromDataFrame(df_f, seqnames.field = "Chrom", start.field = "Start", end.field = "End", keep.extra.columns = T)
}

# Collect categories in priority order
r_gr <- GenomicRanges::reduce(make_rna_gr(bind_rows(rrna_raw, irna_raw), "rRNA"))
t_gr <- GenomicRanges::reduce(make_rna_gr(bind_rows(trna_raw, irna_raw), "tRNA"))
s_gr <- GenomicRanges::reduce(make_rna_gr(irna_raw, "snRNA"))
o_gr <- GenomicRanges::reduce(make_rna_gr(irna_raw, "Other RNA"))

# Priority-based subtraction (cumulative union)
mask <- r_gr
if (length(t_gr) > 0) t_gr <- GenomicRanges::setdiff(t_gr, mask)
mask <- GenomicRanges::union(mask, t_gr)
if (length(s_gr) > 0) s_gr <- GenomicRanges::setdiff(s_gr, mask)
mask <- GenomicRanges::union(mask, s_gr)
if (length(o_gr) > 0) o_gr <- GenomicRanges::setdiff(o_gr, mask)

# Re-assemble as Factored categories
rna_list <- list()
if (length(r_gr) > 0) rna_list$rRNA <- as_tibble(r_gr) %>% mutate(type = "rRNA")
if (length(t_gr) > 0) rna_list$tRNA <- as_tibble(t_gr) %>% mutate(type = "tRNA")
if (length(s_gr) > 0) rna_list$snRNA <- as_tibble(s_gr) %>% mutate(type = "snRNA")
if (length(o_gr) > 0) rna_list$Other_RNA <- as_tibble(o_gr) %>% mutate(type = "Other RNA")

rna_all <- bind_rows(rna_list) %>%
  dplyr::rename(Chrom = seqnames, Start = start, End = end) %>%
  mutate(
    chr = sub("SUPER_", "Chr ", as.character(Chrom)),
    type = factor(type, levels = c("rRNA", "tRNA", "snRNA", "Other RNA"))
  )

# 3. TR Selection (Relaxed Winners)
message("Calculating TR winners...")
members <- read_tsv(MEMBERS_FILE) %>% filter(Chrom %in% seq_sizes$Sequence)
members_gr <- makeGRangesFromDataFrame(members, seqnames.field = "Chrom", start.field = "Start", end.field = "End", keep.extra.columns = TRUE)
olaps <- findOverlaps(grs_gr, members_gr)
olap_df <- data.frame(
  GRS_ID = grs_gr$GRS_ID[queryHits(olaps)],
  Family_ID = members_gr$Family_ID[subjectHits(olaps)],
  Length = width(pintersect(grs_gr[queryHits(olaps)], members_gr[subjectHits(olaps)]))
)
priority_info <- olap_df %>%
  group_by(GRS_ID, Family_ID) %>%
  summarise(family_span = sum(Length), .groups = "drop_last") %>%
  group_by(GRS_ID) %>%
  filter(family_span == max(family_span)) %>%
  ungroup()

# Categorize GRS and families for color intercalation
grs_order <- grs_bed %>%
  group_by(Chrom) %>%
  mutate(grs_type = case_when(row_number() == 2 ~ "central", TRUE ~ "terminal")) %>%
  ungroup()
priority_info <- priority_info %>% left_join(grs_order %>% select(GRS_ID, grs_type), by = "GRS_ID")

# Known distinct terminal anchors (Autosome ends and X terminal GRS)
anchors <- c("P176_F001", "P351_F001", "P342_F001", "P347_F001")

# Group remaining families by their primary location
central_fams <- priority_info %>%
  filter(grs_type == "central") %>%
  pull(Family_ID) %>%
  unique()
terminal_fams <- priority_info %>%
  filter(grs_type == "terminal") %>%
  pull(Family_ID) %>%
  unique()

# Pools (excluding anchors)
central_only <- setdiff(central_fams, anchors)
terminal_only <- setdiff(terminal_fams, c(anchors, central_only))

# Sort pools by total genomic span
span_rank <- olap_df %>%
  group_by(Family_ID) %>%
  summarise(total = sum(Length))
central_only <- span_rank %>%
  filter(Family_ID %in% central_only) %>%
  arrange(desc(total)) %>%
  pull(Family_ID)
terminal_only <- span_rank %>%
  filter(Family_ID %in% terminal_only) %>%
  arrange(desc(total)) %>%
  pull(Family_ID)

# Intercalate: Separating central families with terminal ones in the color space
intercalated <- c()
for (i in 1:max(length(central_only), length(terminal_only))) {
  if (i <= length(central_only)) intercalated <- c(intercalated, central_only[i])
  if (i <= length(terminal_only)) intercalated <- c(intercalated, terminal_only[i])
}

# Final factor levels: Anchors first, then intercalated central/terminal
priority_families <- c(anchors, setdiff(intercalated, anchors))

# Ensure all unique priority families are included
all_priority <- priority_info %>%
  pull(Family_ID) %>%
  unique()
priority_families <- c(priority_families, setdiff(all_priority, priority_families))

total_cov <- GenomicRanges::coverage(makeGRangesFromDataFrame(members, seqnames.field = "Chrom", start.field = "Start", end.field = "End", seqinfo = gnm_gr))
total_binned <- binnedAverage(genome_windows, total_cov, "fraction") %>%
  as_tibble() %>%
  mutate(chr = sub("SUPER_", "Chr ", seqnames))

tr_plot_data <- tibble()
for (fam_id in priority_families) {
  fam_mem <- members %>% filter(Family_ID == fam_id)
  if (nrow(fam_mem) == 0) next
  fam_cov <- GenomicRanges::coverage(makeGRangesFromDataFrame(fam_mem, seqnames.field = "Chrom", start.field = "Start", end.field = "End", seqinfo = gnm_gr))
  binned <- binnedAverage(genome_windows, fam_cov, "fraction") %>%
    as_tibble() %>%
    mutate(Family_ID = fam_id)
  tr_plot_data <- bind_rows(tr_plot_data, binned)
}
tr_plot_data <- tr_plot_data %>% mutate(chr = sub("SUPER_", "Chr ", seqnames), Family_ID = factor(Family_ID, levels = priority_families))

# GRS layers
chr_rects <- seq_sizes %>% mutate(chr = sub("SUPER_", "Chr ", Sequence))
grs_regions_plot <- grs_bed %>% mutate(chr = sub("SUPER_", "Chr ", Chrom))
grs_boundaries <- grs_regions_plot %>% pivot_longer(cols = c(Start, End), names_to = "type", values_to = "position")

# Filter for internal boundaries
internal_grs_boundaries <- grs_boundaries %>%
  left_join(chr_rects, by = "chr") %>%
  filter(position > 10, position < size - 10)

# RNA window-based fractions (10 Kb windows)
message("Calculating RNA fractions per 10 Kb window...")
RNA_WINDOW_SIZE <- 10000
rna_windows <- tileGenome(gnm_gr, tilewidth = RNA_WINDOW_SIZE, cut.last.tile.in.chrom = TRUE)

rna_plot_data <- tibble()
rna_types <- c("rRNA", "tRNA", "snRNA", "Other RNA")

for (rna_type in rna_types) {
  rna_subset <- rna_all %>% filter(type == rna_type)
  if (nrow(rna_subset) == 0) next

  rna_gr_type <- makeGRangesFromDataFrame(rna_subset, seqnames.field = "Chrom", start.field = "Start", end.field = "End", seqinfo = gnm_gr)
  rna_cov <- GenomicRanges::coverage(GenomicRanges::reduce(rna_gr_type))
  binned <- binnedAverage(rna_windows, rna_cov, "fraction") %>%
    as_tibble() %>%
    mutate(RNA_Type = rna_type, fraction = pmax(fraction, 1e-4))
  rna_plot_data <- bind_rows(rna_plot_data, binned)
}

rna_plot_data <- rna_plot_data %>%
  mutate(
    chr = sub("SUPER_", "Chr ", seqnames),
    RNA_Type = factor(RNA_Type, levels = c("rRNA", "tRNA", "snRNA", "Other RNA")),
    fraction = pmin(fraction, 1.0)
  )

# Filter TRs based on rRNA fraction thresholds
message("Filtering TRs based on rRNA content...")
# Calculate rRNA fraction per 10 Kb window
rrna_per_window <- rna_plot_data %>%
  filter(RNA_Type == "rRNA") %>%
  group_by(seqnames, start, end) %>%
  summarise(rrna_fraction = sum(fraction), .groups = "drop") %>%
  mutate(high_rrna = rrna_fraction > 0.05)

# Identify which 10 Kb windows have high rRNA
high_rrna_windows <- rrna_per_window %>%
  filter(high_rrna) %>%
  select(seqnames, window_start = start, window_end = end)

# For each 100 Kb TR window, check if it overlaps with any high-rRNA 10 Kb window
tr_plot_data_filtered <- tr_plot_data %>%
  rowwise() %>%
  mutate(
    overlaps_high_rrna = any(
      high_rrna_windows$seqnames == seqnames &
        high_rrna_windows$window_start < end &
        high_rrna_windows$window_end > start
    )
  ) %>%
  ungroup() %>%
  filter(!overlaps_high_rrna) %>%
  select(-overlaps_high_rrna)

message("TR windows filtered (rRNA overlap): ", nrow(tr_plot_data) - nrow(tr_plot_data_filtered), " excluded")

# 4. Plots
theme_panel <- theme_bw() + theme(
  panel.grid = element_blank(),
  panel.border = element_blank(),
  strip.background = element_blank(),
  strip.text.y = element_text(angle = 0, face = "bold"),
  axis.title.x = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  plot.margin = margin(b = 2, t = 4, l = 5, r = 5)
)

# Panel A: TRs
pA <- ggplot() +
  geom_area(data = total_binned, aes(x = start, y = fraction), fill = "grey90", alpha = 0.8) +
  geom_rect(data = grs_regions_plot, aes(xmin = Start, xmax = End, ymin = 0, ymax = 1), fill = "wheat", alpha = 0.2) +
  geom_bar(data = tr_plot_data_filtered, aes(x = start, y = fraction, fill = Family_ID), stat = "identity", position = "stack", width = WINDOW_SIZE) +
  geom_vline(data = internal_grs_boundaries, aes(xintercept = position), linetype = "dotted", color = "black", alpha = 0.5) +
  geom_point(data = internal_grs_boundaries, aes(x = position, y = 1.05), shape = 25, fill = "black", size = 1, color = "black") +
  geom_rect(data = chr_rects, aes(xmin = 0, xmax = size, ymin = 0, ymax = 1), fill = NA, color = "black", linewidth = 0.5) +
  facet_grid(chr ~ .) +
  scale_x_continuous(labels = scales::label_number(scale = 1e-6, suffix = "M"), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0, 1), limits = c(0, 1.1), expand = c(0, 0)) +
  labs(y = "TR Fraction", x = "Position", fill = "Major TR Families") +
  theme_panel +
  theme(legend.position = "top", legend.text = element_text(size = 7), legend.key.size = unit(0.3, "cm"), axis.title.x = element_text(), axis.text.x = element_text(), axis.ticks.x = element_line()) +
  guides(fill = guide_legend(nrow = 2))

# Panel B: Genes
pB <- ggplot() +
  geom_rect(data = grs_regions_plot, aes(xmin = Start, xmax = End, ymin = 0, ymax = 1), fill = "wheat", alpha = 0.2) +
  geom_rect(data = genes, aes(xmin = Start, xmax = End, ymin = 0.1, ymax = 0.9), fill = "darkgreen") +
  geom_vline(data = internal_grs_boundaries, aes(xintercept = position), linetype = "dotted", color = "black", alpha = 0.5) +
  geom_point(data = internal_grs_boundaries, aes(x = position, y = 1.05), shape = 25, fill = "black", size = 1, color = "black") +
  geom_rect(data = chr_rects, aes(xmin = 0, xmax = size, ymin = 0, ymax = 1), fill = NA, color = "black", linewidth = 0.5) +
  facet_grid(chr ~ .) +
  scale_x_continuous(labels = scales::label_number(scale = 1e-6, suffix = "M"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1.1), expand = c(0, 0), breaks = NULL) +
  labs(y = "Genes", x = "Position") +
  theme_panel +
  theme(axis.title.x = element_text(), axis.text.x = element_text(), axis.ticks.x = element_line())

# Panel C: RNA fractions (per 10 Kb window, log scale)
pC <- ggplot() +
  geom_rect(data = grs_regions_plot, aes(xmin = Start, xmax = End, ymin = 1e-4, ymax = 0.95), fill = "wheat", alpha = 0.2) +
  geom_line(data = rna_plot_data, aes(x = start, y = fraction, color = RNA_Type), linewidth = 0.5) +
  geom_vline(data = internal_grs_boundaries, aes(xintercept = position), linetype = "dotted", color = "black", alpha = 0.5) +
  geom_point(data = internal_grs_boundaries, aes(x = position, y = 1), shape = 25, fill = "black", size = 1, color = "black") +
  geom_rect(data = chr_rects, aes(xmin = 0, xmax = size, ymin = 1e-4, ymax = 2), fill = NA, color = "black", linewidth = 0.5) +
  facet_grid(chr ~ .) +
  scale_x_continuous(labels = scales::label_number(scale = 1e-6, suffix = "M"), expand = c(0, 0)) +
  scale_y_log10(labels = scales::label_number()) +
  scale_color_brewer(palette = "Set1", drop = FALSE) +
  labs(y = "RNA Fraction (log10)", x = "Position", color = "RNA Type") +
  guides(color = guide_legend(override.aes = list(linewidth = 2))) +
  theme_panel +
  theme(legend.position = "bottom", axis.title.x = element_text(), axis.text.x = element_text(), axis.ticks.x = element_line())

# Combine
message("Combining panels...")
combined <- (pA) / (pB | pC) +
  plot_layout(heights = c(2, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 16))

ggsave(OUTPUT_PDF, combined, width = 16, height = 15)
message("Done!")
