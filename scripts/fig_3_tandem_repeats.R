library(tidyverse)
library(GenomicRanges)
library(patchwork)
library(ggnewscale)

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
    TRUE ~ "Other ncRNAs"
  ))

# Resolve overlaps based on priority: rRNA > tRNA > snRNA > Other ncRNAs
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
o_gr <- GenomicRanges::reduce(make_rna_gr(irna_raw, "Other ncRNAs"))

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
if (length(o_gr) > 0) rna_list$Other_RNA <- as_tibble(o_gr) %>% mutate(type = "Other ncRNAs")

rna_all <- bind_rows(rna_list) %>%
  dplyr::rename(Chrom = seqnames, Start = start, End = end) %>%
  mutate(
    chr = sub("SUPER_", "Chr ", as.character(Chrom)),
    type = factor(type, levels = c("rRNA", "tRNA", "snRNA", "Other ncRNAs"))
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
priority_families <- unique(c(priority_families, all_priority, "P348_F003", "P291_F001"))

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
tr_plot_data <- tr_plot_data %>%
  mutate(
    chr = sub("SUPER_", "Chr ", seqnames),
    Family_ID = factor(Family_ID, levels = priority_families)
  )

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
rna_types <- c("rRNA", "tRNA", "snRNA", "Other ncRNAs")

for (rna_type in rna_types) {
  rna_subset <- rna_all %>% filter(type == rna_type)
  if (nrow(rna_subset) == 0) next

  rna_gr_type <- makeGRangesFromDataFrame(rna_subset, seqnames.field = "Chrom", start.field = "Start", end.field = "End", seqinfo = gnm_gr)
  rna_cov <- GenomicRanges::coverage(GenomicRanges::reduce(rna_gr_type))
  binned <- binnedAverage(rna_windows, rna_cov, "fraction") %>%
    as_tibble() %>%
    mutate(RNA_Type = rna_type)
  rna_plot_data <- bind_rows(rna_plot_data, binned)
}

rna_plot_data <- rna_plot_data %>%
  mutate(
    chr = sub("SUPER_", "Chr ", seqnames),
    RNA_Type = factor(RNA_Type, levels = c("rRNA", "tRNA", "snRNA", "Other ncRNAs")),
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

# 4. TR naming scheme
tr_name_map <- c(
  # Autosomal terminal repeats (conserved across all autosomes)
  "P351_F001" = "AT-1",
  "P176_F001" = "AT-2",
  # Autosomal central repeats (chromosome-specific)
  "P348_F001" = "AC-I",
  "P399_F001" = "AC-II",
  "P167_F002" = "AC-III",
  "P332_F002" = "AC-IV",
  "P231_F004" = "AC-V",
  "P412_F002" = "AC-VI",
  # X chromosome repeats
  "P347_F001" = "XT-1",
  "P342_F001" = "XT-2",
  "P248_F001" = "XC",
  "P348_F003" = "XR-1",
  "P291_F001" = "XR-2"
)

tr_recode_map <- setNames(names(tr_name_map), unname(tr_name_map))
tr_plot_data_filtered <- tr_plot_data_filtered %>%
  mutate(Family_ID = forcats::fct_recode(Family_ID, !!!tr_recode_map))

# Reorder factor levels: named families in biological order, unnamed last
elim_named_order <- c(
  "AT-1", "AT-2",
  "AC-I", "AC-II", "AC-III", "AC-IV", "AC-V", "AC-VI",
  "XT-1", "XT-2", "XC"
)
soma_names <- c("XR-1", "XR-2")
all_levels <- levels(tr_plot_data_filtered$Family_ID)
unnamed <- setdiff(all_levels, c(elim_named_order, soma_names))
tr_plot_data_filtered <- tr_plot_data_filtered %>%
  mutate(Family_ID = factor(Family_ID,
    levels = c(intersect(elim_named_order, all_levels), unnamed, soma_names)
  ))

# 5. Plots
# Custom colors for Chr X somatic families
x_soma_colors <- c(
  "XR-1" = "#1331F5",
  "XR-2" = "#9E1F87"
)

# Split TR data for dual legends
tr_soma <- tr_plot_data_filtered %>% filter(Family_ID %in% names(x_soma_colors))
tr_elim <- tr_plot_data_filtered %>% filter(!(Family_ID %in% names(x_soma_colors)))

theme_panel <- theme_bw() + theme(
  panel.grid = element_blank(),
  panel.border = element_blank(),
  strip.background = element_blank(),
  strip.text.y = element_text(angle = 0, face = "plain", size = 7),
  axis.text.y = element_text(size = 7),
  axis.title = element_text(size = 8),
  axis.title.x = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  plot.margin = margin(b = 8, t = 4, l = 5, r = 5)
)

# Panel A: TRs
pA <- ggplot() +
  geom_area(data = total_binned, aes(x = start, y = fraction), fill = "grey85", alpha = 0.8) +
  geom_rect(data = grs_regions_plot, aes(xmin = Start, xmax = End, ymin = 0, ymax = 1), fill = "grey70", alpha = 0.3) +
  # Layer 1: Eliminated regions
  geom_bar(data = tr_elim, aes(x = start, y = fraction, fill = Family_ID), stat = "identity", position = "stack", width = WINDOW_SIZE) +
  scale_fill_hue() +
  guides(fill = guide_legend(title = "Major TR families in eliminated regions", nrow = 2, order = 1, title.position = "top")) +
  # Layer 2: Somatic-retained (Chr X)
  new_scale_fill() +
  geom_bar(data = tr_soma, aes(x = start, y = fraction, fill = Family_ID), stat = "identity", position = "stack", width = WINDOW_SIZE) +
  scale_fill_manual(
    values = x_soma_colors,
    labels = c("XR-1", "XR-2")
  ) +
  guides(fill = guide_legend(title = "Major somatic-retained\nTR families", nrow = 1, order = 2, title.position = "top")) +
  # Annotations
  geom_vline(data = internal_grs_boundaries, aes(xintercept = position), linetype = "dashed", color = "black", alpha = 0.8, linewidth = 0.5) +
  geom_point(data = internal_grs_boundaries, aes(x = position, y = 1.05), shape = 25, fill = "black", size = 2, color = "black") +
  geom_rect(data = chr_rects, aes(xmin = 0, xmax = size, ymin = 0, ymax = 1), fill = NA, color = "black", linewidth = 0.5) +
  facet_grid(chr ~ .) +
  scale_x_continuous(labels = scales::label_number(scale = 1e-6), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0, 1), limits = c(0, 1.1), expand = c(0, 0)) +
  labs(y = "TR Fraction", x = "Position (Mb)") +
  theme_panel +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.4, "cm"),
    axis.title.x = element_text(),
    axis.text.x = element_text(),
    axis.ticks.x = element_line()
  )

# Panel B: Genes
pB <- ggplot() +
  geom_rect(data = grs_regions_plot, aes(xmin = Start, xmax = End, ymin = 0, ymax = 1), fill = "grey70", alpha = 0.3) +
  geom_rect(data = genes, aes(xmin = Start, xmax = End, ymin = 0.1, ymax = 0.9), fill = "darkgreen") +
  geom_vline(data = internal_grs_boundaries, aes(xintercept = position), linetype = "dashed", color = "black", alpha = 0.8, linewidth = 0.5) +
  geom_point(data = internal_grs_boundaries, aes(x = position, y = 1.05), shape = 25, fill = "black", size = 2, color = "black") +
  geom_rect(data = chr_rects, aes(xmin = 0, xmax = size, ymin = 0, ymax = 1), fill = NA, color = "black", linewidth = 0.5) +
  facet_grid(chr ~ .) +
  scale_x_continuous(labels = function(x) ifelse(x == 0, "", scales::label_number(scale = 1e-6)(x)), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1.1), expand = c(0, 0), breaks = NULL) +
  coord_cartesian(clip = "off") +
  labs(y = "Protein-coding genes", x = "Position (Mb)") +
  theme_panel +
  theme(axis.title.x = element_text(), axis.text.x = element_text(), axis.ticks.x = element_line())

# Panel C: RNA fractions (per 10 Kb window, stacked histogram)
pC <- ggplot() +
  geom_rect(data = grs_regions_plot, aes(xmin = Start, xmax = End, ymin = 0, ymax = 1), fill = "grey70", alpha = 0.3) +
  geom_bar(data = rna_plot_data, aes(x = start, y = fraction, fill = RNA_Type), stat = "identity", position = "stack", width = RNA_WINDOW_SIZE) +
  geom_vline(data = internal_grs_boundaries, aes(xintercept = position), linetype = "dashed", color = "black", alpha = 0.8, linewidth = 0.5) +
  geom_point(data = internal_grs_boundaries, aes(x = position, y = 1.05), shape = 25, fill = "black", size = 2, color = "black") +
  geom_rect(data = chr_rects, aes(xmin = 0, xmax = size, ymin = 0, ymax = 1), fill = NA, color = "black", linewidth = 0.5) +
  facet_grid(chr ~ .) +
  scale_x_continuous(labels = function(x) ifelse(x == 0, "", scales::label_number(scale = 1e-6)(x)), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0, 1), limits = c(0, 1.1), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  scale_fill_brewer(palette = "Set1", drop = FALSE) +
  labs(y = "RNA Fraction", x = "Position (Mb)", fill = "RNA Type") +
  guides(fill = guide_legend(nrow = 1, title.position = "top")) +
  theme_panel +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.4, "cm"),
    legend.spacing.x = unit(0.2, "cm"),
    plot.margin = margin(b = 8, t = 4, l = 5, r = 5),
    axis.title.x = element_text(),
    axis.text.x = element_text(),
    axis.ticks.x = element_line()
  )

# Combine
message("Combining panels...")
bottom_row <- (pB | pC) + plot_layout(widths = c(1, 1.5))
combined <- (pA) / bottom_row +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 10))

ggsave(OUTPUT_PDF, combined, width = 170, height = 230, units = "mm", dpi = 300, device = "pdf")
message("Done!")
