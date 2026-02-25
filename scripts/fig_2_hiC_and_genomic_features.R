library(ggplot2)
library(tidyverse)
library(gtools)
library(scales)
library(ggpubr)
library(cowplot)
library(HiCExperiment)
library(HiContacts)
library(InteractionSet)

# Nigon color palette
cols <- c(
  "A" = "#af0e2b", "B" = "#e4501e",
  "C" = "#4caae5", "D" = "#f3ac2c",
  "E" = "#57b741", "N" = "#8880be",
  "X" = "#81008b", "-" = "#aaaaaa"
)

# Journal-standard theme (Current Biology / Cell Press)
theme_journal <- function(base_size = 7, base_family = "Helvetica") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = 8, face = "plain"),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 7),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      plot.title = element_text(size = 9, face = "bold", hjust = 0),
      plot.margin = margin(0, 0, 0, 0, "mm")
    )
}

#### PANEL A: Full Genome Hi-C ####
cat("Loading Hi-C data...\n")

hic_full <- import("analyses/hi-C/cooltools/nxAuaRhod1_1.mcool", format = "mcool", resolution = 32000)
main_chroms <- c("SUPER_1", "SUPER_2", "SUPER_3", "SUPER_4", "SUPER_5", "SUPER_6", "SUPER_X")

# Filter interactions to main chromosomes
hic_gi <- interactions(hic_full)
anchor1_chr <- as.character(seqnames(anchors(hic_gi, "first")))
anchor2_chr <- as.character(seqnames(anchors(hic_gi, "second")))
keep <- anchor1_chr %in% main_chroms & anchor2_chr %in% main_chroms
hic <- hic_full[keep]

# Calculate chromosome boundaries for labeling
chrom_lengths <- seqlengths(seqinfo(hic))[main_chroms]
chrom_ends <- cumsum(as.numeric(chrom_lengths))
chrom_starts <- c(0, chrom_ends[-length(chrom_ends)])
chrom_centers <- (chrom_starts + chrom_ends) / 2

# Panel A: Full genome Hi-C
panel_a <- plotMatrix(hic,
  show_grid = FALSE,
  caption = FALSE
) +
  scale_x_continuous(
    breaks = chrom_centers,
    labels = c("Chr 1", "Chr 2", "Chr 3", "Chr 4", "Chr 5", "Chr 6", "Chr X")
  ) +
  scale_y_reverse(
    breaks = chrom_centers,
    labels = c("Chr 1", "Chr 2", "Chr 3", "Chr 4", "Chr 5", "Chr 6", "Chr X")
  ) +
  labs(x = NULL, y = NULL, title = NULL) +
  theme_journal() +
  theme(
    legend.position = c(0.93, 0.8),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.key.size = unit(3, "mm"),
    panel.border = element_rect(color = "black", fill = NA),
    axis.ticks = element_blank()
  )

#### PANEL B: Chr 5 Hi-C Zoom ####
cat("Loading Chr 5 Hi-C zoom...\n")

# Read break sites first (needed for Panel B)
break_sites <- read_tsv("analyses/genome_features/elim_coords/nxAuaRhod1_1.break_sites.tsv") %>%
  mutate(chr = sub("SUPER_", "Chr ", chrom))

hic_chr5 <- import("analyses/hi-C/cooltools/nxAuaRhod1_1.mcool", format = "mcool", resolution = 32000, focus = "SUPER_5")

# Get Chr 5 GRS boundaries
chr5_break_sites <- break_sites %>%
  filter(chr == "Chr 5")

panel_b <- plotMatrix(hic_chr5,
  show_grid = FALSE,
  caption = FALSE
) +
  geom_vline(
    data = chr5_break_sites, aes(xintercept = coordinate),
    linetype = "dashed", color = "black", linewidth = 0.35
  ) +
  geom_hline(
    data = chr5_break_sites, aes(yintercept = coordinate),
    linetype = "dashed", color = "black", linewidth = 0.35
  ) +
  scale_x_continuous(
    breaks = chr5_break_sites$coordinate,
    labels = sprintf("%.2f", chr5_break_sites$coordinate / 1e6),
    guide = guide_axis(check.overlap = TRUE)
  ) +
  scale_y_reverse(
    breaks = chr5_break_sites$coordinate,
    labels = sprintf("%.2f", chr5_break_sites$coordinate / 1e6),
    guide = guide_axis(check.overlap = TRUE)
  ) +
  labs(x = "Position (Mb)", y = NULL, title = NULL) +
  theme_journal() +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA),
    axis.ticks = element_blank()
  )

#### GENOMIC FEATURES PANELS ####
cat("Loading genomic features data...\n")

# Read Nigon dictionary
nigon_dict <- read_tsv("metadata/gene2Nigon_busco20200927.tsv.gz",
  col_types = c(col_character(), col_character())
)

# Read BUSCO full table
busco <- read_tsv("analyses/genome_features/nemaChromQC/busco/nxAuaRhod1_1.curated_primary_nematoda_odb10_full_table.tsv",
  comment = "#",
  col_names = c(
    "busco_id", "status", "sequence", "gene_start", "gene_end",
    "strand", "score", "length", "orthodb_url", "description"
  )
) %>%
  filter(status == "Complete") %>%
  filter(sequence %in% paste0("SUPER_", c(1:6, "X"))) %>%
  left_join(nigon_dict, by = c("busco_id" = "Orthogroup")) %>%
  mutate(
    nigon = ifelse(is.na(nigon), "-", nigon),
    position = gene_start,
    chr = sub("SUPER_", "Chr ", sequence)
  )

# Read additional data files
repsWind <- read_tsv("analyses/genome_features/nemaChromQC/red/nxAuaRhod1_1.curated_primary.red.bed.gz",
  col_names = c("contig", "wStart", "wEnd", "value", "feat")
) %>%
  mutate(chr = sub("SUPER_", "Chr ", contig))

gcWind <- read_tsv("analyses/genome_features/nemaChromQC/gc/nxAuaRhod1_1.curated_primary.gc.bed.gz",
  col_names = c("contig", "wStart", "wEnd", "value")
) %>%
  mutate(chr = sub("SUPER_", "Chr ", contig))

cov <- read_tsv("analyses/genome_features/nemaChromQC/nxAuaRhod1_1.regions.bed.gz",
  col_names = c("contig", "wStart", "wEnd", "value")
) %>%
  mutate(chr = sub("SUPER_", "Chr ", contig))

# Read telomere data
# Functions to read PAFcopied from github.com/thackl/thacklr/R/read.R
read_paf <- function(file, max_tags = 20) {
  col_names <- c(
    "query_name", "query_length", "query_start",
    "query_end", "strand", "target_name", "target_length",
    "target_start", "target_end", "map_match", "map_length",
    "map_quality"
  )
  col_types <- "ciiicciiiiin"

  if (max_tags > 0) {
    col_names <- c(col_names, paste0("tag_", seq_len(max_tags)))
    col_types <- paste0(col_types, paste(rep("?", max_tags), collapse = ""))
  }

  read_tsv(file, col_names = col_names, col_types = col_types) %>%
    tidy_paf_tags()
}

tidy_paf_tags <- function(.data) {
  tag_df <- tibble(.rows = nrow(.data))
  tag_types <- c()
  seen_empty_tag_col <- FALSE

  for (x in select(.data, starts_with("tag_"))) {
    tag_mx <- str_split(x, ":", 3, simplify = T)
    tag_mx_nr <- na.omit(unique(tag_mx[, 1:2]))
    if (nrow(tag_mx_nr) == 0) {
      seen_empty_tag_col <- TRUE
      break
    }
    tags <- tag_mx_nr[, 1]
    tag_type <- tag_mx_nr[, 2]
    names(tag_type) <- tags
    tag_types <- c(tag_types, tag_type)
    tag_types <- tag_types[unique(names(tag_types))]
    for (tag in tags) {
      if (!has_name(tag_df, tag)) {
        tag_df[[tag]] <- NA
      }
      tag_idx <- tag_mx[, 1] %in% tag
      tag_df[[tag]][tag_idx] <- tag_mx[tag_idx, 3]
    }
  }

  tag_df <- tag_df %>%
    mutate_at(names(tag_types)[tag_types == "i"], as.integer) %>%
    mutate_at(names(tag_types)[tag_types == "f"], as.numeric)

  bind_cols(select(.data, -starts_with("tag_")), tag_df)
}

getMode <- function(x) {
  keys <- unique(x)
  keys[which.max(tabulate(match(x, keys)))]
}

block_mappings <- function(teloMappings, windwSize = 5e5) {
  hqTeloMappings <- filter(
    teloMappings,
    tp == "P",
    map_length > query_length * 0.8,
    target_length > windwSize * 2
  )
  duplicateIds <- filter(hqTeloMappings, duplicated(query_name)) %>%
    pull(query_name) %>%
    unique()
  unduplicated <- filter(hqTeloMappings, query_name %in% duplicateIds) %>%
    group_by(query_name) %>%
    mutate(onlyInFrag = max(target_length) < windwSize * 2) %>%
    filter((target_length < windwSize * 2 & onlyInFrag) | target_length >= windwSize * 2) %>%
    slice_min(query_start) %>%
    slice_max(map_quality) %>%
    slice_min(target_start) %>%
    ungroup()
  onlyOneMap <- filter(hqTeloMappings, !query_name %in% duplicateIds) %>%
    bind_rows(unduplicated)

  teloBlocks <- mutate(onlyOneMap, rstart = ifelse(strand == "+",
    target_start, target_end
  )) %>%
    arrange(target_name, rstart) %>%
    group_by(target_name, strand) %>%
    mutate(block = cumsum(c(1, diff(rstart) > 100))) %>%
    group_by(target_name, strand, block) %>%
    summarise(
      strand = unique(strand),
      regionStart = ifelse(strand == "-", getMode(target_start),
        getMode(target_start)
      ),
      regionEnd = ifelse(strand == "+", getMode(target_end),
        getMode(target_end)
      ),
      teloPos = ifelse(strand == "+", regionStart, regionEnd),
      target_length = unique(target_length),
      regSupport = n(),
      .groups = "drop"
    ) %>%
    arrange(target_name, teloPos) %>%
    select(target_name, target_length, strand, block, teloPos, regSupport)
  return(teloBlocks)
}

teloMappings <- read_paf("analyses/genome_features/nemaChromQC/teloMaps/nxAuaRhod1_1.curated_primary.teloMapped.paf.gz") %>%
  mutate(chr = sub("SUPER_", "Chr ", target_name))

# Chromosome lengths
chrom_lengths_bp <- c(24729150, 24040024, 23577664, 23111030, 21638500, 21466504, 26933942)
names(chrom_lengths_bp) <- c("Chr 1", "Chr 2", "Chr 3", "Chr 4", "Chr 5", "Chr 6", "Chr X")

# Get max length for x-axis
max_length <- max(chrom_lengths_bp)

# Filter data to main chromosomes
main_contigs <- c("Chr 1", "Chr 2", "Chr 3", "Chr 4", "Chr 5", "Chr 6", "Chr X")
allRepsForPlot <- filter(repsWind, chr %in% main_contigs)
gcForPlot <- filter(gcWind, chr %in% main_contigs)
covForPlot <- filter(
  cov, chr %in% main_contigs,
  value < quantile(value, 0.98)
)

# Create chr_rects for boundaries
chr_rects <- tibble(
  chr = names(chrom_lengths_bp),
  size = chrom_lengths_bp
)

# Filter break sites for internal boundaries only
internal_break_sites <- break_sites %>%
  left_join(chr_rects, by = "chr") %>%
  filter(coordinate > 10, coordinate < size - 10)

# Create windowed Nigon count
windwSize <- 5e5

busco_windowed <- busco %>%
  group_by(chr) %>%
  mutate(
    ints = as.numeric(as.character(cut(position,
      breaks = seq(0, max(position), windwSize),
      labels = seq(windwSize, max(position), windwSize)
    ))),
    ints = ifelse(is.na(ints), max(ints, na.rm = T) + windwSize, ints)
  ) %>%
  ungroup() %>%
  as.data.frame()

consUsco <- dplyr::count(busco_windowed, chr, ints, nigon) %>%
  filter(nigon != "-") %>%
  mutate(chr = factor(chr, levels = mixedsort(unique(chr))))

# Panel C: Repeat density (moved up)
plReps <- mutate(allRepsForPlot, chr = factor(chr, levels = mixedsort(unique(chr)))) %>%
  ggplot(aes(x = wStart, y = value)) +
  geom_rect(
    data = chr_rects, aes(xmin = 0, xmax = size, ymin = 0, ymax = Inf),
    fill = NA, color = "black", linewidth = 0.5, inherit.aes = FALSE
  ) +
  geom_point(alpha = 0.1, color = "gray") +
  geom_vline(
    data = internal_break_sites, aes(xintercept = coordinate),
    linetype = "dashed", color = "black", alpha = 0.5, linewidth = 0.35, inherit.aes = FALSE
  ) +
  facet_grid(. ~ chr) +
  theme_journal() +
  scale_y_continuous("Repeat density",
    labels = scales::percent,
    position = "left",
    limits = c(0, 1.05),
    expand = c(0, 0)
  ) +
  scale_x_continuous(
    breaks = seq(0, 25e6, 10e6),
    labels = label_number(scale = 1e-6),
    limits = c(0, max_length),
    expand = c(0, 0)
  ) +
  ggtitle("") +
  theme(
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# Panel D: Nigon plot (moved down)
plNigon <- ggplot(consUsco, aes(fill = nigon, y = n, x = ints - (windwSize / 2))) +
  geom_rect(
    data = chr_rects, aes(xmin = 0, xmax = size, ymin = 0, ymax = Inf),
    fill = NA, color = "black", linewidth = 0.5, inherit.aes = FALSE
  ) +
  geom_bar(position = "stack", stat = "identity") +
  geom_vline(
    data = internal_break_sites, aes(xintercept = coordinate),
    linetype = "dashed", color = "black", alpha = 0.5, linewidth = 0.35, inherit.aes = FALSE
  ) +
  facet_grid(. ~ chr) +
  theme_bw() +
  scale_y_continuous("Nigon loci",
    breaks = scales::pretty_breaks(4),
    position = "left"
  ) +
  scale_x_continuous(
    breaks = seq(0, 25e6, 10e6),
    labels = label_number(scale = 1e-6),
    limits = c(0, max_length),
    expand = c(0, 0)
  ) +
  scale_fill_manual(values = cols) +
  guides(fill = guide_legend(nrow = 1, title = "Nigon")) +
  ggtitle("") +
  theme_journal() +
  theme(
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(-2, 0, 0, 0)
  )

# Panel E: GC content - REMOVED

# Panel F: Telomere reads
windwSize <- 5e5
teloBlocks <- block_mappings(teloMappings)

longSeqTeloMappings <- filter(
  teloMappings,
  tp == "P",
  map_length > query_length * 0.8,
  target_length > windwSize * 2
)

mappedTelo <- mutate(longSeqTeloMappings,
  frac_target_start = (target_start / target_length)
) %>%
  group_by(chr) %>%
  mutate(tReads = n()) %>%
  ungroup() %>%
  filter(tReads > 0.1 * max(tReads)) %>%
  select(chr, target_start, target_end, strand, frac_target_start)

plTelo <- mutate(mappedTelo, chr = factor(chr, levels = mixedsort(unique(chr)))) %>%
  ggplot(aes(x = target_start, color = strand)) +
  geom_rect(
    data = chr_rects, aes(xmin = 0, xmax = size, ymin = 0, ymax = Inf),
    fill = NA, color = "black", linewidth = 0.5, inherit.aes = FALSE
  ) +
  geom_histogram(aes(y = after_stat(ifelse(count < 2, NA, count))), bins = 100, alpha = 1, position = "identity", fill = "white") +
  geom_vline(
    data = internal_break_sites, aes(xintercept = coordinate),
    linetype = "dashed", color = "black", alpha = 0.5, linewidth = 0.35, inherit.aes = FALSE
  ) +
  facet_grid(. ~ chr) +
  theme_bw() +
  scale_y_continuous("Telomeric\nreads", position = "left", expand = c(0, 0.1)) +
  scale_x_continuous(
    "Position (Mb)",
    breaks = seq(0, 25e6, 10e6),
    labels = label_number(scale = 1e-6),
    limits = c(0, max_length),
    expand = c(0, 0)
  ) +
  guides(color = guide_legend(nrow = 1, title = "Telomere strand")) +
  ggtitle("") +
  theme_journal() +
  theme(
    legend.position = "bottom",
    legend.justification = "left",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(-2, 0, 0, 0)
  )

# Panel G: Coverage
plCov <- mutate(covForPlot, chr = factor(chr, levels = mixedsort(unique(chr)))) %>%
  ggplot(aes(x = wStart, y = value)) +
  geom_rect(
    data = chr_rects, aes(xmin = 0, xmax = size, ymin = 0, ymax = Inf),
    fill = NA, color = "black", linewidth = 0.5, inherit.aes = FALSE
  ) +
  geom_point(alpha = 0.1, color = "gray") +
  geom_vline(
    data = internal_break_sites, aes(xintercept = coordinate),
    linetype = "dashed", color = "black", alpha = 0.5, linewidth = 0.35, inherit.aes = FALSE
  ) +
  facet_grid(. ~ chr) +
  theme_journal() +
  scale_y_continuous("Coverage", position = "left") +
  scale_x_continuous(
    breaks = seq(0, 25e6, 10e6),
    labels = label_number(scale = 1e-6),
    limits = c(0, max_length),
    expand = c(0, 0)
  ) +
  ggtitle("") +
  theme(
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  )

#### COMBINE ALL PANELS ####
cat("Combining all panels...\n")

# Combine Hi-C panels
hic_panels <- ggarrange(
  panel_a + theme(plot.margin = margin(r = 5, unit = "mm")),
  panel_b + theme(plot.margin = margin(l = 5, unit = "mm")),
  nrow = 1,
  labels = c("A", "B"),
  align = "h",
  font.label = list(size = 10, face = "bold", family = "Helvetica")
)

# Combine genomic features panels (Reordered: Nigon C, Reps D, Cov E, Telo F)
genomic_panels <- plot_grid(
  plNigon, plReps, plCov, plTelo,
  ncol = 1,
  align = "v",
  axis = "lr",
  labels = c("C", "D", "E", "F"),
  label_size = 10,
  label_fontfamily = "Helvetica",
  label_fontface = "bold",
  rel_heights = c(1.3, 1, 1, 1.3)
)

# Combine all
combined_plot <- plot_grid(
  hic_panels, genomic_panels,
  ncol = 1,
  rel_heights = c(1, 1.2)
) +
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"))

# Display and save
# Using cairo_pdf for better font embedding if needed, or standard pdf
ggsave("report/figures/fig_2_hiC_and_genomic_features.pdf",
  combined_plot,
  width = 170, height = 225, units = "mm", device = "pdf"
)

ggsave("report/figures/fig_2_hiC_and_genomic_features.png",
  combined_plot,
  width = 170, height = 225, units = "mm", dpi = 300
)

cat("\nFigure 2 saved to:\n")
cat("  report/figures/fig_2_hiC_and_genomic_features.pdf\n")
cat("  report/figures/fig_2_hiC_and_genomic_features.png\n")
