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
cols <- c("A" = "#af0e2b", "B" = "#e4501e",
          "C" = "#4caae5", "D" = "#f3ac2c",
          "E" = "#57b741", "N" = "#8880be",
          "X" = "#81008b", "-" = "#aaaaaa")

#### PANEL A: Full Genome Hi-C ####
cat("Loading Hi-C data...\n")

hic_full <- import("analyses/hi-C/cooltools/nxAuaRhod1_1.mcool", format = 'mcool', resolution = 32000)
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
                     caption = FALSE) +
  scale_x_continuous(
    breaks = chrom_centers,
    labels = c("Chr 1", "Chr 2", "Chr 3", "Chr 4", "Chr 5", "Chr 6", "Chr X")
  ) +
  scale_y_reverse(
    breaks = chrom_centers,
    labels = c("Chr 1", "Chr 2", "Chr 3", "Chr 4", "Chr 5", "Chr 6", "Chr X")
  ) +
  labs(x = NULL, y = NULL, title = NULL) +
  theme(
    legend.position = "right",
    legend.justification = "center",
    axis.text = element_text(size = 10),
    axis.ticks = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  )

#### PANEL B: Chr 5 Hi-C Zoom ####
cat("Loading Chr 5 Hi-C zoom...\n")

# Read break sites first (needed for Panel B)
break_sites <- read_tsv("analyses/genome_features/elim_coords/nxAuaRhod1_1.break_sites.tsv") %>%
  mutate(chr = sub("SUPER_", "Chr ", chrom))

hic_chr5 <- import("analyses/hi-C/cooltools/nxAuaRhod1_1.mcool", format = 'mcool', resolution = 32000, focus = 'SUPER_5')

# Get Chr 5 GRS boundaries
chr5_break_sites <- break_sites %>%
  filter(chr == "Chr 5")

panel_b <- plotMatrix(hic_chr5,
                     show_grid = FALSE,
                     caption = FALSE) +
  geom_vline(data = chr5_break_sites, aes(xintercept = coordinate),
             linetype = "dotted", color = "black", linewidth = 0.8) +
  geom_hline(data = chr5_break_sites, aes(yintercept = coordinate),
             linetype = "dotted", color = "black", linewidth = 0.8) +
  scale_x_continuous(breaks = chr5_break_sites$coordinate,
                     labels = paste0(sprintf("%.2f", chr5_break_sites$coordinate / 1e6), "M")) +
  scale_y_reverse(breaks = chr5_break_sites$coordinate,
                  labels = paste0(sprintf("%.2f", chr5_break_sites$coordinate / 1e6), "M")) +
  labs(x = NULL, y = NULL, title = NULL) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 8.5),
    axis.ticks = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  )

#### GENOMIC FEATURES PANELS ####
cat("Loading genomic features data...\n")

# Read Nigon dictionary
nigon_dict <- read_tsv("metadata/gene2Nigon_busco20200927.tsv.gz",
                       col_types = c(col_character(), col_character()))

# Read BUSCO full table
busco <- read_tsv("analyses/genome_features/nemaChromQC/busco/nxAuaRhod1_1.curated_primary_nematoda_odb10_full_table.tsv",
                  comment = "#",
                  col_names = c("busco_id", "status", "sequence", "gene_start", "gene_end",
                                "strand", "score", "length", "orthodb_url", "description")) %>%
  filter(status == "Complete") %>%
  filter(sequence %in% paste0("SUPER_", c(1:6, "X"))) %>%
  left_join(nigon_dict, by = c("busco_id" = "Orthogroup")) %>%
  mutate(nigon = ifelse(is.na(nigon), "-", nigon),
         position = gene_start,
         chr = sub("SUPER_", "Chr ", sequence))

# Read additional data files
repsWind <- read_tsv("analyses/genome_features/nemaChromQC/red/nxAuaRhod1_1.curated_primary.red.bed.gz",
                     col_names = c("contig", "wStart", "wEnd", "value", "feat")) %>%
  mutate(chr = sub("SUPER_", "Chr ", contig))

gcWind <- read_tsv("analyses/genome_features/nemaChromQC/gc/nxAuaRhod1_1.curated_primary.gc.bed.gz",
                   col_names = c("contig", "wStart", "wEnd", "value")) %>%
  mutate(chr = sub("SUPER_", "Chr ", contig))

cov <- read_tsv("analyses/genome_features/nemaChromQC/nxAuaRhod1_1.regions.bed.gz",
                col_names = c("contig", "wStart", "wEnd", "value")) %>%
  mutate(chr = sub("SUPER_", "Chr ", contig))

# Read telomere data
# Functions to read PAFcopied from github.com/thackl/thacklr/R/read.R
read_paf <- function (file, max_tags = 20){
  col_names <- c("query_name", "query_length", "query_start", 
                 "query_end", "strand", "target_name", "target_length", 
                 "target_start", "target_end", "map_match", "map_length", 
                 "map_quality")
  col_types <- "ciiicciiiiin"
  
  if(max_tags > 0){
    col_names <- c(col_names, paste0("tag_", seq_len(max_tags)))
    col_types <- paste0(col_types, paste(rep("?", max_tags), collapse=""))
  }
  
  read_tsv(file, col_names = col_names, col_types = col_types) %>%
    tidy_paf_tags
}

tidy_paf_tags <- function(.data){
  tag_df <- tibble(.rows=nrow(.data))
  tag_types <- c()
  seen_empty_tag_col <- FALSE
  
  for (x in select(.data, starts_with("tag_"))){
    tag_mx <- str_split(x, ":", 3, simplify=T)
    tag_mx_nr <- na.omit(unique(tag_mx[,1:2]))
    if(nrow(tag_mx_nr) == 0){
      seen_empty_tag_col <- TRUE
      break;
    }
    tags <- tag_mx_nr[,1]
    tag_type <- tag_mx_nr[,2]
    names(tag_type) <- tags
    tag_types <- c(tag_types, tag_type)
    tag_types <- tag_types[unique(names(tag_types))]
    for (tag in tags){
      if(!has_name(tag_df, tag)){
        tag_df[[tag]] <- NA
      }
      tag_idx <- tag_mx[,1] %in% tag
      tag_df[[tag]][tag_idx] <- tag_mx[tag_idx,3]
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

block_mappings <- function(teloMappings, windwSize = 5e5){
  hqTeloMappings <- filter(teloMappings,
                           tp == "P",
                           map_length > query_length * 0.8,
                           target_length > windwSize * 2)
  duplicateIds <- filter(hqTeloMappings, duplicated(query_name)) %>%
    pull(query_name) %>% unique
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
                                                   target_start, target_end)) %>%
    arrange(target_name, rstart) %>%
    group_by(target_name, strand) %>%
    mutate(block = cumsum(c(1, diff(rstart) > 100))) %>%
    group_by(target_name, strand, block) %>%
    summarise(strand = unique(strand),
              regionStart = ifelse(strand == "-", getMode(target_start),
                                   getMode(target_start)),
              regionEnd = ifelse(strand == "+", getMode(target_end),
                                 getMode(target_end)),
              teloPos = ifelse(strand == "+", regionStart, regionEnd),
              target_length = unique(target_length),
              regSupport = n(),
              .groups = "drop") %>%
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
covForPlot <- filter(cov, chr %in% main_contigs,
                     value < quantile(value, 0.98))

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
  mutate(ints = as.numeric(as.character(cut(position,
                                            breaks = seq(0, max(position), windwSize),
                                            labels = seq(windwSize, max(position), windwSize)))),
         ints = ifelse(is.na(ints), max(ints, na.rm = T) + windwSize, ints)) %>%
  ungroup() %>%
  as.data.frame()

consUsco <- dplyr::count(busco_windowed, chr, ints, nigon) %>%
  mutate(chr = factor(chr, levels = mixedsort(unique(chr))))

# Panel C: Repeat density (moved up)
plReps <- mutate(allRepsForPlot, chr = factor(chr, levels = mixedsort(unique(chr)))) %>%
  ggplot(aes(x=wStart, y=value)) +
  geom_rect(data = chr_rects, aes(xmin = 0, xmax = size, ymin = 0, ymax = Inf), 
            fill = NA, color = "black", linewidth = 0.5, inherit.aes = FALSE) +
  geom_point(alpha=0.1, color = "gray") +
  geom_vline(data = internal_break_sites, aes(xintercept = coordinate),
             linetype = "dotted", color = "black", alpha = 0.5, inherit.aes = FALSE) +
  geom_point(data = internal_break_sites, aes(x = coordinate, y = Inf),
             shape = 25, fill = "black", size = 1, color = "black", inherit.aes = FALSE) +
  facet_grid(. ~ chr) +
  theme_bw() +
  scale_y_continuous("Repeat density",
                     labels = scales::percent,
                     position = "left",
                     limits = c(0, 1.05),
                     expand = c(0, 0)) +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale()),
                     limits = c(0, max_length),
                     expand = c(0, 0)) +
  ggtitle("") +
  theme(axis.title.x = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank(),
        plot.margin = margin(5, 75, 5, 10))

# Panel D: Nigon plot (moved down)
plNigon <- ggplot(consUsco, aes(fill=nigon, y=n, x=ints-(windwSize/2))) +
  geom_rect(data = chr_rects, aes(xmin = 0, xmax = size, ymin = 0, ymax = Inf), 
            fill = NA, color = "black", linewidth = 0.5, inherit.aes = FALSE) +
  geom_bar(position="stack", stat="identity") +
  geom_vline(data = internal_break_sites, aes(xintercept = coordinate),
             linetype = "dotted", color = "black", alpha = 0.5, inherit.aes = FALSE) +
  geom_point(data = internal_break_sites, aes(x = coordinate, y = Inf),
             shape = 25, fill = "black", size = 1, color = "black", inherit.aes = FALSE) +
  facet_grid(. ~ chr) +
  theme_bw() +
  scale_y_continuous("Nigon loci",
                     breaks = scales::pretty_breaks(4),
                     position = "left") +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale()),
                     limits = c(0, max_length),
                     expand = c(0, 0)) +
  scale_fill_manual(values = cols) +
  guides(fill = guide_legend(nrow = 1, title = "Nigon")) +
  ggtitle("") +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom",
        panel.grid = element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank(),
        plot.margin = margin(5, 75, 5, 10))

# Panel E: GC content
GC_upper_lim <- 0.45
plGC <- mutate(gcForPlot, chr = factor(chr, levels = mixedsort(unique(chr)))) %>%
  ggplot(aes(x=wStart, y=value)) +
  geom_rect(data = chr_rects, aes(xmin = 0, xmax = size, ymin = 0.20, ymax = GC_upper_lim), 
            fill = NA, color = "black", linewidth = 0.5, inherit.aes = FALSE) +
  geom_point(alpha=0.1, color = "gray") +
  geom_vline(data = internal_break_sites, aes(xintercept = coordinate),
             linetype = "dotted", color = "black", alpha = 0.5, inherit.aes = FALSE) +
  geom_point(data = internal_break_sites, aes(x = coordinate, y = GC_upper_lim),
             shape = 25, fill = "black", size = 1, color = "black", inherit.aes = FALSE) +
  facet_grid(. ~ chr) +
  theme_bw() +
  scale_y_continuous("GC", labels = scales::percent, position = "left",
                     limits = c(0.2, GC_upper_lim)) +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale()),
                     limits = c(0, max_length),
                     expand = c(0, 0)) +
  ggtitle("") +
  theme(axis.title.x = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank(),
        plot.margin = margin(5, 75, 5, 10))

# Panel F: Telomere reads
windwSize <- 5e5
teloBlocks <- block_mappings(teloMappings)

longSeqTeloMappings <- filter(teloMappings,
                              tp == "P",
                              map_length > query_length * 0.8,
                              target_length > windwSize * 2)

mappedTelo <- mutate(longSeqTeloMappings,
                     frac_target_start = (target_start / target_length)) %>%
  group_by(chr) %>%
  mutate(tReads = n()) %>%
  ungroup() %>%
  filter(tReads > 0.1 * max(tReads)) %>%
  select(chr, target_start, target_end, strand, frac_target_start)

plTelo <- mutate(mappedTelo, chr = factor(chr, levels = mixedsort(unique(chr)))) %>%
  ggplot(aes(x = target_start, color = strand)) +
  geom_rect(data = chr_rects, aes(xmin = 0, xmax = size, ymin = 0, ymax = Inf), 
            fill = NA, color = "black", linewidth = 0.5, inherit.aes = FALSE) +
  geom_histogram(bins = 100, alpha=1, position="identity", fill = "white") +
  geom_vline(data = internal_break_sites, aes(xintercept = coordinate),
             linetype = "dotted", color = "black", alpha = 0.5, inherit.aes = FALSE) +
  geom_point(data = internal_break_sites, aes(x = coordinate, y = Inf),
             shape = 25, fill = "black", size = 1, color = "black", inherit.aes = FALSE) +
  facet_grid(. ~ chr) +
  theme_bw() +
  scale_y_continuous("Telomeric reads", position = "left", expand = c(0, 0.1)) +
  scale_x_continuous("Position",
                     labels = label_number(scale_cut = cut_short_scale()),
                     limits = c(0, max_length),
                     expand = c(0, 0)) +
  guides(color = guide_legend(nrow = 1, title = "Telomere strand")) +
  ggtitle("") +
  theme(legend.position = "bottom",
        legend.justification = "left",
        panel.grid = element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank(),
        plot.margin = margin(5, 75, 5, 10))

# Panel G: Coverage
plCov <- mutate(covForPlot, chr = factor(chr, levels = mixedsort(unique(chr)))) %>%
  ggplot(aes(x=wStart, y=value)) +
  geom_rect(data = chr_rects, aes(xmin = 0, xmax = size, ymin = 0, ymax = Inf), 
            fill = NA, color = "black", linewidth = 0.5, inherit.aes = FALSE) +
  geom_point(alpha=0.1, color = "gray") +
  geom_vline(data = internal_break_sites, aes(xintercept = coordinate),
             linetype = "dotted", color = "black", alpha = 0.5, inherit.aes = FALSE) +
  geom_point(data = internal_break_sites, aes(x = coordinate, y = Inf),
             shape = 25, fill = "black", size = 1, color = "black", inherit.aes = FALSE) +
  facet_grid(. ~ chr) +
  theme_bw() +
  scale_y_continuous("Coverage", position = "left") +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale()),
                     limits = c(0, max_length),
                     expand = c(0, 0)) +
  ggtitle("") +
  theme(axis.title.x = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank(),
        plot.margin = margin(5, 75, 5, 10))

#### COMBINE ALL PANELS ####
cat("Combining all panels...\n")

# Combine Hi-C panels
hic_panels <- ggarrange(panel_a, panel_b,
                       nrow = 1,
                       labels = c("A", "B"),
                       align = "h")

# Combine genomic features panels (Telomeres moved to last)
genomic_panels <- ggarrange(plReps, plNigon, plGC, plCov, plTelo,
                           nrow = 5,
                           labels = c("C", "D", "E", "F", "G"),
                           align = "v")

# Combine all
combined_plot <- ggarrange(hic_panels, genomic_panels,
                          nrow = 2,
                          heights = c(1, 1.5))

# Display and save
print(combined_plot)

ggsave("report/figures/nxAuaRhod1_1_combined_figure.pdf",
       combined_plot, width = 16, height = 18, units = "in")

ggsave("report/figures/nxAuaRhod1_1_combined_figure.png",
       combined_plot, width = 16, height = 18, units = "in", dpi = 300)

cat("\nCombined figure saved to:\n")
cat("  report/figures/nxAuaRhod1_1_combined_figure.pdf\n")
cat("  report/figures/nxAuaRhod1_1_combined_figure.png\n")
