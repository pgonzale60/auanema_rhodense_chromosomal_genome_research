library(ggseqlogo)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(rtracklayer)
library(motifStack)

# Paths
MEME_OT <- "analyses/diminution/meme/meme_out/nxOscTipu1.1.meme.txt"
MEME_AR <- "analyses/diminution/meme/meme_out/meme.txt"
COV_BED <- "analyses/diminution/grs_visualization/SUPER_5_GRS_left_border.bed"
TELO_TSV <- "analyses/genome_features/elim_coords/nxAuaRhod1.pb.miltel.telomeric.tsv.gz"
FIMO_TSV <- "analyses/diminution/fimo_out/fimo.tsv.gz"
GENOME_GRS <- "analyses/diminution/nxAuaRhod1_1.GRS.bed"
GENOME_FA <- "nxAuaRhod1_1.primary.fa.gz"

#### DATA PREPARATION ####

# Load Regions
grs <- read.table(GENOME_GRS, col.names = c("chr", "start", "end"))
grs_gr <- GRanges(seqnames = grs$chr, ranges = IRanges(start = grs$start, end = grs$end))
# Left breaks (1bp at start)
breaks_left_gr <- GRanges(seqnames = grs$chr, ranges = IRanges(start = grs$start, end = grs$start + 1))
breaks_right_gr <- GRanges(seqnames = grs$chr, ranges = IRanges(start = grs$end - 1, end = grs$end))
breaks_all_gr <- c(breaks_left_gr, breaks_right_gr)

# Load FIMO hits
fimo_raw <- read_tsv(FIMO_TSV, comment = "#", show_col_types = FALSE) %>%
    filter(!is.na(score)) %>%
    arrange(desc(score))

fimo_gr_all <- GRanges(
    seqnames = fimo_raw$sequence_name,
    ranges = IRanges(start = fimo_raw$start, end = fimo_raw$stop),
    score = fimo_raw$score
)

#### PANEL B COORDINATE SELECTION (Fixed) ####

# Best motif hit at the SUPER_5 left GRS boundary (score 39.4, highest genome-wide)
# Somatic (high coverage) is LEFT of boundary; eliminated (low coverage) is RIGHT.
# GRS boundary: 16502174; motif hit: 16502162-16502190
target_chrom <- "SUPER_5"
target_boundary <- 16502174
motif_hit_start <- 16502162
motif_hit_end <- 16502190
x_start <- 16502153
x_end <- 16502197

#### PANEL A: ALIGNED MOTIF COMPARISON ####

motif_ot_stack <- motifStack::importMatrix(MEME_OT, format = "meme")
motif_ar_stack <- motifStack::importMatrix(MEME_AR, format = "meme")

spnames <- c("A. rhodense", "O. tipulae")
m_list <- list(motif_ar_stack[[1]], motif_ot_stack[[1]])

ord_motifs <- list()
for (i in 1:length(m_list)) {
    ord_motifs[[i]] <- motifStack::trimMotif(m_list[[i]], t = 0.4)
    ord_motifs[[i]]$name <- spnames[i]
}

pfmsAligned <- motifStack::DNAmotifAlignment(ord_motifs, rcpostfix = "")

motifs_aligned <- list()
for (i in 1:length(pfmsAligned)) {
    motifs_aligned[[pfmsAligned[[i]]$name]] <- as.matrix(as.data.frame(pfmsAligned[[i]]))
}
motifs_aligned <- motifs_aligned[spnames]

# Compute break position within aligned motif for panel A annotation
ar_orig_mat <- as.matrix(as.data.frame(motif_ar_stack[[1]]))
ar_trim_mat <- as.matrix(as.data.frame(ord_motifs[[1]])) # A. rhodense trimmed
ar_aligned_mat <- motifs_aligned[["A. rhodense"]]

left_trim_ar <- (which(sapply(
    seq_len(ncol(ar_orig_mat) - ncol(ar_trim_mat) + 1),
    function(ci) all(abs(ar_orig_mat[, ci] - ar_trim_mat[, 1]) < 0.001)
))[1]) - 1
if (is.na(left_trim_ar)) left_trim_ar <- 0

align_pad_ar <- (which(sapply(
    seq_len(ncol(ar_aligned_mat) - ncol(ar_trim_mat) + 1),
    function(ci) all(abs(ar_aligned_mat[, ci] - ar_trim_mat[, 1]) < 0.001)
))[1]) - 1
if (is.na(align_pad_ar)) align_pad_ar <- 0

# Break is at genomic position target_boundary; motif starts at motif_hit_start
break_pos_in_logo <- (target_boundary - motif_hit_start + 1.5) - left_trim_ar + align_pad_ar

# Color scheme for nucleotides
cs1 <- make_col_scheme(
    chars = c("A", "C", "G", "T"),
    cols = c("#009E73", "#0072B2", "#E69F00", "#D55E00")
)

panel_a <- suppressWarnings(ggseqlogo(motifs_aligned, ncol = 1, col_scheme = cs1))

# Move the vline to the first layer (behind the letters)
panel_a$layers <- c(
    geom_vline(xintercept = break_pos_in_logo, linetype = "dashed", color = "black", alpha = 0.8, linewidth = 0.5),
    panel_a$layers
)

panel_a <- panel_a +
    scale_x_continuous(breaks = seq(5, 30, by = 5)) +
    annotate("point", x = break_pos_in_logo, y = 2.3, shape = 25, fill = "black", size = 2, color = "black") +
    coord_cartesian(clip = "off") +
    labs(x = "Position (bp)") +
    theme_bw() +
    theme(
        strip.text = element_text(face = "italic", size = 11),
        strip.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size = 10, face = "plain"),
        axis.title.y = element_text(size = 10, face = "plain")
    )

#### PANEL B: BREAKSITE PRECISION ####

# Load coverage - EXACT window
df_cov <- read.table(COV_BED, col.names = c("chr", "start", "stop", "cov")) %>%
    filter(chr == target_chrom & stop >= x_start & start <= x_end)

# Load telomeres from TSV - breakage_score > 0.1
df_telo <- read_tsv(TELO_TSV,
    col_names = c(
        "chr", "start", "end", "score", "orientation",
        "gap_avg", "senses", "gaps", "avg_cov"
    ),
    show_col_types = FALSE
) %>%
    filter(chr == target_chrom & start >= x_start & start <= x_end) %>%
    filter(score > 0.1) %>%
    mutate(
        count = as.integer(sub("[-\\+]\\*", "", senses)),
        start = start - 1,
        end = end - 1
    ) # Minus 1 because this being 1-based whereas bed files being 0-indexed

# Motif hit coordinates
df_motif_hit <- data.frame(
    chr = target_chrom,
    start = motif_hit_start,
    stop = motif_hit_end,
    label = "Motif"
)

# Fetch nucleotide sequence using samtools
seq_cmd <- paste0("samtools faidx ", GENOME_FA, " ", target_chrom, ":", x_start, "-", x_end)
fasta_lines <- system(seq_cmd, intern = TRUE)
raw_seq <- paste(fasta_lines[-1], collapse = "")
seq_chars <- strsplit(raw_seq, "")[[1]]
df_seq <- data.frame(
    pos = x_start:(x_start + length(seq_chars) - 1),
    base = seq_chars
)

# Legend and Scales
panel_b <- ggplot() +
    # Motif hit highlight (sequence row only)
    geom_rect(data = df_motif_hit, aes(xmin = start - 0.5, xmax = stop + 0.5, ymin = -25, ymax = -5, fill = "Motif"), alpha = 0.5) +
    # Background coverage
    geom_rect(data = df_cov, aes(xmin = start - 0.5, xmax = stop + 0.5, ymin = 0, ymax = cov, fill = "All reads")) +
    # Telomere bars
    geom_rect(data = df_telo, aes(xmin = start - 0.5, xmax = start + 0.5, ymin = 0, ymax = count, fill = "Reads with\nsoft-clipped\ntelomeric\nrepeat")) +
    # Nucleotide sequence
    geom_text(data = df_seq, aes(x = pos, y = -15, label = base, color = base), size = 2, fontface = "bold") +
    # "Eliminated DNA" annotation (right of boundary = eliminated region)
    annotate("text", x = target_boundary + 13, y = 110, label = "Eliminated DNA", color = "grey40", fontface = "italic", size = 3) +
    annotate("segment", x = target_boundary, xend = x_end, y = 90, yend = 90, color = "grey40", arrow = arrow(ends = "both", length = unit(0.2, "cm"))) +
    # Legend and Scales
    scale_fill_manual(values = c(
        "All reads" = "grey85",
        "Motif" = "#166eb7",
        "Reads with\nsoft-clipped\ntelomeric\nrepeat" = "#E31A1C"
    ), breaks = c("All reads", "Reads with\nsoft-clipped\ntelomeric\nrepeat", "Motif")) +
    scale_color_manual(values = c(
        "A" = "#009E73", "C" = "#0072B2", "G" = "#E69F00", "T" = "#D55E00"
    ), guide = "none") +
    scale_x_continuous(expand = c(0, 0), breaks = target_boundary) +
    scale_y_continuous(expand = c(0, 0), limits = c(-30, 300)) +
    coord_cartesian(xlim = c(x_start, x_end)) +
    labs(x = "Position (bp)", y = "Coverage", fill = "") +
    theme_bw() +
    theme(
        legend.position = c(0.75, 0.75),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.spacing.y = unit(0, "pt"),
        legend.key.size = unit(2.2, "mm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        panel.grid = element_blank(),
        axis.title = element_text(size = 10, face = "plain")
    )


#### PANEL C: SINGLE-STRAND MOTIF SPECIFICITY ####

keep_idx <- c()
if (length(fimo_gr_all) > 0) {
    rem_gr <- fimo_gr_all
    while (length(rem_gr) > 0) {
        best_hit_c <- rem_gr[1]
        overlaps <- rem_gr %over% best_hit_c
        match_idx <- which(fimo_raw$sequence_name == seqnames(best_hit_c) &
            fimo_raw$start == start(best_hit_c) &
            fimo_raw$stop == end(best_hit_c) &
            fimo_raw$score == score(best_hit_c))[1]
        keep_idx <- c(keep_idx, match_idx)
        rem_gr <- rem_gr[!overlaps]
    }
}
fimo_filtered <- fimo_raw[keep_idx, ]

fimo_gr_final <- GRanges(
    seqnames = fimo_filtered$sequence_name,
    ranges = IRanges(start = fimo_filtered$start, end = fimo_filtered$stop),
    score = fimo_filtered$score
)

BREAK_SITES <- "analyses/genome_features/elim_coords/nxAuaRhod1_1.break_sites.tsv"
bps <- read.table(BREAK_SITES, header = TRUE)
bps_gr <- GRanges(seqnames = bps$chrom, ranges = IRanges(start = bps$coordinate - 50, end = bps$coordinate + 50))

fimo_gr_final$location <- "Retained Genome"
fimo_gr_final$location[fimo_gr_final %over% grs_gr] <- "Eliminated Region"
fimo_gr_final$location[fimo_gr_final %over% bps_gr] <- "Break Site"

df_spec <- as.data.frame(fimo_gr_final) %>%
    mutate(location = factor(location, levels = c("Break Site", "Eliminated Region", "Retained Genome")))

loc_counts <- df_spec %>%
    group_by(location) %>%
    summarise(n = n(), .groups = "drop")

# Consistent color palette across compartments (clearly distinguishing all three)
compartment_colors <- c(
    "Break Site" = "#D9381E",        # Vivid red/vermillion
    "Eliminated Region" = "#2171B5", # Medium royal blue
    "Retained Genome" = "#238B45"    # Forest green
)

min_bp_single <- min(df_spec$score[df_spec$location == "Break Site"])

panel_c <- ggplot(df_spec, aes(x = location, y = score)) +
    geom_jitter(aes(color = location), width = 0.28, alpha = 0.55, size = 1.4) +
    geom_hline(yintercept = min_bp_single, linetype = "dotted", color = "black", linewidth = 0.6) +
    geom_text(data = loc_counts, aes(x = location, y = 43, label = paste0("n=", n)), size = 2.8) +
    scale_color_manual(values = compartment_colors) +
    scale_x_discrete(labels = c(
        "Break Site" = "Break\nSite",
        "Eliminated Region" = "Eliminated\nRegion",
        "Retained Genome" = "Retained\nGenome"
    )) +
    scale_y_continuous(limits = c(0, 46), breaks = seq(0, 40, by = 10)) +
    labs(x = "Genomic Location", y = "Single-strand motif score (bits)", color = "") +
    theme_bw() +
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8.5),
        axis.text.y = element_text(size = 8.5),
        axis.title = element_text(size = 9.5, face = "plain")
    )


#### PANEL D: DUAL-STRAND PALINDROMIC SPECIFICITY (plus vs minus strand) ####

chroms <- paste0("SUPER_", c(1:6, "X"))
fimo_chr <- fimo_raw %>% filter(sequence_name %in% chroms)

pos_hits <- fimo_chr %>% filter(strand == "+")
neg_hits <- fimo_chr %>% filter(strand == "-")

# Canonical palindromic pairs: opposite strand, exact 2 bp coordinate stagger
pairs <- inner_join(pos_hits, neg_hits, by = "sequence_name", suffix = c("_plus", "_minus"), relationship = "many-to-many") %>%
    mutate(offset = start_plus - start_minus) %>%
    filter(offset == 2) %>%
    mutate(
        center = (start_plus + stop_plus + start_minus + stop_minus) / 4.0
    )

pairs_gr <- GRanges(pairs$sequence_name, IRanges(start = as.integer(pairs$center), width = 1))

pairs$comp <- "Retained Genome"
pairs$comp[pairs_gr %over% grs_gr] <- "Eliminated Region"
pairs$comp[pairs_gr %over% bps_gr] <- "Break Site"
pairs$comp <- factor(pairs$comp, levels = c("Break Site", "Eliminated Region", "Retained Genome"))

pair_counts <- pairs %>%
    group_by(comp) %>%
    summarise(n = n(), .groups = "drop")

poly_pass <- data.frame(
    x = c(0, 37.5, 42, 42, 0),
    y = c(37.5, 0, 0, 42, 42)
)

panel_d <- ggplot() +
    # Shaded grey passing zone (dual strand sum: score_plus + score_minus >= 37.5 bits)
    geom_polygon(data = poly_pass, aes(x = x, y = y), fill = "grey90", alpha = 0.7) +
    # Diagonal threshold line (score_plus + score_minus = 37.5)
    geom_abline(intercept = 37.5, slope = -1, linetype = "dotted", color = "grey40", linewidth = 0.5) +
    # Data points (consistent size = 1.4, alpha = 0.55 across all compartments)
    geom_point(
        data = pairs %>% filter(comp == "Retained Genome"),
        aes(x = score_plus, y = score_minus, color = comp),
        alpha = 0.55, size = 1.4
    ) +
    geom_point(
        data = pairs %>% filter(comp == "Eliminated Region"),
        aes(x = score_plus, y = score_minus, color = comp),
        alpha = 0.55, size = 1.4
    ) +
    geom_point(
        data = pairs %>% filter(comp == "Break Site"),
        aes(x = score_plus, y = score_minus, color = comp),
        alpha = 0.55, size = 1.4
    ) +
    scale_color_manual(
        values = compartment_colors,
        labels = c(
            "Break Site" = paste0("Break Sites (n=", pair_counts$n[pair_counts$comp == "Break Site"], ")"),
            "Eliminated Region" = paste0("Eliminated Region (n=", pair_counts$n[pair_counts$comp == "Eliminated Region"], ")"),
            "Retained Genome" = paste0("Retained Genome (n=", pair_counts$n[pair_counts$comp == "Retained Genome"], ")")
        )
    ) +
    scale_x_continuous(limits = c(0, 42), breaks = seq(0, 40, by = 10)) +
    scale_y_continuous(limits = c(0, 42), breaks = seq(0, 40, by = 10)) +
    labs(
        x = "Plus strand score (bits)",
        y = "Minus strand score (bits)",
        color = ""
    ) +
    theme_bw() +
    theme(
        legend.position = c(0.30, 0.90),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.key = element_rect(fill = "transparent", color = NA),
        legend.key.size = unit(2.5, "mm"),
        legend.text = element_text(size = 7.5),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8.5),
        axis.text.y = element_text(size = 8.5),
        axis.title = element_text(size = 9.5, face = "plain")
    )


#### COMBINE AND SAVE (2x2) ####

# Top row: coverage (panel_b) as A, motif alignment (panel_a) as B
# Bottom row: single-strand (panel_c) as C, dual-strand (panel_d) as D
final_plot <- ((panel_b | panel_a) / (panel_c | panel_d)) +
    plot_layout(heights = c(1, 1.1)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 11))

ggsave("report/figures/Figure_3.pdf", final_plot,
    width = 175, height = 150, units = "mm", dpi = 300, device = "pdf"
)
