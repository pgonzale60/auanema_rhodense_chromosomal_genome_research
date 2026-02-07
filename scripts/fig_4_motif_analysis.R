library(ggseqlogo)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(rtracklayer)
library(motifStack)

# Paths
MEME_OT <- "analyses/diminution/meme/meme_out/nxOscTipu1.1.meme.txt"
MEME_AR <- "analyses/diminution/meme/meme_out/meme.txt"
COV_BED <- "analyses/diminution/grs_visualization/SUPER_X_GRS_borders.bed"
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

# Pick the best hit at a SUPER_X GRS boundary FOR WHICH WE HAVE COVERAGE DATA
# We know SUPER_X_GRS_borders.bed has data for 10706804 boundary
# Let's find motif hits overlapping ANY break site on SUPER_X
fimo_superx <- fimo_raw %>% filter(sequence_name == "SUPER_X")
fimo_superx_gr <- GRanges(fimo_superx$sequence_name, IRanges(fimo_superx$start, fimo_superx$stop))

# We specifically want a hit at a boundary that exists in our coverage file.
# The user wants "the left end of the GRS".
# 10706804 is a left end.
target_boundary <- 10706804
target_break_gr <- GRanges("SUPER_X", IRanges(target_boundary, target_boundary + 1))

hits_idx <- which(fimo_superx_gr %over% target_break_gr)

if (length(hits_idx) == 0) {
    # If no exact overlap, find the closest high scoring hit within 100bp
    hits_idx <- which(fimo_superx_gr %over% (target_break_gr + 100))
}

if (length(hits_idx) == 0) {
    stop("Could not find a motif hit near the specified GRS boundary on SUPER_X")
}

target_hits <- fimo_superx[hits_idx, ] %>% arrange(desc(score))
best_hit <- target_hits[1, ]

x_start <- 10706780
x_end <- 10706830
motif_hit_start <- best_hit$start
motif_hit_end <- best_hit$stop

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

# Color scheme for nucleotides
cs1 <- make_col_scheme(
    chars = c("A", "C", "G", "T"),
    cols = c("#009E73", "#0072B2", "#E69F00", "#D55E00")
)

panel_a <- ggseqlogo(motifs_aligned, ncol = 1, col_scheme = cs1) +
    theme_bw() +
    theme(
        strip.text = element_text(face = "italic", size = 12),
        strip.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line.y = element_line(color = "black")
    )

#### PANEL B: BREAKSITE PRECISION ####

# Load coverage - EXACT window
df_cov <- read.table(COV_BED, col.names = c("chr", "start", "stop", "cov")) %>%
    filter(chr == "SUPER_X" & stop >= x_start & start <= x_end)

# Load telomeres from TSV - breakage_score > 0.1
df_telo <- read_tsv(TELO_TSV,
    col_names = c(
        "chr", "start", "end", "score", "orientation",
        "gap_avg", "senses", "gaps", "avg_cov"
    ),
    show_col_types = FALSE
) %>%
    filter(chr == "SUPER_X" & start >= x_start & start <= x_end) %>%
    filter(score > 0.1) %>%
    mutate(
        count = as.integer(sub("[-\\+]\\*", "", senses)),
        start = start - 1,
        end = end - 1
    ) # Minus 1 because this being 1-based whereas bed files being 0-indexed

# Motif hit coordinates
df_motif_hit <- data.frame(
    chr = "SUPER_X",
    start = motif_hit_start,
    stop = motif_hit_end,
    label = "Motif"
)

# Fetch nucleotide sequence using samtools
seq_cmd <- paste0("samtools faidx ", GENOME_FA, " SUPER_X:", x_start, "-", x_end)
fasta_lines <- system(seq_cmd, intern = TRUE)
raw_seq <- paste(fasta_lines[-1], collapse = "")
seq_chars <- strsplit(raw_seq, "")[[1]]
df_seq <- data.frame(
    pos = x_start:(x_start + length(seq_chars) - 1),
    base = seq_chars
)

# Legend and Scales
panel_b <- ggplot() +
    # Motif hit highlight (background)
    geom_rect(data = df_motif_hit, aes(xmin = start - 0.5, xmax = stop + 0.5, ymin = -Inf, ymax = Inf, fill = "Motif"), alpha = 0.3) +
    # Background coverage
    geom_rect(data = df_cov, aes(xmin = start - 0.5, xmax = stop + 0.5, ymin = 0, ymax = cov, fill = "All reads")) +
    # Telomere bars
    geom_rect(data = df_telo, aes(xmin = start - 0.5, xmax = start + 0.5, ymin = 0, ymax = count, fill = "Reads with softclipped\ntelomeric repeat")) +
    # Nucleotide sequence
    geom_text(data = df_seq, aes(x = pos, y = -15, label = base, color = base), size = 3, fontface = "bold") +
    # "Eliminated DNA" annotation (to the right for left boundary)
    annotate("text", x = 10706815, y = 140, label = "Eliminated DNA", color = "grey40", fontface = "italic", size = 4) +
    annotate("segment", x = 10706805, xend = 10706830, y = 125, yend = 125, color = "grey40", arrow = arrow(ends = "both", length = unit(0.2, "cm"))) +
    # Legend and Scales
    scale_fill_manual(values = c(
        "All reads" = "grey85",
        "Motif" = "#166eb7",
        "Reads with softclipped\ntelomeric repeat" = "#E31A1C"
    ), breaks = c("All reads", "Reads with softclipped\ntelomeric repeat", "Motif")) +
    scale_color_manual(values = c(
        "A" = "#009E73", "C" = "#0072B2", "G" = "#E69F00", "T" = "#D55E00"
    ), guide = "none") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), limits = c(-30, 320)) +
    coord_cartesian(xlim = c(x_start, x_end)) +
    labs(x = "Position on chr X (bp)", y = "Count / Coverage", fill = "") +
    theme_bw() +
    theme(
        legend.position = c(0.75, 0.75),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_blank(),
        panel.grid = element_blank()
    )


#### PANEL C: MOTIF SPECIFICITY (FILTERED) ####

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

fimo_gr_final$location <- "Retained"
fimo_gr_final$location[fimo_gr_final %over% grs_gr] <- "Eliminated"
fimo_gr_final$location[fimo_gr_final %over% breaks_all_gr] <- "Telomere addition site"

df_spec <- as.data.frame(fimo_gr_final) %>%
    mutate(location = factor(location, levels = c("Telomere addition site", "Eliminated", "Retained")))

loc_counts <- df_spec %>%
    group_by(location) %>%
    summarise(n = n(), .groups = "drop")

panel_c <- ggplot(df_spec, aes(x = location, y = score)) +
    geom_jitter(aes(color = location), width = 0.3, alpha = 0.6) +
    geom_text(data = loc_counts, aes(x = location, y = 45, label = paste0("n=", n)), size = 3) +
    labs(x = "Location", y = "Motif match score", color = "") +
    theme_bw() +
    theme(legend.position = "none", panel.grid = element_blank())


#### COMBINE AND SAVE ####

final_plot <- (panel_a / panel_b / panel_c) +
    plot_layout(heights = c(1, 1, 1)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 16))

ggsave("report/figures/fig_4_motif_analysis.pdf", final_plot, width = 8, height = 12)
