library(ggseqlogo)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(rtracklayer)
library(motifStack)
library(ape)

# Paths
TREE_NWK <- "phylogenetic_analysis/data/ordered_cladogram_13taxa.nwk"
MEME_OT <- "analyses/diminution/meme/meme_out/nxOscTipu1.1.meme.txt"
MEME_AR <- "analyses/diminution/meme/meme_out/nxAuaRhod1_1.expanded_33bp_pal.meme.txt"
COV_BED <- "analyses/diminution/grs_visualization/SUPER_5_GRS_left_border.bed"
TELO_TSV <- "analyses/genome_features/elim_coords/nxAuaRhod1.pb.miltel.telomeric.tsv.gz"
FIMO_TSV <- "analyses/diminution/fimo_out/nxAuaRhod1_1_expanded_33bp_pal_fimo.tsv.gz"
GENOME_GRS <- "analyses/diminution/nxAuaRhod1_1.GRS.bed"
GENOME_FA <- "nxAuaRhod1_1.primary.fa.gz"
BREAK_SITES <- "analyses/genome_features/elim_coords/nxAuaRhod1_1.break_sites.tsv"

#### DATA PREPARATION (PANELS B, C, D) ####

# Load Regions
grs <- read.table(GENOME_GRS, col.names = c("chr", "start", "end"))
grs_gr <- GRanges(seqnames = grs$chr, ranges = IRanges(start = grs$start, end = grs$end))
# Left breaks (1bp at start)
breaks_left_gr <- GRanges(seqnames = grs$chr, ranges = IRanges(start = grs$start, end = grs$start + 1))
breaks_right_gr <- GRanges(seqnames = grs$chr, ranges = IRanges(start = grs$end - 1, end = grs$end))
breaks_all_gr <- c(breaks_left_gr, breaks_right_gr)

# Load FIMO hits (expanded 33-bp palindromic motif)
chroms <- paste0("SUPER_", c(1:6, "X"))
fimo_raw <- read_tsv(FIMO_TSV, comment = "#", show_col_types = FALSE) %>%
    filter(!is.na(score)) %>%
    filter(sequence_name %in% chroms) %>%
    arrange(desc(score))

fimo_gr_all <- GRanges(
    seqnames = fimo_raw$sequence_name,
    ranges = IRanges(start = fimo_raw$start, end = fimo_raw$stop),
    score = fimo_raw$score
)


#### PANEL A: CLADOGRAM OF PROGRAMMED DNA ELIMINATION ####

if (file.exists(TREE_NWK)) {
    tree <- read.tree(TREE_NWK)
} else {
    tree_text <- "(((((Auanema_rhodense:1,Oscheius_tipulae:1):1,((Caenorhabditis_auriculariae:1,Caenorhabditis_monodelphis:1):1,(Caenorhabditis_parvicauda:1,(Caenorhabditis_briggsae:1,Caenorhabditis_elegans:1):1):1):1):1,Mesorhabditis_belari:1):1,(Toxocara_canis:1,(Parascaris_univalens:1,(Ascaris_lumbricoides:1,Ascaris_suum:1):1):1):1):0.5,Trichinella_spiralis:0.5);"
    tree <- read.tree(text = tree_text)
}

n_tips <- length(tree$tip.label)
tip_order <- tree$tip.label

# Vertical spread between tips (step = 1.88 to align T. spiralis with Genomic Location)
step_y <- 1.88
y_vals <- seq(from = (n_tips - 1) * step_y, to = 0.0, by = -step_y)
y_tip <- setNames(y_vals, tip_order)

n_nodes <- tree$Nnode
node_y <- numeric(n_tips + n_nodes)
node_y[1:n_tips] <- y_tip[tree$tip.label]

tree_post <- reorder(tree, "postorder")
edge <- tree_post$edge

for (i in 1:nrow(edge)) {
    parent <- edge[i, 1]
    children <- edge[edge[, 1] == parent, 2]
    if (all(node_y[children] != 0)) {
        node_y[parent] <- mean(node_y[children])
    }
}
for (i in 1:nrow(edge)) {
    parent <- edge[i, 1]
    children <- edge[edge[, 1] == parent, 2]
    node_y[parent] <- mean(node_y[children])
}

root_node <- n_tips + 1
node_step <- numeric(n_tips + n_nodes)

get_depths <- function(node, cur_depth) {
    node_step[node] <<- cur_depth
    children <- tree$edge[tree$edge[, 1] == node, 2]
    for (ch in children) {
        get_depths(ch, cur_depth + 1)
    }
}
get_depths(root_node, 0)

max_depth <- max(node_step)
tip_x <- max_depth + 0.8

horiz_segs <- data.frame()
vert_segs <- data.frame()

for (p in (n_tips + 1):(n_tips + n_nodes)) {
    children <- tree$edge[tree$edge[, 1] == p, 2]
    y_ch <- node_y[children]
    vert_segs <- rbind(vert_segs, data.frame(
        x = node_step[p],
        xend = node_step[p],
        y = min(y_ch),
        yend = max(y_ch)
    ))
    for (ch in children) {
        x_target <- if (ch <= n_tips) tip_x else node_step[ch]
        horiz_segs <- rbind(horiz_segs, data.frame(
            x = node_step[p],
            xend = x_target,
            y = node_y[ch],
            yend = node_y[ch]
        ))
    }
}

label_map <- c(
    "Trichinella_spiralis" = "T. spiralis",
    "Parascaris_univalens" = "P. univalens",
    "Toxocara_canis" = "T. canis",
    "Ascaris_lumbricoides" = "A. lumbricoides",
    "Ascaris_suum" = "A. suum",
    "Mesorhabditis_belari" = "M. belari",
    "Auanema_rhodense" = "Auanema rhodense",
    "Oscheius_tipulae" = "O. tipulae",
    "Caenorhabditis_auriculariae" = "C. auriculariae",
    "Caenorhabditis_monodelphis" = "C. monodelphis",
    "Caenorhabditis_parvicauda" = "C. parvicauda",
    "Caenorhabditis_briggsae" = "C. briggsae",
    "Caenorhabditis_elegans" = "C. elegans"
)

tips_df <- data.frame(
    id = 1:n_tips,
    raw_name = tree$tip.label,
    label = label_map[tree$tip.label],
    x = tip_x,
    y = node_y[1:n_tips],
    is_focal = (tree$tip.label == "Auanema_rhodense"),
    stringsAsFactors = FALSE
)

half_step <- step_y / 2
shade_df <- data.frame(
    xmin = -0.5,
    xmax = tip_x + 6.6,
    ymin = c(y_vals[5] - half_step, y_vals[7] - half_step, y_vals[8] - half_step, y_vals[12] - half_step),
    ymax = c(y_vals[1] + half_step, y_vals[6] + half_step, y_vals[8] + half_step, y_vals[9] + half_step),
    category = c("Precise PDE", "No PDE", "Precise PDE", "Imprecise PDE"),
    stringsAsFactors = FALSE
)

bg_colors <- c(
    "Precise PDE" = "#edf7ee",
    "Imprecise PDE" = "#fdf2e9",
    "No PDE" = "#f3f0f7"
)

panel_cladogram <- ggplot() +
    geom_rect(data = shade_df,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = category),
              alpha = 0.95, inherit.aes = FALSE) +
    geom_segment(data = vert_segs, aes(x = x, xend = xend, y = y, yend = yend),
                 color = "#2b2b2b", linewidth = 0.6) +
    geom_segment(data = horiz_segs, aes(x = x, xend = xend, y = y, yend = yend),
                 color = "#2b2b2b", linewidth = 0.6) +
    geom_text(data = tips_df,
              aes(x = x + 0.25, y = y, label = label),
              fontface = "italic", hjust = 0, size = 3.05, color = "#222222") +
    scale_fill_manual(name = NULL, values = bg_colors,
                      breaks = c("Precise PDE", "Imprecise PDE", "No PDE"),
                      guide = guide_legend(
                          nrow = 1,
                          byrow = TRUE,
                          override.aes = list(fill = c("#c7e9c0", "#fdd0a2", "#dadaeb"))
                      )) +
    scale_x_continuous(limits = c(-0.5, tip_x + 6.6), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-0.3, y_vals[1] + step_y * 0.95), expand = c(0, 0)) +
    theme_void() +
    theme(
        legend.position = c(0.50, 0.98),
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.text = element_text(size = 7.6, color = "#1a1a1a"),
        legend.key.size = unit(3.0, "mm"),
        legend.key = element_rect(color = "#bbbbbb", linewidth = 0.3),
        legend.background = element_rect(fill = alpha("white", 0.95), color = "#cccccc", linewidth = 0.4),
        legend.margin = margin(t = 2, r = 4, b = 2, l = 4),
        plot.margin = margin(t = 1, r = 8, b = 0, l = 3)
    )


#### PANEL B: BREAKSITE PRECISION (LONG-READ COVERAGE) ####

target_chrom <- "SUPER_5"
target_boundary <- 16502174
motif_hit_start <- 16502159
motif_hit_end <- 16502191
x_start <- 16502150
x_end <- 16502198

df_cov <- read.table(COV_BED, col.names = c("chr", "start", "stop", "cov")) %>%
    filter(chr == target_chrom & stop >= x_start & start <= x_end)

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
    )

df_motif_hit <- data.frame(
    chr = target_chrom,
    start = motif_hit_start,
    stop = motif_hit_end,
    label = "Motif"
)

seq_cmd <- paste0("samtools faidx ", GENOME_FA, " ", target_chrom, ":", x_start, "-", x_end)
fasta_lines <- system(seq_cmd, intern = TRUE)
raw_seq <- paste(fasta_lines[-1], collapse = "")
seq_chars <- strsplit(raw_seq, "")[[1]]
df_seq <- data.frame(
    pos = x_start:(x_start + length(seq_chars) - 1),
    base = seq_chars
)

panel_b <- ggplot() +
    geom_rect(data = df_motif_hit, aes(xmin = start - 0.5, xmax = stop + 0.5, ymin = -25, ymax = -5, fill = "Motif"), alpha = 0.5) +
    geom_rect(data = df_cov, aes(xmin = start - 0.5, xmax = stop + 0.5, ymin = 0, ymax = cov, fill = "All reads")) +
    geom_rect(data = df_telo, aes(xmin = start - 0.5, xmax = start + 0.5, ymin = 0, ymax = count, fill = "Reads with\nsoft-clipped\ntelomeric\nrepeat")) +
    geom_text(data = df_seq, aes(x = pos, y = -15, label = base, color = base), size = 2, fontface = "bold") +
    annotate("text", x = target_boundary + 13, y = 110, label = "Eliminated DNA", color = "grey40", fontface = "italic", size = 3) +
    annotate("segment", x = target_boundary, xend = x_end, y = 90, yend = 90, color = "grey40", arrow = arrow(ends = "both", length = unit(0.2, "cm"))) +
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
        axis.title = element_text(size = 9.5, face = "plain")
    )


#### PANEL C: ALIGNED MOTIF COMPARISON ####

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

ar_orig_mat <- as.matrix(as.data.frame(motif_ar_stack[[1]]))
ar_trim_mat <- as.matrix(as.data.frame(ord_motifs[[1]]))
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

break_pos_in_logo <- (target_boundary - motif_hit_start + 0.5) - left_trim_ar + align_pad_ar

cs1 <- make_col_scheme(
    chars = c("A", "C", "G", "T"),
    cols = c("#009E73", "#0072B2", "#E69F00", "#D55E00")
)

panel_a <- suppressWarnings(ggseqlogo(motifs_aligned, ncol = 1, col_scheme = cs1))
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
        strip.text = element_text(face = "italic", size = 10),
        strip.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        axis.text.x = element_text(size = 8.5),
        axis.text.y = element_text(size = 8.5),
        axis.title.x = element_text(size = 9.5, face = "plain"),
        axis.title.y = element_text(size = 9.5, face = "plain")
    )


#### PANEL D: SINGLE-STRAND MOTIF SPECIFICITY ####

if (length(fimo_gr_all) > 0) {
    hits <- findOverlaps(fimo_gr_all, fimo_gr_all)
    has_earlier_overlap <- split(subjectHits(hits) < queryHits(hits), queryHits(hits))
    keep_hits <- !sapply(has_earlier_overlap, any)
    fimo_filtered <- fimo_raw[keep_hits, ]
} else {
    fimo_filtered <- fimo_raw
}

fimo_gr_final <- GRanges(
    seqnames = fimo_filtered$sequence_name,
    ranges = IRanges(start = fimo_filtered$start, end = fimo_filtered$stop),
    score = fimo_filtered$score
)

bps <- read.table(BREAK_SITES, header = TRUE)
bps_gr <- GRanges(seqnames = bps$chrom, ranges = IRanges(start = bps$coordinate - 25, end = bps$coordinate + 25))

fimo_gr_final$location <- "Retained Genome"
fimo_gr_final$location[fimo_gr_final %over% grs_gr] <- "Eliminated Region"
fimo_gr_final$location[fimo_gr_final %over% bps_gr] <- "Break Site"

df_spec <- as.data.frame(fimo_gr_final) %>%
    mutate(location = factor(location, levels = c("Break Site", "Eliminated Region", "Retained Genome")))

loc_counts <- df_spec %>%
    group_by(location) %>%
    summarise(n = n(), .groups = "drop")

compartment_colors <- c(
    "Break Site" = "#D9381E",
    "Eliminated Region" = "#2171B5",
    "Retained Genome" = "#238B45"
)

min_bp_single <- min(df_spec$score[df_spec$location == "Break Site"])

panel_c <- ggplot(df_spec, aes(x = location, y = score)) +
    geom_jitter(aes(color = location), width = 0.28, alpha = 0.55, size = 1.3) +
    geom_hline(yintercept = min_bp_single, linetype = "dotted", color = "black", linewidth = 0.6) +
    geom_text(data = loc_counts, aes(x = location, y = 43, label = paste0("n=", n)), size = 2.6) +
    scale_color_manual(values = compartment_colors) +
    scale_x_discrete(labels = c(
        "Break Site" = "Break\nSite",
        "Eliminated Region" = "Eliminated\nRegion",
        "Retained Genome" = "Retained\nGenome"
    )) +
    scale_y_continuous(limits = c(0, 46), breaks = seq(0, 40, by = 10)) +
    labs(x = "Genomic Location", y = "Motif score (bits)", color = "") +
    theme_bw() +
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 9, face = "plain")
    )


#### COMBINE AND SAVE (FOUR PANELS: A, B, C, D) ####

right_stack <- (panel_b / panel_a / panel_c) +
    plot_layout(heights = c(1, 1.05, 1.05))

final_plot <- (panel_cladogram | right_stack) +
    plot_layout(widths = c(0.95, 1.05)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 11))

ggsave("report/figures/Figure_3.pdf", final_plot,
    width = 185, height = 185, units = "mm", dpi = 300, device = "pdf"
)
ggsave("report/figures/Figure_3.png", final_plot,
    width = 185, height = 185, units = "mm", dpi = 300
)

cat("Successfully generated complete Figure 3 with four panels (A, B, C, D)!\n")
