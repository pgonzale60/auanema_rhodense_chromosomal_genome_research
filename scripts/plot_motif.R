# motifStack is hard to install but enables motif alignment
# https://www.bioconductor.org/packages/release/bioc/vignettes/motifStack/inst/doc/motifStack_HTML.html
# mamba create -n motifStack -c bioconda -c conda-forge bioconductor-motifstack
# Plotting is better done ggseqlogo because motifStack changes the bitscore
# ggseqlogo https://doi.org/10.1093/bioinformatics/btx469
library(ggseqlogo)
library(ggplot2)
library(rtracklayer)
library(ape)
library(phytools)
library(scales)
library(tidyverse)
library(cowplot)
library(motifStack, lib.loc = "/Users/pg17/miniconda3/envs/motifStack_w_ape/lib/R/library/")


meme_files <-  dir("../multiOscheius/analyses/paper2/prelim_20220527/motif/meme/raw","nx.*meme.txt$", 
                   full.names = TRUE) %>% c(
  dir("../multiOscheius/analyses/paper1/prelim_to_20220220/motif/meme/curated_break_sites", "nxOscTipu1.1.meme.txt$", 
      full.names = TRUE) )
motifs <- motifStack::importMatrix(meme_files, format = "meme")
# Remove the additional APS4 motifs (these have less than 1e-10
# and show no unique occurrence in each break site)
motifs <- motifs[c(1,4)]

# Extract and assign metadata to the motifs
meme_names <- names(motifs)
meme_n_sites <- sub(".+sites[\\.]+([0-9]+)\\..+", "\\1", meme_names)
meme_evals <- sub(".+E.value\\.+([0-9e\\.]+)", "\\1", meme_names)
strains <- sub(".+nx(.+).meme.txt", "\\1", meme_files)
spnames <- c("A. rhodense", "O. tipulae")

# create trimmed motifs in desired order for alignment
ord_motifs <- list()
for (strain_num in 1:length(motifs)) {
  strain <- strains[strain_num]
  ord_motifs[[strain]] <- motifStack::trimMotif(motifs[[strain_num]], t=0.4)
  ord_motifs[[strain]]$name <- spnames[strain_num]
  ord_motifs[[strain]]$color <- motifStack::colorset(colorScheme='blindnessSafe')
  
}

pfmsAligned <- motifStack::DNAmotifAlignment(ord_motifs, rcpostfix = "")


# Create custom colour scheme
cs1 = make_col_scheme(chars=c('A', 'C', 'G', 'T'),
                      cols=c('#009E73', '#0072B2', '#E69F00', '#D55E00'))

# Extract the frequency matrices for ggseqlogo
motifs <- list()
for (i in 1:2) {
  motifs[[pfmsAligned[[i]]$name]] <- as.matrix(as.data.frame(pfmsAligned[[i]]))
}

# plot logos
main_motif_plot <- ggseqlogo(motifs, ncol=1, col_scheme=cs1) +
  theme(strip.text = element_text(face = "italic") #,
        # text = element_text(family = "Arial")
        )

# ggsave(filename = "report/figures/Arhod_motif.jpg", width = 6,
#        height = 4)




#### Plot motif specificity ####
### Core tibble creation ####

# TeloCoords
aua_telomil <- read_tsv("analyses/genome_features/elim_coords/nxAuaRhod1.pb.miltel.telomeric.tsv.gz",
                        col_names = c("Sequence", "start", "end", "score",
                                      "orientation", "telomere_gap_average", 
                                      "telomere_senses", "telomere_gaps",
                                      "average_coverage")) %>%
  mutate(telocov = as.integer(sub("[-\\+]\\*", "", telomere_senses)))

# Sequence sizes
fai_files <- list.files("analyses/genome_features/sequence_sizes/",
                        full.names = T, pattern = "1_1.primary.fa.gz.fai")
names(fai_files) <- make.names(sub(".+//(.+).primary.fa.gz.fai", "\\1", fai_files))
seq_sizes <- map_df(fai_files, read_tsv,
                    col_names = c("Sequence", "size", "L1", "L2", "L3"),
                    col_types = c("ciiii"),
                    .id = "assembly") %>%
  select(assembly, Sequence, size) %>%
  filter(!grepl("MT", Sequence)) %>%
  select(-assembly)

# GRS files
GRS_files <- list.files("analyses/genome_features/elim_coords/",
                        full.names = T, pattern = "nxAua.*.GRS.bed$")
names(GRS_files) <- make.names(sub(".+//(.+).GRS.bed", "\\1", GRS_files))
GRS <- map_df(GRS_files, read_tsv,
              col_names = c("Sequence", "start", "end"),
              col_types = c("cii"),
              .id = "assembly") %>%
  select(-assembly)

lGRS <- pivot_longer(GRS, !Sequence, names_to = "varia", values_to = "POS")


# Fimo sites
fimofiles <- list.files("../multiOscheius/analyses/paper2/prelim_20220527/motif/fimo/",
                        full.names = T, pattern = "fimo.tsv.gz")
names(fimofiles) <- make.names(sub(".+//(.+).fimo.tsv", "\\1", fimofiles))
fimo <- map_df(fimofiles, read_tsv,
               comment = "#",
               .id = "assembly") %>%
  filter(!grepl("MT", sequence_name)) %>%
  mutate(telo_only = ifelse( grepl("TTAGGCTTAGGCTTAGGCTTAGGCTT", matched_sequence,
                                   fixed = T), T, F),
         motif_id = ifelse(is.na(motif_id), "N", motif_id),
         motif_alt_id = ifelse(is.na(motif_alt_id), "MEME-1", motif_alt_id)) #%>%
# Make sure only highest scoring strand, when more than one is present, is kept
# group_by(assembly, motif_id, sequence_name, start) %>%
# slice_min(`p-value`, with_ties = F) %>%
# ungroup()

count(fimo, matched_sequence) %>%
  arrange(desc(n)) # No clear outliers


#### Process tibbles ####
elim_seqs <- filter(seq_sizes, grepl("loc|scaffold", Sequence)) %>%
  mutate(start = 1, end = size-1) %>%
  select(-size) %>%
  rbind(GRS)

### Genomic ranges ####
gnm_gr <- Seqinfo(seqnames = seq_sizes$Sequence, seqlengths=seq_sizes$size,
                  isCircular = rep(F, nrow(seq_sizes)), genome = "Auanema")
# GRS regions
GRS_gr <- GRanges(seqnames = Rle(elim_seqs$Sequence), 
                  ranges = IRanges(start = elim_seqs$start,
                                   end = elim_seqs$end), 
                  strand = "*")

# Rename sequences to include assembly name
fimo_gr <- mutate(fimo, chrom=sequence_name,
                  fimo_ID = paste0(chrom, start)) %>%
  select(chrom,
         start,
         end = stop,
         score = score,
         pval = `p-value`,
         qval = `q-value`,
         fimo_ID,
         motif_alt_id = motif_alt_id,
         telo_only) %>% 
  makeGRangesFromDataFrame(seqinfo = gnm_gr, keep.extra.columns = T)

# 
break_sites_gr <- mutate(lGRS, chrom=Sequence,
                         end = POS + 1) %>%
  select(chrom,
         start = POS,
         end) %>%
  makeGRangesFromDataFrame(seqinfo = gnm_gr, keep.extra.columns = T) %>%
  GenomicRanges::trim()

telocoord_gr <- rename(aua_telomil, chrom=Sequence) %>%
  select(chrom,
         start,
         end,
         orientation,
         telocov) %>%
  makeGRangesFromDataFrame(seqinfo = gnm_gr, keep.extra.columns = T) 

#### Intersect ranges ####

# Many break sites have more than one motif hit
# Sort by score; keep the matches with highest score
fimo_gr <- fimo_gr[order(fimo_gr$score, decreasing = T)]
fimo_m1_gr <- fimo_gr[fimo_gr$motif_alt_id == "MEME-1"]
highest_scoring_non_ov_ranges_idx <- findOverlaps(fimo_m1_gr, fimo_m1_gr) %>% as_tibble() %>%
  arrange(queryHits, subjectHits) %>%
  filter(!duplicated(queryHits)) %>% pull(subjectHits) %>%
  unique()
highest_fimo_gr <- fimo_m1_gr[highest_scoring_non_ov_ranges_idx]

# tmp_gr <- highest_fimo_gr[highest_fimo_gr$fimo_ID %in% OscSper_merge_motifs$highest_m1_fimo_gr$fimo_ID] 
# OscSper_merge_motifs$highest_m1_fimo_gr$score[match(tmp_gr$fimo_ID, OscSper_merge_motifs$highest_m1_fimo_gr$fimo_ID)] %>%
#   as_tibble() %>% View
highest_fimo_gr$atBreakSite <- highest_fimo_gr %over% break_sites_gr
highest_fimo_gr$atBreakSite <- ifelse(highest_fimo_gr$atBreakSite,
                                      "breaksite",
                                      ifelse(highest_fimo_gr %over% GRS_gr,
                                             "eliminated", "core")
)
highest_fimo_gr$teloSupport <- highest_fimo_gr %over% telocoord_gr



# Data for plots
pd_fimo <- as_tibble(highest_fimo_gr) %>%
  filter(score > 0)

filter(pd_fimo, atBreakSite == "core", score > 20)

# filter(!telo_only)
tb_n <- group_by(pd_fimo, atBreakSite) %>%
  summarise(n = n(), .groups = "drop")


# Boxplots comparison
motif_specificity <- ggplot(pd_fimo,
                            aes(y = score, x = atBreakSite)) + 
  geom_jitter(aes(color = teloSupport), width = 0.4, size = 0.8) +
  geom_text(data = tb_n, aes(y = 45, label = n), size = 3,
            angle = 90) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 50),
                     breaks = scales::pretty_breaks(5)) +
  ylab("Motif match score") +
  scale_color_manual(name = "Telomere\nsupport",
                     breaks = c("TRUE", "FALSE"),
                     labels = c("Yes", "No"),
                     values = c("#FF7F0E", "#1F77B4")) +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        # text = element_text(family = "Arial"),
        legend.position = "bottom")

# motif_specificity
#### Combine into single figure ####
# Align all images vertically for the top half of the figure
top_half <- align_plots(main_motif_plot, motif_specificity,
                        align = 'v', 
                        axis = 'l')
names(top_half) <- letters[1:2]

fig4 <- plot_grid(top_half$a, top_half$b,
                  labels = LETTERS[1:2],
                  label_size = 11,
                  rel_heights = c(1, 1.5),
                  ncol = 1)

ggsave("report/figures/figure_motif.pdf", fig4,
       units = "in", 
       width = 8.5*1.0*0.8,
       height = 11*0.8*0.8,
       dpi = 400
)


### Export table ####
filter(pd_fimo, score > 23) %>%
  arrange(atBreakSite, seqnames, start) %>%
  select(seqnames, start, end, score, atBreakSite) %>%
  write_tsv("report/tables/high_scoring_motif_positions.tsv")




