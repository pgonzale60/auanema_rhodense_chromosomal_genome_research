library(GenomicRanges)
library(rtracklayer)
library(scales)
library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(patchwork)

# In this figure I want to depict what's in the genome
# Both in eliminated and retained regions.
# Originally I planned to show the occurrence the distribution
# of genetic diversity but I decided to remove it
# because it closely reflects the distribution of predicted genes
# likely due to unmappability in the repetitive regions.
# Hence panel A will combine both gene densitty and tandem repeats


# Create a custom color palette with 11 distinguishable colors
custom_palette <- brewer.pal(12, "Paired")
three_custom_palette <- brewer.pal(3, "Set1")
six_custom_palette <- brewer.pal(6, "Set1")
thirteen_custom_palette <- c(custom_palette, "#a1a1a1")


### plotting Function ####
genom_dens_feat_plot <- function(dimi_positions, elim_regs,
                                 plData, window_size,
                                 ytext = "Density", maxy = 1,
                                 xtext = "Physical position (Mb)",
                                 mpalette = NULL, 
                                 pcategories = NULL,
                                 legnd_title = NULL){
  comb_data <- bind_rows(select(dimi_positions,
                                seqnames = Sequence,
                                diminution_pos = POS),
                         plData) %>%
    bind_rows(select(elim_regs,
                     seqnames = Sequence,
                     estarts = start,
                     eend = end)) %>%
    mutate(seqnames = paste0("chr", sub("SUPER_", "", seqnames)),
           binned_cov = ifelse(binned_cov > maxy, maxy, binned_cov))
  
  if(is.null(pcategories)){
    baseg <- ggplot(comb_data,
                    aes(x=start/1e6, y=binned_cov, diminution_pos)) +
      theme_bw()
  } else {
    baseg <- ggplot(comb_data,
                    aes(x=start/1e6, y=binned_cov, diminution_pos,
                        fill = .data[[pcategories]]))
  }
  if(!is.null(mpalette) && !is.null(pcategories) &&
     !is.null(legnd_title)){
    baseg <- baseg +
      scale_fill_manual(name = legnd_title,
                        values = mpalette) +
      theme_bw() +
      theme(legend.position = "bottom")
  }
  
  outplot <- baseg + 
    facet_grid(seqnames ~ .) +
    geom_rect(aes(xmin = estarts/1e6, xmax = eend/1e6), 
              ymin = 0, ymax = Inf, alpha = 0.5) +
    geom_vline(aes(xintercept = diminution_pos/1e6), linetype="dotted",
               color = "black", linewidth = 0.5) +
    geom_bar(position="stack", stat="identity",
             width = window_size/1e6, alpha = 0.8) +
    scale_x_continuous(n.breaks = 6, labels = scales::comma) +
    scale_y_continuous(breaks = scales::pretty_breaks(3), limits = c(0, maxy)) +
    ylab(ytext) +
    xlab(xtext) 
}


### Load data ####
# Sequence sizes
fai_files <- list.files("analyses/genome_features/sequence_sizes/",
                        full.names = T, pattern = ".primary.fa.gz.fai")
names(fai_files) <- make.names(sub(".+//(.+).primary.fa.gz.fai", "\\1", fai_files))
seq_sizes <- map_df(fai_files, read_tsv,
                    col_names = c("Sequence", "size", "L1", "L2", "L3"),
                    col_types = c("ciiii"),
                    .id = "assembly") %>%
  select(assembly, Sequence, size) %>%
  filter(!grepl("unloc|MT|scaff", Sequence)) %>%
  select(-assembly)

# GRS files
GRS_files <- list.files("analyses/genome_features/elim_coords/",
                        full.names = T, pattern = "nx.*.GRS.bed$")
names(GRS_files) <- make.names(sub(".+//(.+).GRS.bed", "\\1", GRS_files))
GRS <- map_df(GRS_files, read_tsv,
              col_names = c("Sequence", "start", "end"),
              col_types = c("cii"),
              .id = "assembly") %>%
  select(-assembly) #%>%

lGRS <- pivot_longer(GRS, !Sequence, names_to = "varia", values_to = "POS")

# Protein files
protfiles <- list.files("analyses/genes/GTFs/",
                        full.names = T, pattern = ".gz")
names(protfiles) <- make.names(sub(".+//(.+).gff.gz", "\\1", protfiles))
genes_df <- map_df(protfiles, read_tsv,
                   comment = "#", quote = "\"",
                   col_names = c("Sequence", "source", "type", "start", "end", "score",
                                 "strand", "phase", "ID"),
                   .id = "assembly") %>%
  filter(!grepl("unloc|MT|scaff", Sequence),
         type == "gene"
         # !type %in% c("intron", "start_codon", "stop_codon")
         )  %>%
  mutate(seqnames = Sequence) %>%
  select(-assembly, -Sequence)


# Genome-wide diversity
# rad <- read_tsv("analyses/genome_features/RAD_map/populations.snps.vcf.gz",
#                 comment = "#") %>%
#   select(CHROM, POS)  %>%
#   filter(!grepl("unloc|MT|scaff", CHROM))


# RNAmmer files
rrnafiles <- list.files("analyses/genome_features/repeats/rnammer/",
                        full.names = T, pattern = "gff3.gz")
names(rrnafiles) <- make.names(sub(".+//(.+).rnammer.gff3.gz", "\\1", rrnafiles))
rrna_df <- map_df(rrnafiles, read_tsv,
                  .id = "assembly", col_names = c("Sequence", "source", "type", "start", "end", "score",
                                                  "strand", "phase", "ID")) %>%
  mutate(seqnames = Sequence) %>%
  filter(!grepl("MT|scaf|unloc", Sequence)) %>%
  select(-assembly, -Sequence, -phase)

# tRNAscan files
trnafiles <- list.files("analyses/genome_features/repeats/tRNAscan/",
                        full.names = T, pattern = "gff.gz")
names(trnafiles) <- make.names(sub(".+//(.+).trnas.gff.gz", "\\1", trnafiles))
trna_df <- map_df(trnafiles, read_tsv,
                  .id = "assembly", col_names = c("Sequence", "source", "type", "start", "end", "score",
                                                  "strand", "phase", "ID")) %>% 
  mutate(seqnames = Sequence) %>%
  filter(!grepl("MT|scaf|unloc", Sequence)) %>%
  select(-assembly, -Sequence, -phase)


# infernal files
miscrnafiles <- list.files("analyses/genome_features/repeats/infernal/",
                           full.names = T, pattern = "gff.gz")
names(miscrnafiles) <- make.names(sub(".+//(.+).infernal.gff.gz", "\\1", miscrnafiles))
miscrna_df <- map_df(miscrnafiles, read_tsv,
                     .id = "assembly", col_names = c("Sequence", "source", "type", "start", "end", "score",
                                                     "strand", "phase", "ID")) %>% 
  filter(!grepl("rRNA|tRNA|Protozoa", type),
         !type %in% c("K_chan_RES", "RNase_MRP", "Fluoride", "GlsR7", "RAGATH-21", "snosnR60_Z15")) %>%
  mutate(seqnames = Sequence) %>%
  filter(!grepl("MT|scaf|unloc", Sequence))

# Tandem repeats
# include the rRNA repeat
rrna_clust <- filter(rrna_df, grepl("18s|28s", ID)) %>%
  mutate(qaccver = "cluster33", saccver = seqnames,
         pident = 100, length = end - start,
         qlen = 100, slen = 100, qstart = 0,
         qend = 100, sstart = start, send = end,
         evalue = 1e-4, bitscore = 100, qcovs = 100) %>%
  select(-c(source:seqnames))
aua_rep_blast <- read_tsv("analyses/genome_features/repeats/TRF/arhod_reps_clustrd_plus_chr5C.tsv.gz", # arhod_reps_clustrd_plus_chr5C.tsv
                          col_names = c("qaccver", "saccver", "pident",
                                        "length", "qlen", "slen", 
                                        "qstart", "qend", "sstart", 
                                        "send", "evalue", "bitscore", 
                                        "qcovs") )%>%
  filter(!grepl("MT|scaf|unloc", saccver),
         qaccver != "cluster17") %>% # cluster17 is redundant with cluster13
  rbind(rrna_clust) %>% 
  mutate(qaccver = sub("cluster", "", qaccver)) 
  
  

### Preprocess tibbles ####

elim_seqs <- filter(seq_sizes, grepl("loc|scaffold|MT", Sequence)) %>%
  mutate(start = 1, end = size-1) %>%
  select(-size) %>%
  rbind(GRS)


#### Tandem repeats #####
large_TR <- filter(aua_rep_blast, !qaccver %in%
                     c("13", "18", # these are dispersed across the genome
                       "6")) %>% # cluster6 is semi compact, but I prefer more compact
  group_by(qaccver) %>%
  summarise(tot_hit_span = sum(length)) %>%
  ungroup() %>%
  slice_max(tot_hit_span, n = 11) %>%
  pull(qaccver) %>% c("other")

rep_for_plot <- mutate(aua_rep_blast, chrom=saccver,
         start = ifelse(sstart < send, sstart, send),
         end = ifelse(sstart > send, sstart, send),
         rep_id = 
           # qaccver
         ifelse(qaccver %in% large_TR,
                         qaccver, "other")
         ) %>%
  select(chrom, start, end, rep_id)


#### RAD
# rad_fg <- mutate(rad, chrom=CHROM, 
#                  start = POS, end = POS+1) %>%
#   select(chrom, start, end)

#### ncRNAs ####
# RNAmmer 
rrna_df <- mutate(rrna_df, rRNA_type = sub("ID=([0-9]+)s.+", "\\1S", ID),
                  gtypeB = ifelse(rRNA_type == "8S", "5S",
                                  rRNA_type),
                  gtypeA = "rRNA")
abund_rrna <- count(rrna_df, gtypeB) %>% 
  arrange(desc(n)) %>%
  pull(gtypeB)

# tRNAscan 
trna_df <- mutate(trna_df, tRNA_type = sub("ID=tRNA-([^0-9]+)[0-9].+", "\\1", ID),
                  anticodon = sub(".*anticodon=(...).*", "\\1", ID),
                  gtypeA = "tRNA")
abund_trna <- count(trna_df, tRNA_type) %>% 
  filter(n > 800) %>% arrange(desc(n)) %>%
  pull(tRNA_type) %>% c("other")
trna_df <- mutate(trna_df, 
                  gtypeB = ifelse(tRNA_type %in% abund_trna,
                                  tRNA_type, "other"))

# infernal files
miscrna_df <- mutate(miscrna_df, evalue = sub("evalue=([^;]+).+", "\\1", ID),
                     description = sub(".+desc=([^;]+)", "\\1", ID),
                     misc_type = sub("([^_]+)_.*", "\\1", description),
                     gtypeA = "misc")
abund_miscrna <- count(miscrna_df, misc_type) %>% 
  filter(n > 116) %>% arrange(desc(n)) %>%
  pull(misc_type) %>% c("other")
miscrna_df <- mutate(miscrna_df, gtypeB = ifelse(misc_type %in% abund_miscrna,
                                     misc_type, "other"))




### Genomic ranges ####
gnm_gr <- Seqinfo(seqnames = seq_sizes$Sequence, seqlengths=seq_sizes$size,
                  isCircular = rep(F, nrow(seq_sizes)), genome = "Auanema")
# GRS regions
GRS_gr <- GRanges(seqnames = Rle(elim_seqs$Sequence), 
                  ranges = IRanges(start = elim_seqs$start,
                                   end = elim_seqs$end), 
                  strand = "*")
# genes_gr <- makeGRangesFromDataFrame(genes_df, keep.extra.columns=T,
#                                      seqinfo=gnm_gr)

gtrf_gr <- select(genes_df, chrom = seqnames,
                  start, end, rep_id = type) %>%
  rbind(rep_for_plot) %>%
  makeGRangesFromDataFrame(seqinfo = gnm_gr, keep.extra.columns = T)
# rad_gr <- makeGRangesFromDataFrame(rad_fg, seqinfo = gnm_gr, keep.extra.columns = T)
# trf_gr <- makeGRangesFromDataFrame(rep_for_plot, seqinfo = gnm_gr, keep.extra.columns = T)
# RNA granges
rrna_gr <- makeGRangesFromDataFrame(rrna_df, keep.extra.columns=T, seqinfo=gnm_gr)
trna_gr <- makeGRangesFromDataFrame(trna_df, keep.extra.columns=T, seqinfo=gnm_gr)
miscrna_gr <- makeGRangesFromDataFrame(miscrna_df, keep.extra.columns=T, seqinfo=gnm_gr)
all_nc <- c(rrna_gr, trna_gr, miscrna_gr)


### Coverages ####
window_size <- 100000
aua_1kb <- tileGenome(gnm_gr, tilewidth = window_size, cut.last.tile.in.chrom = T)

#### RAD
# rad_cov <- GenomicRanges::coverage(rad_gr)
# rad_abundance <- GenomicRanges::binnedAverage(aua_1kb, rad_cov, "binned_cov") %>%
#   as_tibble()

#### Genes and Tandem ####
# gene_cov <- GenomicRanges::coverage(genes_gr)
# gene_abundance <- GenomicRanges::binnedAverage(aua_1kb, gene_cov, "binned_cov") %>%
#   as_tibble()
#### Tandem 
# tr_span <- tibble(clust_id = character(), span = integer())
# for (clust_id in unique(rep_for_plot$rep_id)) {
#   tr_span <- tibble(clust_id = clust_id,
#          span = sum(width(GenomicRanges::reduce(trf_gr[trf_gr$rep_id == clust_id])))) %>%
#     rbind(tr_span)
# }
# 
# trf_gr[trf_gr[trf_gr$rep_id != "cluster13"] %over% trf_gr[trf_gr$rep_id == "cluster13"]] %>%
#   as_tibble() %>%
#   count(rep_id)
# View(tr_span)
# filter(rep_for_plot, rep_id %in% c("cluster13", "cluster16")) %>% View

seq_factors <- c(large_TR, "gene")
gtr_abundance <- tibble(seqnames = character(), start = integer(), end = integer(),
                       width = integer(), strand = factor(), binned_cov = double(),
                       clust_id = character())
gtcumRang <- gtrf_gr[0]
for (clust_id in seq_factors) {
  clust_gr <- GenomicRanges::reduce(gtrf_gr[gtrf_gr$rep_id == clust_id])
  uclust_gr <- unlist(GenomicRanges::subtract(clust_gr, gtcumRang))
  gtcumRang <- c(gtcumRang, uclust_gr)
  gtr_cov <- GenomicRanges::coverage(uclust_gr)
  gtr_abundance <- GenomicRanges::binnedAverage(aua_1kb, gtr_cov, "binned_cov") %>%
    as_tibble() %>%
    mutate(clust_id = clust_id) %>%
    rbind(gtr_abundance)
}

gtr_abundance <- mutate(gtr_abundance, 
                        clust_id = factor(clust_id,
                                          levels = seq_factors))

# 
# seq_factors <- large_TR
# seq_factors <- "cluster33"
# tr_abundance <- tibble(seqnames = character(), start = integer(), end = integer(),
#                           width = integer(), strand = factor(), binned_cov = double(),
#                           clust_id = character())
# for (clust_id in seq_factors) {
#   tr_cov <- GenomicRanges::coverage(GenomicRanges::reduce(trf_gr[trf_gr$rep_id == clust_id]))
#   tr_abundance <- GenomicRanges::binnedAverage(aua_1kb, tr_cov, "binned_cov") %>%
#     as_tibble() %>%
#     mutate(clust_id = clust_id) %>%
#     rbind(tr_abundance)
# }
# 
# tr_abundance <- mutate(tr_abundance, clust_id = factor(clust_id,
#                                                        levels = seq_factors))



#### ncRNA ####
seq_factors <- rev(c("tRNA", "rRNA", "misc"))
var_to_check <- "gtypeA"

ncRNA_abundance <- tibble(seqnames = character(), start = integer(), end = integer(),
                          width = integer(), strand = factor(), binned_cov = double(),
                          ncRNA = character())
gtcumRang <- all_nc[0]
for (clust_id in seq_factors) {
  clust_gr <- GenomicRanges::reduce(all_nc[mcols(all_nc)[, var_to_check] == clust_id])
  uclust_gr <- unlist(GenomicRanges::subtract(clust_gr, gtcumRang))
  gtcumRang <- c(gtcumRang, uclust_gr)
  ncRNA_cov <- GenomicRanges::coverage(uclust_gr)
  ncRNA_abundance <- GenomicRanges::binnedAverage(aua_1kb, ncRNA_cov, "binned_cov") %>%
    as_tibble() %>%
    mutate(ncRNA = clust_id) %>%
    rbind(ncRNA_abundance)
}

ncRNA_abundance <- mutate(ncRNA_abundance, 
                          ncRNA = factor(ncRNA,
                                         levels = seq_factors))


#### rRNA ####
seq_factors <- abund_rrna
var_to_check <- "gtypeB"

rRNA_abundance <- tibble(seqnames = character(), start = integer(), end = integer(),
                          width = integer(), strand = factor(), binned_cov = double(),
                          rRNA = character())
gtcumRang <- all_nc[0]
for (clust_id in seq_factors) {
  clust_gr <- GenomicRanges::reduce(all_nc[mcols(all_nc)[, var_to_check] == clust_id])
  uclust_gr <- unlist(GenomicRanges::subtract(clust_gr, gtcumRang))
  gtcumRang <- c(gtcumRang, uclust_gr)
  rRNA_cov <- GenomicRanges::coverage(uclust_gr)
  rRNA_abundance <- GenomicRanges::binnedAverage(aua_1kb, rRNA_cov, "binned_cov") %>%
    as_tibble() %>%
    mutate(rRNA = clust_id) %>%
    rbind(rRNA_abundance)
}

rRNA_abundance <- mutate(rRNA_abundance, 
                          rRNA = factor(rRNA,
                                        levels = seq_factors))

#### tRNA ####
seq_factors <- abund_trna
var_to_check <- "gtypeB"

tRNA_abundance <- tibble(seqnames = character(), start = integer(), end = integer(),
                         width = integer(), strand = factor(), binned_cov = double(),
                         tRNA = character())
gtcumRang <- all_nc[0]
for (clust_id in seq_factors) {
  clust_gr <- GenomicRanges::reduce(all_nc[mcols(all_nc)[, var_to_check] == clust_id])
  uclust_gr <- unlist(GenomicRanges::subtract(clust_gr, gtcumRang))
  gtcumRang <- c(gtcumRang, uclust_gr)
  tRNA_cov <- GenomicRanges::coverage(uclust_gr)
  tRNA_abundance <- GenomicRanges::binnedAverage(aua_1kb, tRNA_cov, "binned_cov") %>%
    as_tibble() %>%
    mutate(tRNA = clust_id) %>%
    rbind(tRNA_abundance)
}

tRNA_abundance <- mutate(tRNA_abundance, 
                         tRNA = factor(tRNA,
                                       levels = seq_factors))


#### miscRNA ####
seq_factors <- abund_miscrna
var_to_check <- "gtypeB"

miscRNA_abundance <- tibble(seqnames = character(), start = integer(), end = integer(),
                         width = integer(), strand = factor(), binned_cov = double(),
                         miscRNA = character())
gtcumRang <- all_nc[0]
for (clust_id in seq_factors) {
  clust_gr <- GenomicRanges::reduce(all_nc[mcols(all_nc)[, var_to_check] == clust_id])
  uclust_gr <- unlist(GenomicRanges::subtract(clust_gr, gtcumRang))
  gtcumRang <- c(gtcumRang, uclust_gr)
  miscRNA_cov <- GenomicRanges::coverage(uclust_gr)
  miscRNA_abundance <- GenomicRanges::binnedAverage(aua_1kb, miscRNA_cov, "binned_cov") %>%
    as_tibble() %>%
    mutate(miscRNA = clust_id) %>%
    rbind(miscRNA_abundance)
}

miscRNA_abundance <- mutate(miscRNA_abundance, 
                            miscRNA = factor(miscRNA,
                                       levels = seq_factors))



### Plots ####
#### RAD 
# rad_plt <- genom_dens_feat_plot(lGRS,
#                      rad_abundance, window_size,
#                      ytext = "SNP density per 100 Kb window",
#                      maxy = 0.002,
#                      xtext = "Physical position (Mb)"
# )

# gene_plt <- genom_dens_feat_plot(lGRS,
#                                  gene_abundance, window_size,
#                                 ytext = "Gene density per 100 Kb window",
#                                 maxy = 0.9,
#                                 xtext = "Physical position (Mb)"
# )

#### Gene and tandem ####

gtrf_plt <- genom_dens_feat_plot(lGRS, GRS,
                                 gtr_abundance, window_size,
                                ytext = "Tandem repeat density per 100 Kb window",
                                maxy = 1,
                                xtext = "Physical position (Mb)",
                                mpalette = thirteen_custom_palette,
                                pcategories = "clust_id",
                                legnd_title = "Tandem\nrepeat"
)


  
#### Tandem ####
# trf_plt <- genom_dens_feat_plot(lGRS,
#                                 tr_abundance, window_size,
#                                 ytext = "Tandem repeat density per 100 Kb window",
#                                 maxy = 1,
#                                 xtext = "Physical position (Mb)",
#                                 mpalette = custom_palette,
#                                 pcategories = "clust_id",
#                                 legnd_title = "Tandem\nrepeat"
# )


#### rRNA ####
rRNA_plt <- genom_dens_feat_plot(lGRS, GRS,
                                rRNA_abundance, window_size,
                                ytext = "rRNA density per 100 Kb window",
                                maxy = 1,
                                xtext = "Physical position (Mb)",
                                mpalette = three_custom_palette,
                                pcategories = "rRNA",
                                legnd_title = "rRNA"
)


#### tRNA ####
tRNA_plt <- genom_dens_feat_plot(lGRS, GRS,
                                 tRNA_abundance, window_size,
                                 ytext = "tRNA density per 100 Kb window",
                                 maxy = 0.2,
                                 xtext = "Physical position (Mb)",
                                 mpalette = custom_palette,
                                 pcategories = "tRNA",
                                 legnd_title = "tRNA"
)


#### misc ####
misc_plt <- genom_dens_feat_plot(lGRS, GRS,
                                 miscRNA_abundance, window_size,
                                 ytext = "tRNA density per 100 Kb window",
                                 maxy = 0.2,
                                 xtext = "Physical position (Mb)",
                                 mpalette = six_custom_palette,
                                 pcategories = "miscRNA",
                                 legnd_title = "miscRNA"
)



#### Combine into single figure ####

# rad_plt, gene_plt, trf_plt, rRNA_plt, tRNA_plt, misc_plt
# rad_plt, gene_plt, trf_plt, tRNA_plt,
# Align all images vertically for the top half of the figure
all_plt <- align_plots(gtrf_plt, tRNA_plt,
                       rRNA_plt, misc_plt,
                       align = 'v', 
                       axis = 'l')

names(all_plt) <- letters[1:4]

top_row <- plot_grid(all_plt$a, all_plt$b,
                     labels = LETTERS[1:2],
                     nrow = 1)
mid_row <- plot_grid(all_plt$c, all_plt$d,
                     labels = LETTERS[3:4],
                     nrow = 1)
# bottom_row <- plot_grid(all_plt$e, all_plt$f,
#                         labels = LETTERS[5:6],
#                         nrow = 1)

fig4.4 <- plot_grid(top_row, mid_row, # bottom_row, 
                  rel_heights = c(1.2, 1),
                  nrow = 2)


ggsave("report/figures/figure4.4.pdf", fig4.4,
       units = "in", 
       width = 8.5*1.1,
       height = 11*0.8
)

## Tables #####
elt_rrna <- rrna_gr[rrna_gr %over% GRS_gr] %>%
  as_tibble() %>%
  count(gtypeB) %>%
  mutate(type = "rrna")

elt_trna <- trna_gr[trna_gr %over% GRS_gr] %>%
  as_tibble() %>%
  count(gtypeB) %>%
  mutate(type = "trna")

elt_misc <- miscrna_gr[miscrna_gr %over% GRS_gr] %>%
  as_tibble() %>%
  count(gtypeB) %>%
  mutate(type = "miscrna")

# tot
tot_rrna <- rrna_gr %>%
  as_tibble() %>%
  count(gtypeB) %>%
  mutate(type = "rrna")

tot_trna <- trna_gr %>%
  as_tibble() %>%
  count(gtypeB) %>%
  mutate(type = "trna")

tot_misc <- miscrna_gr %>%
  as_tibble() %>%
  count(gtypeB) %>%
  mutate(type = "miscrna")

part <- bind_rows(elt_rrna, elt_trna, elt_misc) %>%
  mutate(elim = "part")

tot <- bind_rows(tot_rrna, tot_trna, tot_misc) %>%
     mutate(elim = "total")

bind_rows(part, tot) %>%
  group_by(type) %>%
  pivot_wider(values_from = n,
              names_from = c(elim),
              values_fill = 0) %>%
  ungroup() %>%
  mutate(perc_elim = round(100 * (part/total), 1),
         core = total - part) %>%
  arrange(desc(perc_elim)) %>%
  select(gtypeB, type, part, core, perc_elim) %>%
  write_tsv("report/tables/eliminated_ncRNAs.tsv")

count(trna_df, tRNA_type) %>%
  arrange(desc(n))

#### Alex SINEs test ####
atrna_gr <- trna_gr
atrna_gr$eliminated <- atrna_gr %over% GRS_gr
atrna_df <- as_tibble(atrna_gr)
# count(anticodon, eliminated) %>%
atrna_smr <- group_by(atrna_df, eliminated, anticodon) %>%
  summarise(tRNA_type = unique(tRNA_type),
            minscore = min(score),
            maxscore = max(score),
            meanscore = mean(score),
            n = n()) %>%
  arrange(desc(n)) #%>%
View(atrna_smr)
View(atrna_df)

abund_antic <- filter(atrna_smr, n > 10) %>%
  pull(anticodon) %>% unique

# filter(atrna_df, anticodon %in% abund_antic) %>%
#   ggplot(aes(x = anticodon, y = score, fill = eliminated)) +
#   geom_boxplot(alpha = 0.6, outlier.shape = NA) +  # Boxplot, remove outliers as these will be covered by the jitter plot
#   geom_jitter(aes(color = eliminated), alpha = 0.4, width = 0.2) +  # Points as jitter, set width as required
#   scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "red")) +  # Set fill colors as required
#   scale_color_manual(values = c("TRUE" = "lightblue", "FALSE" = "pink")) +  # Set point colors as required
#   theme(axis.text.x = element_text(angle = 90)) +  # Rotate x-axis labels, might be needed if labels are long
#   labs(title = "Boxplot of Scores by Anticodon", x = "Anticodon", y = "Score")  # Add labels


# First, we define the codon table
# codon_table <- c("AAA"="Lys", "AAC"="Asn", "AAG"="Lys", "AAU"="Asn", "ACA"="Thr", 
#                  "ACC"="Thr", "ACG"="Thr", "ACU"="Thr", "AGA"="Arg", "AGC"="Ser", 
#                  "AGG"="Arg", "AGU"="Ser", "AUA"="Ile", "AUC"="Ile", "AUG"="Met", 
#                  "AUU"="Ile", "CAA"="Gln", "CAC"="His", "CAG"="Gln", "CAU"="His", 
#                  "CCA"="Pro", "CCC"="Pro", "CCG"="Pro", "CCU"="Pro", "CGA"="Arg", 
#                  "CGC"="Arg", "CGG"="Arg", "CGU"="Arg", "CUA"="Leu", "CUC"="Leu", 
#                  "CUG"="Leu", "CUU"="Leu", "GAA"="Glu", "GAC"="Asp", "GAG"="Glu", 
#                  "GAU"="Asp", "GCA"="Ala", "GCC"="Ala", "GCG"="Ala", "GCU"="Ala", 
#                  "GGA"="Gly", "GGC"="Gly", "GGG"="Gly", "GGU"="Gly", "GUA"="Val", 
#                  "GUC"="Val", "GUG"="Val", "GUU"="Val", "UAA"="Stop", "UAC"="Tyr", 
#                  "UAG"="Stop", "UAU"="Tyr", "UCA"="Ser", "UCC"="Ser", "UCG"="Ser", 
#                  "UCU"="Ser", "UGA"="Stop", "UGC"="Cys", "UGG"="Trp", "UGU"="Cys", 
#                  "UUA"="Leu", "UUC"="Phe", "UUG"="Leu", "UUU"="Phe")

# Now let's say you have a list of anticodo

# Here's a function to reverse complement a DNA sequence
# reverse_complement <- function(sequence) {
#   return(paste(rev(chartr("ATGC", "UACG", strsplit(toupper(sequence), "")[[1]])), collapse = ""))
# }


# You can reverse complement these and find the corresponding amino acids like this:
# Create summary table
# amino_top20row <- select(atrna_smr, anticodon, tRNA_type, eliminated, n) %>%
#   mutate(eliminated = ifelse(eliminated, "GRS", "core")) %>%
#   pivot_wider(id_cols = c("anticodon", "tRNA_type"), names_from = "eliminated",
#               values_from = "n", values_fill = 0) %>%
#   arrange(desc(GRS+core)) %>%
#   # mutate(codon = sapply(anticodon, reverse_complement),
#   #        aminoacid = codon_table[codon]) %>%
#   # select(anticodon, aminoacid, GRS, core) %>%
#   head(20) %>% arrange(tRNA_type)

amino_top20row <- select(atrna_smr, anticodon, tRNA_type, eliminated, n) %>%
  mutate(eliminated = ifelse(eliminated, "eliminated", "core")) %>%
  pivot_wider(id_cols = c("anticodon", "tRNA_type"), names_from = "eliminated",
              values_from = "n", values_fill = 0) %>%
  filter(eliminated + core > 10) %>%
  group_by(tRNA_type) %>%
  mutate(ntype = sum(eliminated) + sum(core)) %>%
  ungroup() %>%
  arrange(desc(ntype), anticodon) %>%
  mutate(aaa = paste(anticodon, tRNA_type, sep = "_"),
         aaa = factor(aaa,
                            levels = rev(unique(aaa))))

freq_cod_plot <-
  select(amino_top20row, -ntype) %>%
pivot_longer(cols = c("eliminated", "core")) %>%
  filter(value > 0) %>%
  ggplot(aes(y = aaa, x = value, color = name)) +
  geom_point() +
  labs(x = "", y = "Anticodon", color = "Localization") +
  theme_minimal() +
  theme(legend.position = c(0.95, 0.05),
        legend.justification = c(1, 0),
        legend.box.background = element_rect(fill = "white", color = "black"),
        legend.margin = margin(0.1, 0.1, 0.0, 0.1, "cm"),
        plot.margin = margin(0.0, 0.0, 0.0, 0.0, "cm")) 
  


# filter(atrna_df, anticodon %in% amino_top20row$anticodon) %>%
trna_score_plot <- mutate(atrna_df,
       eliminated = ifelse(eliminated, "eliminated", "core")) %>%
ggplot(aes(x = eliminated, y = score)) +
  coord_cartesian(ylim = c(0, 90)) +
  geom_jitter(alpha = 0.4, width = 0.2) +  # Points as jitter, set width as required
  geom_boxplot(alpha = 0.6, outlier.shape = NA, color = "blue") +  # Boxplot, remove outliers as these will be covered by the jitter plot
  theme_minimal() +
  # theme(axis.text.x = element_text(angle = 90)) +  # Rotate x-axis labels, might be needed if labels are long
  labs(x = "", y = "Cove score")  # Add labels

# Combine the plots into a single composite figure
composite_figure <- trna_score_plot + freq_cod_plot +
  plot_annotation(tag_levels = 'A')

# Display the composite figure
composite_figure

ggsave(plot = composite_figure, filename = "report/figures/tRNA_anticodons.pdf",
       width = 8, height = 8)


