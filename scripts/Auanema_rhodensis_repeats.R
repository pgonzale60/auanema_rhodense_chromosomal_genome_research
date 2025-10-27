library(Gviz)


library(GenomicRanges)   # Load the GenomicRanges library for genomic range operations
library(tidyverse)       # Load the tidyverse library for data manipulation

# Read a TSV file containing data on AuaRhod1.1 TRF (tandem repeat finder) results
aua_reps <- read_tsv("analyses/genome_features/repeats/TRF/windowed100K_nxAuaRhod1_1.trf.tsv.gz",
                     col_names = c("chrom", "start", "end", "unit_size",
                                   "num_copies", "perc_ident", "perc_indels",
                                   "score", "entropy", "sequence")) %>%
  mutate(prestart = as.integer(sub(".+:(.+)-.+", "\\1", chrom)),
         start = start + prestart,
         end = end + prestart,
         chrom = sub(":.+", "", chrom))


# Filter AuaRhod1.1 TRF data for rows where the score is greater than 80000
# Group the filtered data by chromosome and select the row with the minimum unit_size for each group
# Ungroup the data and create a new column seq_for_clust by concatenating the sequence with itself
# Select the chrom and seq_for_clust columns and write them to a TSV file

aua_rep_blast <- read_tsv("analyses/genome_features/repeats/TRF/arhod_reps_clustrd.tsv.gz", # arhod_reps_clustrd_plus_chr5C.tsv
                          col_names = c("qaccver", "saccver", "pident",
                                        "length", "qlen", "slen", 
                                        "qstart", "qend", "sstart", 
                                        "send", "evalue", "bitscore", 
                                        "qcovs"))


coord_aua_rep_blast <- group_by(aua_rep_blast, qaccver, saccver) %>%
  arrange(sstart) %>%
  mutate(block = cumsum(c(1, diff(sstart) > 10000))) %>%
  ungroup() %>%
  group_by(qaccver, saccver, block) %>%
  summarise(
    usstart = min(c(sstart, send)),
    usend = max(c(sstart, send)),
    avg_ident = round(mean(pident), 2),
    span = usend-usstart,
    n_hits = n(),
    slen = unique(slen),
    .groups = "drop") %>%
  arrange(qaccver, saccver, usstart) 

frac_per_aua_chr <- group_by(coord_aua_rep_blast, qaccver, saccver) %>%
  summarise(
    len_per_chrom = sum(usend-usstart),
    slen = unique(slen),
    frac_per_chrom = round((len_per_chrom/slen)*100, 2),
    .groups = "drop") 

big_clusters <- filter(frac_per_aua_chr, !grepl("unloc|scaff", saccver),
       frac_per_chrom > 1) %>% 
  pull(qaccver) %>% unique()

rep_for_plot <- filter(aua_rep_blast, qaccver %in% big_clusters,
       !grepl("unloc|scaff", saccver))
  



library(rtracklayer)
library(tidyverse)
library(scales)

### Core tibble creation ####
# Sequence sizes
fai_files <- list.files("analyses/genome_features/sequence_sizes/",
                        full.names = T, pattern = ".primary.fa.gz.fai")
names(fai_files) <- make.names(sub(".+//(.+).primary.fa.gz.fai", "\\1", fai_files))
seq_sizes <- map_df(fai_files, read_tsv,
                    col_names = c("Sequence", "size", "L1", "L2", "L3"),
                    col_types = c("ciiii"),
                    .id = "assembly") %>%
  select(assembly, Sequence, size) %>%
  filter(assembly == "nxAuaRhod1.1",
         Sequence %in% unique(rep_for_plot$saccver)) %>%
  mutate(multispecies_sequence = paste0(assembly, "_", Sequence))

# Break sites
break_files <- list.files("analyses/curated_GRS_coords/",
                          full.names = T, pattern = "nxAuaRhod1.1.chr_diminutions_sites.tsv")
names(break_files) <- make.names(sub(".+//(.+).chr_diminutions_sites.tsv", "\\1", break_files))
break_sites <- map_df(break_files, read_tsv,
                      col_names = T,
                      .id = "assembly") %>%
  filter(#feature == "telomere_seq_split",
         assembly == "nxAuaRhod1.1") %>%
  mutate(multispecies_sequence = paste0(assembly, "_", chr))



#### Genomic ranges
gnm_gr <- Seqinfo(seqnames = seq_sizes$multispecies_sequence, seqlengths=seq_sizes$size,
                  isCircular = rep(F, nrow(seq_sizes)), genome = "Auanema")

trf_gr <- mutate(rep_for_plot, chrom=paste0("nxAuaRhod1.1_", saccver),
         start = ifelse(sstart < send, sstart, send),
         end = ifelse(sstart > send, sstart, send),
         rep_id = qaccver) %>%
  select(chrom,
         start,
         end,
         rep_id) %>%
  makeGRangesFromDataFrame(seqinfo = gnm_gr, keep.extra.columns = T)

# What is it's distribution in the chr context?
window_size <- 100000
aua_1kb <- tileGenome(gnm_gr, tilewidth = window_size, cut.last.tile.in.chrom = T)

### Get coverage signal as Rle object

### Get average coverage in each bin
# (since the bins are 1-bp wide, this just keeps the original coverage value)

clust_abundance <- tibble(seqnames = character(), start = integer(), end = integer(),
       width = integer(), strand = factor(), binned_cov = double(),
       clust_id = character())
for (clust_id in unique(trf_gr$rep_id)) {
  gr1_cov <- GenomicRanges::coverage(GenomicRanges::reduce(trf_gr[trf_gr$rep_id == clust_id]))
  clust_abundance <- GenomicRanges::binnedAverage(aua_1kb, gr1_cov, "binned_cov") %>%
    as_tibble() %>%
    mutate(clust_id = clust_id) %>%
    rbind(clust_abundance)
  
}

trf_plt <- 
  bind_rows(dplyr::select(break_sites, seqnames = multispecies_sequence,
                          diminution_pos),
            clust_abundance) %>%
  mutate(chr = sub("[^_]+[_\\.]1_SUPER_(.+)", "chr\\1", seqnames),
         clust_id = sub("cluster", "", clust_id),
         clust_id = factor(clust_id, levels = order(clust_id))) %>%
  filter(!grepl("MT|unloc|scaff", seqnames)) %>%
  ggplot(aes(x=start, y=binned_cov, fill = clust_id, diminution_pos)) + 
  facet_grid(chr ~ .) +
  geom_bar(position="stack", stat="identity", width = window_size, alpha = 0.8) +
  geom_vline(aes(xintercept = diminution_pos), linetype="dotted",
             color = "black", size = 1) +
  scale_x_continuous(labels = label_number(scale_cut = cut_si("M"))) +
  scale_y_continuous(breaks = scales::pretty_breaks(3), limits = c(0, 1)) +
  # scale_fill_brewer(palette="Set1") +
  ylab("fraction per 100 Kb window") +
  xlab("Position along chromosome") +
  theme_bw()

ggsave("report/paper2_prelim_20220527/figures/exploratory/repeats/tandem/nxAuaRhod1.1_colored_by_tandem_repeat_many_reps_with_5.trf.pdf",
       trf_plt, width = 9, height = 6)








# Get GRS coordinates
### Core tibble creation ####

# TeloCoords
telocoord_files <- list.files("analyses/curated_GRS_coords/",
                              full.names = T, pattern = ".clippedTeloPos.tsv")
names(telocoord_files) <- make.names(sub(".+//(.+).clippedTeloPos.tsv", "\\1", telocoord_files))
telocoord <- map_df(telocoord_files, read_tsv,
                    col_names = c("Sequence", "diminution_pos", "telo_orient", "mapq"),
                    .id = "assembly") %>%
  filter(mapq > 30) %>%
  filter(!grepl("loc|MT", Sequence))%>%
  mutate(assembly = ifelse(assembly == "oscheius_tipulae.PRJNA644888",
                           "nxOscTipu1.1", assembly),
         multispecies_sequence = paste0(assembly, "_", Sequence)) %>%
  filter(assembly == "nxAuaRhod1.1")


#### Group and filter tibbles ####
# Group positions by exact same position
telocoord_per_pos <- mutate(telocoord, break_id = paste0(multispecies_sequence, "_", diminution_pos)) %>%
  group_by(break_id, telo_orient) %>%
  mutate(precise_support = n()) %>%
  slice(1) %>%
  ungroup()


#### GRanges creation ####
# Whole genome sizes
gnm_gr <- Seqinfo(seqnames = seq_sizes$multispecies_sequence, seqlengths=seq_sizes$size,
                  isCircular = rep(F, nrow(seq_sizes)), genome = "Auaunema")

telocoord_gr <- mutate(telocoord_per_pos, chrom=paste0(assembly, "_", Sequence),
                       end = diminution_pos) %>%
  select(chrom,
         start = diminution_pos,
         end,
         precise_support,
         break_id,
         telo_orient) %>%
  makeGRangesFromDataFrame(seqinfo = gnm_gr, keep.extra.columns = T)






#### Define GRS ####
# Semiautomatic approach
# export_seq_sizes_template <- filter(seq_sizes, multispecies_sequence %in% telocoord_grp$chrom) %>%
#   mutate(assembly = sub("([^_]+)_.+", "\\1", multispecies_sequence),
#          chr = sub("[^_]+_(.+)", "\\1", multispecies_sequence),
#          telo_orient = "R", precise_support = 0, joint_support = 0) %>%
#   select(assembly, chr, diminution_pos = size, telo_orient, precise_support, joint_support) 

export_seq_sizes_template <- filter(seq_sizes, multispecies_sequence %in% telocoord_per_pos$multispecies_sequence) %>%
  mutate(assembly = sub("([^_]+)_.+", "\\1", multispecies_sequence),
         chr = sub("[^_]+_(.+)", "\\1", multispecies_sequence),
         telo_orient = "R", precise_support = 0, joint_support = 0) %>%
  select(assembly, chr, diminution_pos = size, telo_orient, precise_support, joint_support) 

export_seq_sizes <- mutate(export_seq_sizes_template, diminution_pos = 1, telo_orient = "L") %>%
  bind_rows(export_seq_sizes_template) 


# Let's filter by at least 4 joint support

filter(telocoord_per_pos, precise_support > 100) %>%
  select(assembly, chr = Sequence, diminution_pos, telo_orient, precise_support) %>%
  bind_rows(select(export_seq_sizes, -joint_support)) %>%
  arrange(precise_support) %>%
  mutate(coord_id = paste0(assembly, chr, diminution_pos)) %>%
  filter(!duplicated(coord_id)) %>%
  select(-coord_id) %>%
  arrange(assembly, chr, diminution_pos) %>%
  write_tsv("~/Downloads/tmp.tsv")

# Group them by pattern of R and L

filt_telo_pos <- read_tsv("analyses/paper2/prelim_20220527/curated_GRS_coords/all_sp_teloclipped_based_20221007.tsv") %>%
  filter(for_coord != "no") %>%
  arrange(assembly, chr, diminution_pos)


grs_coords <- tibble(assembly = character(), chr = character(), start = double(), end = double())

for (i in 1:(nrow(filt_telo_pos)-1)) {
  cur_row <- filt_telo_pos[i,]
  next_row <- filt_telo_pos[i+1,]
  out_row <- grs_coords[0,]
  if(cur_row$chr != next_row$chr){
    next()
  }
  
  
  if(cur_row$telo_orient == next_row$telo_orient | (cur_row$telo_orient == "R" &
                                                    next_row$telo_orient == "L")){
    grs_coords <- tibble(assembly = cur_row$assembly,
                         chr = cur_row$chr, start = cur_row$diminution_pos,
                         end = next_row$diminution_pos) %>%
      bind_rows(grs_coords)
  }
  
}

grs_coords <- arrange(grs_coords, assembly, chr, start)
# write_tsv(grs_coords, "analyses/paper2/prelim_20220527/curated_GRS_coords/GRS_all_sp_teloclipped_based_20221007.tsv")


filter(break_sites, feature != "chr_end") %>%
  transmute(chr, start = diminution_pos-100, end = diminution_pos+100) %>%
  transmute(out=paste0(chr, ":", start, "-", end)) %>%
  write_tsv("~/Downloads/tmp.tsv", col_names = F)

