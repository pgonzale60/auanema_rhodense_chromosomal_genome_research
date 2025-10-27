library(GenomicRanges)
library(rtracklayer)
library(scales)
library(tidyverse)
library(RColorBrewer)


# Create a custom color palette with 11 distinguishable colors
custom_palette <- brewer.pal(12, "Paired")
three_custom_palette <- brewer.pal(3, "Set1")

### Load data ####
# Sequence sizes
fai_files <- list.files("analyses/genome_features/sequence_sizes/",
                        full.names = T, pattern = "nxAuaRhod1_1.primary.fa.gz.fai")
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
                        full.names = T, pattern = "nxAuaRhod.*.GRS.bed$")
names(GRS_files) <- make.names(sub(".+//(.+).GRS.bed", "\\1", GRS_files))
GRS <- map_df(GRS_files, read_tsv,
              col_names = c("Sequence", "start", "end"),
              col_types = c("cii"),
              .id = "assembly") %>%
  select(-assembly) #%>%

lGRS <- pivot_longer(GRS, !Sequence, names_to = "varia", values_to = "POS")



# Genome-wide diversity
rad <- read_tsv("analyses/genome_features/RAD_map/populations.snps.vcf.gz",
                comment = "#") %>%
  select(CHROM, POS)  %>%
  filter(!grepl("unloc|MT|scaff", CHROM))

# tmp <- filter(rad, CHROM == "SUPER_X") 
# View(tmp)
# table(rad$FILTER)
# table(rad$CHROM)


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
  mutate(seqnames = Sequence)


### Create base genome ranges ####
#### Genomic ranges
gnm_gr <- Seqinfo(seqnames = seq_sizes$Sequence, seqlengths=seq_sizes$size,
                  isCircular = rep(F, nrow(seq_sizes)), genome = "Auanema")

elim_seqs <- filter(seq_sizes, grepl("loc|scaffold|MT", Sequence)) %>%
  mutate(start = 1, end = size-1) %>%
  select(-size) %>%
  rbind(GRS)

# GRS regions
GRS_gr <- GRanges(seqnames = Rle(elim_seqs$Sequence), 
                  ranges = IRanges(start = elim_seqs$start,
                                   end = elim_seqs$end), 
                  strand = "*")

### Plot RAD ####
### Genome wide coverage of RAD SNPs
window_size <- 100000
aua_1kb <- tileGenome(gnm_gr, tilewidth = window_size, cut.last.tile.in.chrom = T)
gr1_cov <- GenomicRanges::coverage(rad_gr)
clust_abundance <- GenomicRanges::binnedAverage(aua_1kb, gr1_cov, "binned_cov") %>%
  as_tibble()

# mutate(seq_sizes, start = 0) %>%
#   pivot_longer(!Sequence, names_to = "varia", values_to = "POS") %>%
#   rename(CHROM = Sequence) %>%
#   select(-varia) %>%
#   rbind(select(rad, CHROM, POS)) %>%
select(rad, CHROM, POS)  %>%
  filter(!grepl("unloc|MT|scaff", CHROM)) %>%
  ggplot(aes(x=POS/1e6, y = CHROM)) +
  scale_x_continuous(n.breaks = 6,
                     labels = scales::comma) +
  geom_point(alpha = 0.5) +
  xlab("Physical position (Mb)") +
  ylab("Chromosome") + 
  theme_classic()

ggsave("report/figures/Auanema_rhodense_RAD-dots.jpg",
       width = 7, height = 6)

rad_gr <- mutate(rad, chrom=CHROM,
                 start = POS,
                 end = POS+1) %>%
  select(chrom,
         start,
         end) %>%
  makeGRangesFromDataFrame(seqinfo = gnm_gr, keep.extra.columns = T)


rad_plt <- 
  bind_rows(select(lGRS,
                   seqnames = Sequence,
                   diminution_pos = POS),
            clust_abundance) %>%
  filter(!grepl("MT|unloc|scaff", seqnames)) %>%
  mutate(seqnames = paste0("chr", sub("SUPER_", "", seqnames))) %>%
  ggplot(aes(x=start/1e6, y=binned_cov, diminution_pos)) + 
  facet_grid(seqnames ~ .) +
  geom_bar(position="stack", stat="identity", width = window_size/1e6, alpha = 0.8) +
  geom_vline(aes(xintercept = diminution_pos/1e6), linetype="dotted",
             color = "black", size = 1) +
  scale_x_continuous(n.breaks = 6,
                     labels = scales::comma) +
  scale_y_continuous(breaks = scales::pretty_breaks(3), limits = c(0, 0.002)) +
  # scale_fill_brewer(palette="Set1") +
  ylab("SNP density per 100 Kb window") +
  xlab("Physical position (Mb)") +
  theme_bw()

ggsave("report/figures/nxAuaRhod1.1_SNP_density.pdf",
       rad_plt, width = 9, height = 6)


### ncRNAs ####
### Separate sub types
# RNAmmer 
rrna_df <- mutate(rrna_df, gtypeA = "rRNA",
                  gtypeB = ifelse(grepl("=8s", ID), "5S",
                                  "28S_and_18S"))
# tRNAscan 
trna_df <- mutate(trna_df, tRNA_type = sub("ID=([^0-9]+)[0-9].+", "\\1", ID),
         gtypeA = "tRNA",
         gtypeB = ifelse(tRNA_type %in% c("tRNA-His", "tRNA-Ser",
                                          "tRNA-Asp", "tRNA-Arg",
                                          "tRNA-Ile", "tRNA-Gly"),
                         tRNA_type, "other_tRNA"))

# infernal files
miscrna_df <- mutate(miscrna_df, evalue = sub("evalue=([^;]+).+", "\\1", ID),
         description = sub(".+desc=([^;]+)", "\\1", ID),
         misc_type = sub("([^_]+)_.*", "\\1", description),
         gtypeA = "misc",
         gtypeB = ifelse(misc_type %in% c("U2", "U1"),
                         misc_type, "other_misc"))

count(miscrna_df, misc_type) %>% arrange(desc(n))
count(trna_df, tRNA_type) %>% arrange(desc(n)) %>% View



# RNA granges
rrna_gr <- makeGRangesFromDataFrame(rrna_df, keep.extra.columns=T, seqinfo=gnm_gr)
trna_gr <- makeGRangesFromDataFrame(trna_df, keep.extra.columns=T, seqinfo=gnm_gr)
miscrna_gr <- makeGRangesFromDataFrame(miscrna_df, keep.extra.columns=T, seqinfo=gnm_gr)

# trna_gr[trna_gr %over% GRS_gr] %>%
#   as_tibble() %>% View

all_nc <- c(rrna_gr, trna_gr, miscrna_gr)

# Plot density of tRNAs and rRNAs
window_size <- 100000
aua_1kb <- tileGenome(gnm_gr, tilewidth = window_size, cut.last.tile.in.chrom = T)

### Get coverage signal as Rle object

### Get average coverage in each bin
# (since the bins are 1-bp wide, this just keeps the original coverage value)
seq_factors <- rev(c("tRNA", "rRNA", "misc"))

seq_factors <- rev(c("5S", "28S_and_18S", "other_tRNA", 
                     "tRNA-Ile", "tRNA-His", "tRNA-Gly",
                     "tRNA-Ser", "tRNA-Arg", "tRNA-Asp",
                     "U2", "U1", "misc"))
# tRNA-Ala

var_to_check <- "gtypeA"

clust_abundance <- tibble(seqnames = character(), start = integer(), end = integer(),
                          width = integer(), strand = factor(), binned_cov = double(),
                          clust_id = character())
for (clust_id in seq_factors) {
  gr1_cov <- GenomicRanges::coverage(GenomicRanges::reduce(all_nc[mcols(all_nc)[, var_to_check] == clust_id]))
  clust_abundance <- GenomicRanges::binnedAverage(aua_1kb, gr1_cov, "binned_cov") %>%
    as_tibble() %>%
    mutate(clust_id = clust_id) %>%
    rbind(clust_abundance)
  
}

clust_abundance <- mutate(clust_abundance, 
                          clust_id = factor(clust_id, 
                                            levels = seq_factors))

upper_dens <- 0.2
plot_pal <- custom_palette # three_custom_palette
bind_rows(dplyr::select(lGRS, seqnames = Sequence,
                        diminution_pos = POS),
          clust_abundance) %>%
  filter(!grepl("MT|unloc|scaff", seqnames)) %>%
  mutate(seqnames = paste0("chr", sub("SUPER_", "", seqnames)),
         binned_cov = ifelse(binned_cov > upper_dens, upper_dens, binned_cov)) %>%
  ggplot(aes(x=start/1e6, y=binned_cov, fill = clust_id, diminution_pos)) + 
  facet_grid(seqnames ~ .) +
  geom_bar(position="stack", stat="identity", width = window_size/1e6, alpha = 0.8) +
  geom_vline(aes(xintercept = diminution_pos/1e6), linetype="dotted",
             color = "black", size = 1) +
  scale_fill_manual(values = plot_pal) +
  # scale_x_continuous(labels = label_number_si()) +
  scale_y_continuous(breaks = scales::pretty_breaks(3), limits = c(0, upper_dens)) +
  # scale_fill_brewer(palette="Set1") +
  ylab("Density per 100 Kb window") +
  xlab("Physical position (Mb)") +
  theme_bw()


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
  scale_x_continuous(labels = label_number_si()) +
  scale_y_continuous(breaks = scales::pretty_breaks(3), limits = c(0, 1)) +
  # scale_fill_brewer(palette="Set1") +
  ylab("fraction per 100 Kb window") +
  xlab("Position along chromosome") +
  theme_bw()

# ggsave("report/paper2_prelim_20220527/figures/exploratory/repeats/tandem/nxAuaRhod1.1_colored_by_tandem_repeat_many_reps_with_5.trf.pdf",
#        trf_plt, width = 9, height = 6)




