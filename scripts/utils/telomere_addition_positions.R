library(GenomicRanges)
library(rtracklayer)
library(scales)
library(tidyverse)
library(cowplot)
library(ggpubr)

aua_telomil <- read_tsv("analyses/genome_features/elim_coords/nxAuaRhod1.pb.miltel.telomeric.tsv.gz",
                    col_names = c("Sequence", "start", "end", "score",
                                  "orientation", "telomere_gap_average", 
                                  "telomere_senses", "telomere_gaps",
                                  "average_coverage"))

ctg_cov <- read_tsv("analyses/telomeres/telo_hifiasm.mosdepth.summary.txt") %>%
  filter(chrom != "total")
ctg_telos <- read_tsv("analyses/telomeres/NCRF_summary.txt") %>%
  left_join(ctg_cov, by = c("seq" = "chrom", "seqLen" = "length"))
reads_telos <- read_tsv("analyses/telomeres/NCRF_reads_summary.txt")

range(reads_telos$querybp)

# Coverage
covFile <- "analyses/genome_features/nemaChromQC/nxAuaRhod1_1.regions.bed.gz"
cov <- read_tsv(covFile,
                col_names = c("seqnames", "start",
                              "end", "median_cov"))%>%
  filter(!grepl("MT", seqnames)) %>%
  mutate(start = start +1)
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


## Genomic ranges ####
gnm_gr <- Seqinfo(seqnames = seq_sizes$Sequence, seqlengths=seq_sizes$size,
                  isCircular = rep(F, nrow(seq_sizes)), genome = "Auanema")

break_sites_gr <- filter(lGRS, POS > 10,
         POS < 20e6) %>%
  mutate(chrom=Sequence,
       end = POS + 1) %>%
  select(chrom,
         start = POS,
         end) %>%
  makeGRangesFromDataFrame(seqinfo = gnm_gr, keep.extra.columns = T) %>%
  GenomicRanges::trim()

telomil_gr <- mutate(aua_telomil, chrom=Sequence) %>%
  select(chrom,
         start,
         end,
         telomere_gap_average,
         average_coverage) %>%
  makeGRangesFromDataFrame(seqinfo = gnm_gr, keep.extra.columns = T) %>%
  GenomicRanges::trim()

cov_gr <- makeGRangesFromDataFrame(cov, keep.extra.columns=T,
                                   seqinfo=gnm_gr)

telomil_gr$atBreakSite <- telomil_gr %over% break_sites_gr

tm_cov <- plyranges::join_overlap_inner(telomil_gr, cov_gr)


## Plots ####

ctgcovplot <- ggplot(ctg_telos, aes(x = querybp, y = mean)) +
  geom_point() +
  theme_bw() +
  xlab("Length of telomere repeat array") +
  ylab("Mean coverage")

readTeloLen <- ggplot(reads_telos, aes(x = querybp)) +
  geom_histogram(binwidth = 200) +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0, 3500)) +
  theme_bw() +
  ylab("Frequency") +
  xlab("Length of telomere repeat array")


clipLenplot <- as_tibble(tm_cov) %>%
  filter(is.finite(telomere_gap_average),
                      !is.na(telomere_gap_average)) %>%
  ggplot(aes(x = median_cov,
             y = telomere_gap_average + 0.5,
             color = atBreakSite)) +
  scale_y_log10(labels = scales::comma,
                expand = c(0,0.1)) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  ylab("Clipped sequence length") +
  xlab("Coverage around soft clip") +
  scale_color_manual(name = "Identified\nbreak site",
                     breaks = c("TRUE", "FALSE"),
                     labels = c("Yes", "No"),
                     values = c("#FF7F0E", "#1F77B4")) +
  guides(color = guide_legend(nrow = 1)) +
  theme_bw() +
  theme(legend.position = "top")

pExpTelo <- ggarrange(ctgcovplot, readTeloLen, 
                      nrow = 2, align = c("v"),
                      labels = LETTERS[1:2])

fexpoTelo <- ggarrange(pExpTelo, clipLenplot, 
                      nrow = 1, labels = c("", "C"))
ggsave("report/figures/telomere_features.pdf",
       fexpoTelo, units = "in", 
       width = 8.5*1.1,
       height = 11*0.6)

filter(ctg_telos, querybp < 2100, mean > 80)
filter(ctg_telos, mean < 80) %>% pull(querybp) %>%
  range()

filter(aua_telomil, average_coverage > 300) %>%
  View



filter(aua_telomil, average_coverage < 200) %>%
  count(round(telomere_gap_average)) %>%
  arrange(desc(n)) %>%
  View

filter(aua_telomil, 
       telomere_gap_average > 1780,
       telomere_gap_average < 1790) %>% View

filter(aua_telomil, 
       telomere_gap_average > 4120,
       telomere_gap_average < 4135) %>% View

filter(aua_telomil, 
       telomere_gap_average > 9478,
       telomere_gap_average < 9489) %>% View

9479-9488


filter(aua_telomil, average_coverage < 200) %>% View

arrange(aua_telomil, Sequence, start) %>%
  View

mutate(aua_telomil, 
       telocov = as.integer(sub("[-\\+]\\*", "", telomere_senses))) %>% 
  # filter(!is.na(cov_at_0), cov_at_0 > 20, grepl("SUPER", Sequence)) %>%
  filter(telocov > 60, grepl("SUPER", Sequence)) %>%
  arrange(Sequence, start) %>%
  View




aua_readtel <- read_tsv("analyses/genome_features/telomere_sizes/NCRF_reads_summary.txt")

hist(aua_readtel$m)




# Create a scatterplot with a log-scale Y axis
p1 <- ggplot(aua_telomil, aes(x = average_coverage,
                              y = telomere_gap_average + 0.5)) +
  scale_y_log10(labels = scales::comma, expand = c(0, 0.1)) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_classic() #+

# Create marginal plots
# Use geom_density/histogram for whatever you plotted on x/y axis 
plot_y <- axis_canvas(p1, axis = "y", coord_flip = TRUE) +
  geom_density(aes(telomere_gap_average + 0.5),
               aua_telomil[aua_telomil$average_coverage < 200, ]) +
  scale_x_log10() +
  coord_flip()
plot_final <- insert_yaxis_grob(p1, plot_y, position = "right")
ggdraw(plot_final)


# Create a histogram of the Y axis variable
p2 <- ggplot(aua_telomil, aes(x = telomere_gap_average + 0.5)) +
  geom_histogram(bins = 30) +
  coord_flip()

# Combine the scatterplot and histogram into a single plot
p_combined <- plot_grid(p1, p2, ncol = 2, align = "v", axis = "tb")

# Display the combined plot
p_combined



x <- read_tsv("~/Downloads/download.tsv")

count(x, family) %>% arrange(desc(n)) %>% View
count(x, class) %>% arrange(desc(n))
