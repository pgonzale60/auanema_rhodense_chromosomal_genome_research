library(GenomicRanges)
library(rtracklayer)
library(scales)
library(tidyverse)
library(cowplot)
library(ggpubr)
library(patchwork)


# tRNA files
t_files <- list.files("analyses/genome_features/repeats/tRNAscan/",
                        full.names = T, pattern = ".hc_trnas.out")
names(t_files) <- make.names(sub(".+//(.+).hc_trnas.out", "\\1", t_files))
hc_atrna_df <- map_df(t_files, readr::read_delim,
                      skip = 3, 
                      trim_ws = T,
                      col_names = c("seqnames", "trna_num", "start", "end", "type", "anticodon", 
                                    "intron_start", "intron_end", "score", "HMM_score", "2pstr_score", "isotype_CM",
                                    "isotype_score","note"),
                    .id = "assembly") %>%
  mutate(multispecies_sequence = paste0(assembly, "_", seqnames),
         tstart = ifelse(start < end, start, end),
         tend = ifelse(start > end, start, end)) %>%
  select(-start, -end, -seqnames) %>%
  rename(start= tstart, end = tend)

# Sequence sizes
fai_files <- list.files("analyses/genome_features/sequence_sizes/",
                        full.names = T, pattern = ".primary.fa.gz.fai")
names(fai_files) <- make.names(sub(".+//(.+).primary.fa.gz.fai", "\\1", fai_files))
seq_sizes <- map_df(fai_files, read_tsv,
                    col_names = c("Sequence", "size", "L1", "L2", "L3"),
                    col_types = c("ciiii"),
                    .id = "assembly") %>%
  select(assembly, Sequence, size) %>%
  filter(!grepl("MT", Sequence)) %>%
  mutate(multispecies_sequence = paste0(assembly, "_", Sequence))

# GRS files
GRS_files <- list.files("analyses/genome_features/elim_coords/",
                        full.names = T, pattern = "nx.*.GRS.bed$")
names(GRS_files) <- make.names(sub(".+//(.+).GRS.bed", "\\1", GRS_files))
GRS <- map_df(GRS_files, read_tsv,
              col_names = c("Sequence", "start", "end"),
              col_types = c("cii"),
              .id = "assembly") %>%
  mutate(multispecies_sequence = paste0(assembly, "_", Sequence))

## Genomic ranges ####
gnm_gr <- Seqinfo(seqnames = seq_sizes$multispecies_sequence, seqlengths=seq_sizes$size,
                  isCircular = rep(F, nrow(seq_sizes)), genome = "Nemas")

GRS_gr <- dplyr::rename(GRS, chrom = multispecies_sequence) %>%
  makeGRangesFromDataFrame(seqinfo = gnm_gr)

trnas_gr <- mutate(hc_atrna_df, seqnames=multispecies_sequence) %>%
  makeGRangesFromDataFrame(seqinfo = gnm_gr, keep.extra.columns = T)

trnas_gr$elim <- trnas_gr %over% GRS_gr

## Filter data ####
# Use the high confidence set
atrna_df <- as_tibble(trnas_gr) #%>%
  filter(note == "high confidence set")

# Explore data ####


count(atrna_df, assembly, anticodon) %>%
  arrange(desc(n))

count(atrna_df, assembly, note == "high confidence set") %>%
    arrange(desc(n))

### Compare species ####
count_data <- atrna_df %>%
  group_by(assembly, type, anticodon) %>%
  summarize(count = n(), .groups = 'drop') %>%
  arrange(type, anticodon) %>%
  mutate(aaa = paste(anticodon, type),
         aaa = factor(aaa,
                      levels = rev(unique(aaa))),
         logdist = log10(5e6)) # log10(1e7) is used as the x-coordinate for the labels

tmp_data <- group_by(atrna_df, anticodon, seqnames) %>%
  arrange(start) %>%
  mutate(dist = c(1, diff(start))) %>%
  ungroup() %>%
  # filter(dist > 1e4) %>%
  filter(dist > 1) %>%
  arrange(type, anticodon) %>%
  mutate(aaa = paste(anticodon, type),
         aaa = factor(aaa,
                      levels = levels(count_data$aaa)))

combined_data <- full_join(tmp_data, count_data, by = c("aaa", "assembly")) 

tRNA_dist_comparison <- ggplot(combined_data, aes(y = aaa)) +
  facet_grid(. ~ assembly) +
  geom_jitter(data = combined_data, aes(x = log10(dist), colour = elim), alpha = 0.1, width = 0.2) +
  geom_text(data = filter(combined_data, !duplicated(paste(aaa, assembly))), 
            aes(x = logdist, label = count), hjust = -0.1, size = 3, color = "#666666") +
  scale_colour_manual(name = "Eliminated", 
                      values = c("TRUE" = "#1f78b4", "FALSE" = "#ff7f00"),
                      labels = c("no", "yes")) +
  ylab("") +
  xlab("Distance between closest copies (log 10)") +
  theme_bw() +
  theme(legend.position = "bottom")


ggsave("report/figures/tRNA_dist_comparison_unfiltered.pdf",
       tRNA_dist_comparison,
       width = 8, height = 13)

### Compare between elimination ####
eatrna_df <- atrna_df %>%
  filter(assembly == "nxAuaRhod1_1") %>%
  mutate(high_conf = note == "high confidence set")
ecount_data <- eatrna_df %>%
  group_by(high_conf, type, anticodon) %>%
  summarize(count = n(), .groups = 'drop') %>%
  arrange(type, anticodon) %>%
  mutate(aaa = paste(anticodon, type),
         aaa = factor(aaa,
                      levels = rev(unique(aaa))),
         logdist = log10(5e6)) # log10(1e7) is used as the x-coordinate for the labels

etmp_data <- group_by(eatrna_df, anticodon, seqnames) %>%
  arrange(start) %>%
  mutate(dist = c(1, diff(start))) %>%
  ungroup() %>%
  filter(dist > 1) %>%
  arrange(type, anticodon) %>%
  mutate(aaa = paste(anticodon, type),
         aaa = factor(aaa,
                      levels = levels(ecount_data$aaa)))


# Define a custom labeller function
confidence_labeller <- as_labeller(c(`TRUE` = "High confidence", 
                             `FALSE` = "Dubious"))


ecombined_data <- full_join(etmp_data, ecount_data, by = c("aaa", "high_conf")) %>%
  mutate(high_conf = factor(high_conf, levels = c("TRUE", "FALSE"))) %>%
  arrange(desc(elim))

etRNA_dist_comparison <- ggplot(ecombined_data, aes(y = aaa)) +
  facet_grid(. ~ high_conf, labeller = as_labeller(confidence_labeller)) +
  geom_jitter(data = ecombined_data, aes(x = log10(dist), colour = elim), alpha = 0.1, width = 0.2) +
  geom_text(data = filter(ecombined_data, !duplicated(paste(aaa, high_conf))), 
            aes(x = logdist, label = count), hjust = -0.1, size = 3, color = "#666666") +
  scale_colour_manual(name = "Localization", 
                      values = c("TRUE" = "#1f78b4", "FALSE" = "#ff7f00"),
                      labels = c("Core", "Eliminated")) +
  ylab("") +
  xlab("Distance between closest copies (log 10)") +
  theme_bw() +
  theme(legend.position = "bottom")


ggsave("report/figures/tRNA_dist_by_confidence.pdf",
       etRNA_dist_comparison,
       width = 8, height = 13)

## Plot data ####

atrna_df_dist <- filter(atrna_df) %>% # , is.na(note)
  group_by(assembly, anticodon, seqnames) %>%
  arrange(start) %>%
  mutate(dist = c(1, diff(start))) %>%
  ungroup() #%>%
  # filter(dist > 1)

# calculate quantiles for each assembly
quantiles <- atrna_df_dist %>%
  group_by(assembly) %>%
  summarise(q90 = quantile(dist, 0.9), 
            q95 = quantile(dist, 0.95),
            q99 = quantile(dist, 0.99))

# join the quantiles back to the main dataframe
atrna_df_dist <- left_join(atrna_df_dist, quantiles, by = "assembly")

p <- mutate(atrna_df_dist, assembly = ifelse(grepl("ipu", assembly),
                                             "O. tipulae", "A. rhodensis")) %>%
  filter(dist > 1) %>%
  ggplot() +
  geom_histogram(aes(x = dist), bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +
  geom_vline(data = . %>% filter(!duplicated(assembly)), aes(xintercept = q90), linetype = "dashed", color = "red") +
  geom_vline(data = . %>% filter(!duplicated(assembly)), aes(xintercept = q95), linetype = "dotted", color = "blue") +
  geom_vline(data = . %>% filter(!duplicated(assembly)), aes(xintercept = q99), linetype = "solid", color = "green") +
  geom_text(data = . %>% filter(!duplicated(assembly)), aes(x = q90, y = 2000, label = round(q90, 0)), angle = 90, size = 3) +
  geom_text(data = . %>% filter(!duplicated(assembly)), aes(x = q95, y = 2000, label = round(q95, 0)), angle = 90, size = 3) +
  geom_text(data = . %>% filter(!duplicated(assembly)), aes(x = q99, y = 2000, label = round(q99, 0)), angle = 90, size = 3) +
  scale_y_log10(expand = c(0, 0.01)) +
  scale_x_continuous(trans = "log10") +
  facet_wrap(. ~ assembly, scales = "fixed", ncol = 2) +
  ylab("Frequency") +
  xlab("Distance between copies of the same isotype (log scale)") +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(face = "italic"))

t_sp_conf_counts <- atrna_df %>%
  mutate(high_conf = ifelse(note == "high confidence set",
                            "High_confidence", 
                            "Dubious"),
         assembly = ifelse(grepl("ipu", assembly),
                           "Oti", "Arh")) %>%
  group_by(assembly, type, anticodon, high_conf) %>%
  summarize(asm_rel_ori = n(), .groups = 'drop') %>%
  group_by(assembly, type, anticodon) %>%
  mutate(tot_ori = sum(asm_rel_ori)) %>%
  ungroup() %>%
  pivot_wider(id_cols = c("assembly", "type", "anticodon", "tot_ori"),
              names_from = "high_conf", values_from = "asm_rel_ori",
              values_fill = 0) %>%
  mutate(perc_hc = round((High_confidence/tot_ori)*100, 1)) %>% 
  select(-Dubious, -High_confidence) %>%
  pivot_wider(id_cols = c("type", "anticodon"),
              names_from = "assembly",
              values_from = c("perc_hc", "tot_ori"),
              values_fill = 0) %>%
  filter(tot_ori_Arh > 10 | tot_ori_Oti > 10) #%>% # remove small types

plot_perc_conf <- pivot_longer(t_sp_conf_counts,
                               cols = perc_hc_Arh:tot_ori_Oti,
                               names_to = c(".value", "assembly"),
                               names_pattern = "(perc_hc|tot_ori)_(Arh|Oti)") %>%  
  ggplot(aes(x = tot_ori, y = perc_hc, color = assembly)) +
  geom_point(alpha = 0.5) +
  scale_y_continuous(expand = c(0, 1)) +
  scale_x_log10(limits = c(1, 2200),
                expand = c(0, 0)
  ) +
  theme_bw() +
  ylab("% of high confidence tRNA predictions") +
  xlab("# of tRNA predictions") +
  theme(legend.position = "bottom")

# Combine plots
p_combined <- p / plot_perc_conf + plot_annotation(tag_levels = 'A')

ggsave("report/figures/compare_tRNA_distances_with_quantiles.jpeg", plot = p_combined,
       width = 7, height = 7)
ggsave("report/figures/compare_tRNA_distances_with_quantiles.pdf", plot = p,
       width = 6, height = 6)




group_by(atrna_df, anticodon, seqnames) %>%
  arrange(start) %>%
  mutate(dist = c(1, diff(start))) %>%
  ungroup() %>% 
  mutate(dist_group = cut(dist,
                          breaks = c(26, 727, 1491, 26840, 11859728),
                          include.lowest = TRUE)) %>%
  filter(!is.na(dist_group)) %>%
  count(anticodon, type, dist_group) %>%
  arrange(type, anticodon) %>%
  mutate(aaa = paste(anticodon, type),
         aaa = factor(aaa,
                      levels = rev(unique(aaa)))) %>%
  ggplot(aes(x = dist_group, y = aaa)) +
  geom_tile(aes(fill = log(n))) +
  scale_fill_gradientn(colors = c("blue", "green", "red")) +
  theme_minimal()
View


## Explore tRNA distances
tRNA_dists <- group_by(atrna_df, anticodon, seqnames) %>%
  arrange(start) %>%
  mutate(block = cumsum(c(0, diff(start) > 1000))) %>%
  ungroup() %>%
  group_by(anticodon, seqnames, block) %>%
  mutate(dist = c(0, diff(start)) ) %>%
  summarise(type = unique(type),
            regionStart = min(start),
            regionEnd = max(start),
            width = regionEnd-regionStart,
            blocksize = n(),
            avgintDist = mean(dist),
            mdnDist = median(dist),
            .groups = "drop") %>%
  arrange(seqnames, regionStart) %>%
  mutate(dist = c(1, diff(regionStart))) %>%
  arrange(type, anticodon) %>%
  mutate(aaa = paste(anticodon, type),
         aaa = factor(aaa,
                      levels = rev(unique(aaa))))

filter(tRNA_dists, blocksize > 12) %>% View

ggplot(tRNA_dists, aes(y = aaa, 
                       x = blocksize)) +
  geom_jitter(alpha = 0.4, width = 0.2)

ggplot(atrna_df_dist, aes(y = anticodon, 
                       x = dist)) +
  geom_jitter(alpha = 0.4, width = 0.2) +
  ylab("") +
  xlab("Distance between closests tRNA isotype copies")


ggplot(tRNA_dists, aes(y = log10(mdnDist), 
                       x = blocksize)) +
  geom_point(alpha = 0.4)


# Main scatter plot
p1 <- tRNA_dists %>%
  filter(blocksize > 12) %>%
  arrange(blocksize) %>%
  ggplot(aes(y = aaa, x = mdnDist, color = blocksize)) +
  geom_jitter(alpha = 0.6, width = 0.2) +
  scale_color_gradientn(colors = c("blue", "green", "red")) +
  ylab("") +
  xlab("Median distance between tRNA copies within cluster") +
  labs(color = "Clustsize") +
  theme_minimal() +
  theme(legend.position = c(0.99, 0.05),
        legend.justification = c(1, 0),
        legend.box.background = element_rect(fill = "white", color = "black"),
        legend.margin = margin(0.1, 0.1, 0.0, 0.1, "cm"))

# Weighted histogram
p2 <- tRNA_dists %>%
  filter(blocksize > 12) %>%
  ggplot(aes(x = mdnDist, weight = blocksize)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +
  theme_minimal() +
  theme(axis.title.x = element_blank()) +
  xlab("")

# Combine plots
p_combined <- p2 / p1 +
  plot_layout(heights = c(1, 4))

print(p_combined)

ggsave(plot = p_combined, filename = "report/figures/tRNA_distance_within_cluster.pdf",
       width = 8, height = 8)



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


# Export data #####
# library(gridExtra)
library(flextable)
# pdf("~/Downloads/tmp.pdf")
mytable <- ecount_data %>%
  select(-aaa, -logdist) %>%
  mutate(high_conf = ifelse(high_conf, "High_confidence", 
                            "Dubious")) %>%
  pivot_wider(id_cols = c("type", "anticodon"),
              names_from = "high_conf", values_from = "count",
              values_fill = 0)

mytable %>%
  mutate(tot_ori = High_confidence + Dubious,
         perc_dub_ori = round((Dubious/tot_ori)*100, 1)) %>%
  filter(tot_ori > 10) %>%
  View

mytable <- mytable %>%
  regulartable() %>%
  fontsize(size = 10)  # change the font size
# Save the table in a pdf file

print(mytable, preview = "pdf")

pdf_table(mytable, filename = "~/Downloads/tmp.pdf")



mytable <- fontsize(mytable, size = 10)
  grid.table()
dev.off()


pdf_document <- 
  kable("latex", booktabs = T) %>%  # Convert table to LaTeX format
  kable_styling(latex_options = "striped", full_width = F)
  
# Save table to a standalone pdf file
cat(c("\\documentclass{standalone}\\usepackage{booktabs}\\begin{document}", 
      pdf_document, 
      "\\end{document}"), 
    sep="\n", 
    file="table.tex")

system("pdflatex table.tex") 





