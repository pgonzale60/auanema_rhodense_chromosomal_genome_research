bln <- read_delim("/Users/pg17/Library/CloudStorage/GoogleDrive-pg17@sanger.ac.uk/My Drive/Documents/Investigacion/A_rhodensis_thesis_chapt/analyses/genome_features/repeats/TE-Aid/blastn.txt", delim = " ") %>%
  mutate(t1 = ifelse(V9 < V10, V9, V10),
         t2 = ifelse(V9 > V10, V9, V10)) %>%
  select(-V9, -V10) %>%
  rename(V9 = t1, V10 = t2)

argT <- read_delim("/Users/pg17/Library/CloudStorage/GoogleDrive-pg17@sanger.ac.uk/My Drive/Documents/Investigacion/A_rhodensis_thesis_chapt/analyses/genome_features/repeats/TE-Aid/tRNA_Arg_ACG.one_copy_tRNA.txt", delim = " ")
subA <- read_delim("/Users/pg17/Library/CloudStorage/GoogleDrive-pg17@sanger.ac.uk/My Drive/Documents/Investigacion/A_rhodensis_thesis_chapt/analyses/genome_features/repeats/TE-Aid/tRNA_Arg_ACG.one_copy_assoc_repeat_subA.txt", delim = " ")
subB <- read_delim("/Users/pg17/Library/CloudStorage/GoogleDrive-pg17@sanger.ac.uk/My Drive/Documents/Investigacion/A_rhodensis_thesis_chapt/analyses/genome_features/repeats/TE-Aid/tRNA_Arg_ACG.one_copy_assoc_repeat_subB.txt", delim = " ")


hc_atrna_df <- readr::read_delim(file = "analyses/genome_features/repeats/tRNAscan/nxAuaRhod1_1.hc_trnas.out",
                              skip = 3, 
                              trim_ws = T,
                              col_names = c("seqnames", "trna_num", "start", "end", "type", "codon", 
                                            "intron_start", "intron_end", "score", "HMM_score", "2pstr_score", "isotype_CM",
                                            "isotype_score","note"))

count(hc_atrna_df, note) %>%
  arrange(desc(n))

atrna_df <- readr::read_delim(file = "analyses/genome_features/repeats/tRNAscan/nxAuaRhod1_1.trnas.txt",
                            skip = 3, 
                            trim_ws = T,
                            col_names = c("seqnames", "trna_num", "start", "end", "type", "codon", 
                                                    "intron_start", "intron_end", "score", "note"))

atrna_df <- filter(hc_atrna_df, note == "high confidence set")


# Remove redundant entries corresponding to introns
atrna_df 



filter(bln, V4 > 1000) %>%
  group_by(V2) %>%
  summarise(minscore = min(V4),
            maxscore = max(V4),
            meanscore = mean(V4),
            q1pos = round(quantile(V9, 0.25)/1e6),
            q2pos = round(quantile(V9, 0.5)/1e6),
            q3pos = round(quantile(V9, 0.75)/1e6),
            # modroundpos = round(V9/1e6),
            n = n())
summary(bln)


filter(bln, V4 < 1000) %>%
  group_by(V2) %>%
  summarise(minscore = min(V4),
            maxscore = max(V4),
            meanscore = mean(V4),
            q1pos = round(quantile(V9, 0.25)/1e6),
            q2pos = round(quantile(V9, 0.5)/1e6),
            q3pos = round(quantile(V9, 0.75)/1e6),
            # modroundpos = round(V9/1e6),
            n = n())

filter(bln, V4 < 1000) %>% count(V4) %>%
  filter(n > 10) %>%
  arrange(desc(n))

filter(bln, V4 < 100000) %>%
  group_by(V4) %>%
  mutate(n_same_len = n()) %>%
  ungroup() %>%
  # filter(n_same_len > 100) %>%
  group_by(V2) %>%
  arrange(V9) %>%
  mutate(block = cumsum(c(1, diff(V9) > 20000)),
         dist = c(1, diff(V9)) ) %>%
  ungroup() %>%
  group_by(V2, block) %>%
  summarise(V2 = unique(V2),
            regionStart = min(V9),
            regionEnd = max(V10),
            width = regionEnd-regionStart,
            blocksize = n(),
            .groups = "drop") %>%
  arrange(V2, regionStart) %>%
  mutate(dist = c(1, diff(regionStart))) %>%
  View
select(target_name, target_length, strand, block, teloPos, regSupport)

arrange(bln, V2, desc(V9)) %>%
  mutate(dist_per_match = c(0, diff(V9))) %>%
  count(dist_per_match) %>%
  filter(n > 10) %>%
  arrange(desc(n)) %>%
  View


## How clustered are tRNAs? ####
group_by(atrna_df, codon, seqnames) %>%
  arrange(start) %>%
  mutate(dist = c(1, diff(start))) %>%
  ungroup() %>%
  filter(dist > 1e4) %>%
  arrange(type, codon) %>%
  mutate(aaa = paste(codon, type),
         aaa = factor(aaa,
                      levels = rev(unique(aaa)))) %>%
  ggplot(aes(y = aaa, x = log10(dist))) +
  geom_jitter(alpha = 0.1, width = 0.2)

group_by(atrna_df, codon, seqnames) %>%
  arrange(start) %>%
  mutate(dist = c(1, diff(start))) %>%
  ungroup() %>%
  count(dist > 1e3)

atrna_df_dist <- filter(atrna_df) %>% # , is.na(note)
  group_by(codon, seqnames) %>%
  arrange(start) %>%
  mutate(dist = c(0, diff(start))) %>%
  ungroup() %>%
  filter(dist > 0)

quantiles <- quantile(atrna_df_dist$dist, c(0.9, 0.95, 0.99))

ggplot(atrna_df_dist, aes(x = dist)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +
  scale_y_log10(expand = c(0, 0.01)) +
  scale_x_log10() +
  geom_vline(aes(xintercept = quantiles[1]), linetype = "dashed") +
  geom_vline(aes(xintercept = quantiles[2]), linetype = "dotted") +
  geom_vline(aes(xintercept = quantiles[3]), linetype = "solid") +
  geom_line(aes(y = 0, linetype = "90th percentile")) +
  geom_line(aes(y = 0, linetype = "95th percentile")) +
  geom_line(aes(y = 0, linetype = "99th percentile")) +
  ylab("Frequency") +
  xlab("Distance between copies of the same isotype") +
  labs(linetype = NULL) +
  theme_minimal() +
  theme(legend.position = c(0.99, 0.7),
        legend.justification = c(1, 0),
        legend.box.background = element_rect(fill = "white", color = "black"),
        legend.margin = margin(0.0, 0.1, 0.0, 0.1, "cm")) +
  guides(linetype = guide_legend(override.aes = list(size = 1))) +
  scale_linetype_manual(values = c("dashed", "dotted", "solid"))


group_by(atrna_df, codon, seqnames) %>%
  arrange(start) %>%
  mutate(dist = c(1, diff(start))) %>%
  ungroup() %>% 
  mutate(dist_group = cut(dist,
                     breaks = c(26, 727, 1491, 26840, 11859728),
                     include.lowest = TRUE)) %>%
  filter(!is.na(dist_group)) %>%
  count(codon, type, dist_group) %>%
  arrange(type, codon) %>%
  mutate(aaa = paste(codon, type),
         aaa = factor(aaa,
                      levels = rev(unique(aaa)))) %>%
  ggplot(aes(x = dist_group, y = aaa)) +
  geom_tile(aes(fill = log(n))) +
  scale_fill_gradientn(colors = c("blue", "green", "red")) +
  theme_minimal()
  View
  
three_reps <- mutate(argT, type = "tRNA") %>%
    bind_rows(mutate(subA, type = "repA")) %>%
    bind_rows(mutate(subB, type = "repB"))  %>%
  mutate(t1 = ifelse(V9 < V10, V9, V10),
         t2 = ifelse(V9 > V10, V9, V10)) %>%
  select(-V9, -V10) %>%
  rename(V9 = t1, V10 = t2)
  
group_by(three_reps, V4) %>%
    group_by(V2) %>%
    arrange(V9) %>%
    mutate(block = cumsum(c(1, diff(V9) > 1000)),
           dist = c(1, diff(V9)) ) %>%
    ungroup() %>%
    group_by(V2, block) %>%
    summarise(V2 = unique(V2),
              regionStart = min(V9),
              regionEnd = max(V10),
              blocksize = n(),
              tRNA_count = sum(type == "tRNA"),
              repA_count = sum(type == "repA"),
              repB_count = sum(type == "repB"),
              .groups = "drop") %>%
    arrange(V2, regionStart) %>%
    mutate(dist = c(1, diff(regionStart))) %>%
    View
  
## Explore tRNA distances
tRNA_dists <- group_by(atrna_df, codon, seqnames) %>%
  arrange(start) %>%
  mutate(block = cumsum(c(0, diff(start) > 1000))) %>%
  ungroup() %>%
  group_by(codon, seqnames, block) %>%
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
  arrange(type, codon) %>%
  mutate(aaa = paste(codon, type),
         aaa = factor(aaa,
                      levels = rev(unique(aaa))))

filter(tRNA_dists, blocksize > 12) %>% View

ggplot(tRNA_dists, aes(y = aaa, 
           x = blocksize)) +
  geom_jitter(alpha = 0.4, width = 0.2)

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


# Now include rRNA and misc
high_freq_misc <- count(miscrna_df, type) %>% arrange(desc(n)) %>%
  filter(n > 10)
misc_dists <- filter(miscrna_df, type %in% high_freq_misc$type) %>%
  group_by(type, seqnames) %>%
  arrange(start) %>%
  mutate(block = cumsum(c(0, diff(start) > 8000))) %>%
  ungroup() %>%
  group_by(type, seqnames, block) %>%
  mutate(dist = c(0, diff(start)) ) %>%
  summarise(regionStart = min(start),
            regionEnd = max(start),
            width = regionEnd-regionStart,
            blocksize = n(),
            avgintDist = mean(dist),
            mdnDist = median(dist),
            .groups = "drop") %>%
  arrange(seqnames, regionStart) %>%
  mutate(dist = c(1, diff(regionStart)))  %>%
  arrange(type)

  
  
rRNA_dists <- select(rrna_df, -type) %>%
rename(type = gtypeB) %>%
  group_by(type, seqnames) %>%
    arrange(start) %>%
  mutate(block = cumsum(c(0, diff(start) > 8000))) %>%
  ungroup() %>%
  group_by(type, seqnames, block) %>%
  mutate(dist = c(0, diff(start)) ) %>%
    summarise(regionStart = min(start),
              regionEnd = max(start),
              width = regionEnd-regionStart,
              blocksize = n(),
              avgintDist = mean(dist),
              mdnDist = median(dist),
              .groups = "drop") %>%
    arrange(seqnames, regionStart) %>%
    mutate(dist = c(1, diff(regionStart)))  %>%
    arrange(type)

p3 <- bind_rows(misc_dists, rRNA_dists) %>%
  filter(blocksize > 12) %>%
  filter(mdnDist < 4000) %>% # Filter out 28 and 18 S rRNAs
  arrange(blocksize) %>%
  ggplot(aes(y = type, x = mdnDist, color = blocksize)) +
  geom_jitter(alpha = 0.6, width = 0.2) +
  scale_color_gradientn(colors = c("blue", "green", "red")) +
  ylab("") +
  xlab("Median distance between gene copies within cluster") +
  labs(color = "Clustsize") +
  theme_minimal() +
  theme(legend.position = c(0.99, 0.05),
        legend.justification = c(1, 0),
        legend.box.background = element_rect(fill = "white", color = "black"),
        legend.margin = margin(0.1, 0.1, 0.0, 0.1, "cm"))

# Weighted histogram
p4 <- bind_rows(misc_dists, rRNA_dists) %>%
  filter(blocksize > 12) %>%
  filter(mdnDist < 4000) %>%
  ggplot(aes(x = mdnDist, weight = blocksize)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +
  theme_minimal() +
  theme(axis.title.x = element_blank()) +
  xlab("")

p_combined2 <- ((p2 + p4) / (p1 + p3 )) +
  plot_layout(heights = c(1, 4))

print(p_combined2)

ggsave(plot = p_combined2, filename = "report/figures/RNA_distances_within_cluster.pdf",
       width = 12, height = 8)

# Explore specific RNAs for TE-Aid
group_by(atrna_df, codon, seqnames) %>%
  arrange(start) %>%
  mutate(block = cumsum(c(0, diff(start) > 1000))) %>%
  ungroup() %>%
  group_by(codon, seqnames, block) %>%
  mutate(dist = c(0, diff(start)) ) %>%
  ungroup %>%
  filter(codon == "AGT") %>%
  arrange(seqnames, start) %>%
  View()

filter(miscrna_df, type %in% high_freq_misc$type) %>%
  group_by(type, seqnames) %>%
  arrange(start) %>%
  mutate(block = cumsum(c(0, diff(start) > 8000))) %>%
  ungroup() %>%
  group_by(type, seqnames, block) %>%
  mutate(dist = c(0, diff(start))) %>%
  ungroup %>%
  filter(type == "U2") %>%
  arrange(seqnames, start) %>%
  View()


nematoda_tRNAs <- read_tsv("~/Downloads/gtrnadb-search168473.out.gz")
count(nematoda_tRNAs, GenomeID)

euk_tRNAs <- read_tsv("~/Downloads/gtrnadb-search168611.out.gz")
count(euk_tRNAs, Genome) %>% arrange(desc(n))






# Nematode-wide search #####
# NCBI assemblies
library(patchwork)
library(tidyverse)
assm_sp <- read_tsv("raw/assembly_species.tsv",
                    col_names = c("ID", "species")) %>%
  bind_rows(tibble(ID = "nxAuaRhod1_1", species = "Auanema_rhodense"))

# tRNA files
t_files <- list.files("analyses/genome_features/repeats/tRNAscan/nematodes/",
                      full.names = T, pattern = ".hc_trnas.out.gz")
names(t_files) <- make.names(sub(".+//((GCA|GCF)_[0-9]+\\.[0-9]+|nxAuaRhod1_1).*\\.hc_trnas.out.gz", "\\1", t_files))
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
  rename(start= tstart, end = tend) %>%
  left_join(assm_sp, by = c("assembly" = "ID"))

count(hc_atrna_df, species, note == "high confidence set") %>%
  arrange(desc(n)) %>%
  View

filter(hc_atrna_df, species == "Allodiplogaster_sudhausi") %>%
  count(type, anticodon, note == "high confidence set") %>%
  arrange(desc(n)) %>%
  View

hc_atrna_rs <- mutate(hc_atrna_df, high_conf = ifelse(note == "high confidence set",
                          "High_confidence", 
                          "Dubious")) %>%
count(species, high_conf) %>%
  arrange(desc(n)) %>%
  pivot_wider(id_cols = c("species"),
              names_from = "high_conf", values_from = "n",
              values_fill = 0) %>%
  mutate(tot_preds = High_confidence + Dubious,
         perc_hc = round((High_confidence/tot_preds)*100, 1)) 

hc_atrna_rs %>%
  ggplot(aes(x = High_confidence, y = perc_hc)) +
  geom_point(alpha = 0.5) +
  scale_y_continuous(expand = c(0, 1)) +
  scale_x_log10(limits = c(1, 25000),
                expand = c(0, 0)
  ) +
  theme_bw() +
  ylab("% of total tRNA predictions") +
  xlab("# of high confidence tRNA predictions")



atrna_df_dist <- filter(hc_atrna_df, species %in%
                          c("Allodiplogaster_sudhausi",
                            "Auanema_rhodense",
                            "Caenorhabditis_elegans",
                            "Oscheius_tipulae"),
                        note == "high confidence set") %>% # , is.na(note)
  group_by(species, anticodon, multispecies_sequence) %>%
  arrange(start) %>%
  mutate(dist = c(1, diff(start))) %>%
  ungroup() #%>%
# filter(dist > 1)

# calculate quantiles for each assembly
quantiles <- atrna_df_dist %>%
  group_by(species) %>%
  summarise(q90 = quantile(dist, 0.9), 
            q95 = quantile(dist, 0.95),
            q99 = quantile(dist, 0.99))

# join the quantiles back to the main dataframe
atrna_df_dist <- left_join(atrna_df_dist, quantiles, by = "species")

p <- mutate(atrna_df_dist, assembly = sub("([A-Za-z])[A-Za-z]*_([A-Za-z]*)", "\\1. \\2", species)) %>%
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



ggsave("report/figures/compare_tRNA_distances_with_quantiles_four_nemas.pdf",
       width = 8.5, height = 5)

ggsave("~/Downloads/tmp.pdf",
       width = 8.5, height = 5)

t_sp_conf_counts <- filter(hc_atrna_df, species %in%
                             c("Allodiplogaster_sudhausi",
                               "Auanema_rhodense",
                               "Caenorhabditis_elegans",
                               "Oscheius_tipulae")) %>%
  mutate(high_conf = ifelse(note == "high confidence set",
                            "High_confidence", 
                            "Dubious"),
         assembly = sub("([A-Za-z])[A-Za-z]*_([A-Za-z]*)", "\\1_\\2", species)) %>%
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
  filter(rowSums(select(., starts_with("tot_ori")) > 10) > 0) #%>% # remove small types

plot_perc_conf <- pivot_longer(t_sp_conf_counts,
                               cols = perc_hc_A_rhodense:tot_ori_O_tipulae,
                               names_to = c(".value", "assembly"),
                               names_pattern = "(perc_hc|tot_ori)_(A_rhodense|A_sudhausi|C_elegans|O_tipulae)") %>%  
  mutate(assembly = sub("_", ". ", assembly)) %>%
  ggplot(aes(x = tot_ori, y = perc_hc, color = assembly)) +
  geom_point(alpha = 0.5) +
  scale_y_continuous(expand = c(0, 1)) +
  scale_x_log10(limits = c(1, 2800),
                expand = c(0, 0)
  ) +
  theme_bw() +
  ylab("% of high confidence tRNA predictions") +
  xlab("# of tRNA predictions") +
  theme(legend.position = "bottom")


# Combine plots
p_combined <- p / plot_perc_conf + plot_annotation(tag_levels = 'A') +
  plot_layout(heights = c(2, 1))

ggsave("report/figures/compare_tRNA_distances_with_quantiles_and_conf_perc.pdf", plot = p_combined,
       width = 6, height = 8)


