library(tidyverse)

# Load BED file
bed <- read_tsv("analyses/genome_features/repeats/modDotPlot/nxAuaRhod1_1_chroms.bed.gz", 
                           col_names = c("query_name", "query_start", "query_end", 
                                         "reference_name", "reference_start", "reference_end", 
                                         "perID_by_events"), 
                           comment = "#") %>%
  mutate(identity = perID_by_events / 100)  # 0-1 scale

# Calculate min and max identity for heatmap range
id_min <- min(bed$identity)
id_max <- max(bed$identity)

# Find faidx index files ".primary.fa.gz.fai" to get sequence sizes
fai_files <- list.files("analyses/genome_features/sequence_sizes/",
                        full.names = TRUE, pattern = "nxAuaRhod1.1.primary.fa.gz.fai")

# Set names for the files based on the extracted names from the file paths, correspond to assembly ID
names(fai_files) <- make.names(sub(".+//(.+).primary.fa.gz.fai", "\\1", fai_files))

# Read the sequence sizes from the files and combine them into a single data frame
seq_sizes <- map_df(fai_files, read_tsv,
                    col_names = c("Sequence", "size", "L1", "L2", "L3"),
                    col_types = c("ciiii"),
                    .id = "assembly") %>%
  select(assembly, Sequence, size) %>%
  filter(!grepl("MT|unloc|scaf", Sequence)) %>%
  mutate(multispecies_sequence = paste0(assembly, "_", Sequence))

# Calculate cumulative positions
chrom_lengths <- seq_sizes %>%
  mutate(cumstart = c(0, cumsum(size)[-n()]), 
         cumend = cumsum(size))

# Break sites
break_files <- list.files("analyses/curated_GRS_coords/",
                          full.names = T, pattern = "nxAuaRhod1.1.chr_diminutions_sites.tsv")
names(break_files) <- make.names(sub(".+//(.+).chr_diminutions_sites.tsv", "\\1", break_files))
break_sites <- map_df(break_files, read_tsv,
                      col_names = T,
                      .id = "assembly") %>%
  filter(assembly == "nxAuaRhod1.1",
         !grepl("MT|unloc|scaf", chr),
         grepl("telomere", feature)) %>%
  mutate(multispecies_sequence = paste0(assembly, "_", chr))

# Load new BED file with ranges
ranges <- read_tsv("analyses/genome_features/elim_coords/nxAuaRhod1_1.GRS.bed", col_names = c("chrom", "start", "end")) %>%
  left_join(chrom_lengths %>% select(Sequence, cumstart), by = c("chrom" = "Sequence")) %>%
  mutate(start_cum = start + cumstart, 
         end_cum = end + cumstart)

# Bin data using actual bin size (165497 bp)
bin_size <- 165497
bed_binned <- bed %>%
  mutate(x_bin = floor(query_start/ bin_size) * bin_size, 
         y_bin = floor(reference_start / bin_size) * bin_size) %>%
  group_by(x_bin, y_bin) %>%
  summarise(identity = mean(identity), .groups = "drop") %>%
  # Mirror the data for upper triangle
  bind_rows(., mutate(., temp = x_bin, x_bin = y_bin, y_bin = temp) %>% select(-temp))

# Create Hi-C-like heatmap with dynamic identity range
p <- ggplot(bed_binned, aes(x = x_bin, y = y_bin, fill = identity)) +
  geom_tile(width = bin_size, height = bin_size) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = (id_min + id_max) / 2, 
                       limits = c(id_min, id_max), 
                       name = "Identity (%)") +
  scale_x_continuous(labels = function(x) round(x / 1e6, 1),  # Mb, rounded to 1 decimal
                     breaks = c(0, chrom_lengths$cumend),              # Breaks at chromosome ends
                     limits = c(0, max(chrom_lengths$cumend)),
                     expand = c(0, 0)) +
  scale_y_reverse(labels = gsub("SUPER_", "", chrom_lengths$Sequence),     # Mb, reversed
                  breaks = chrom_lengths$cumstart,                # Same breaks
                  limits = c(max(chrom_lengths$cumend), 0),
                  expand = c(0, 0)) +
  labs(x = "Position (Mb)", 
       y = "Position (chromosome)") +
  coord_fixed() +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_line(color = "black", linewidth = 0.5))

# Add custom lines (replace with your positions)
# custom_lines <- 
#   left_join(break_sites,
#             chrom_lengths %>% 
#               select(Sequence, cumstart), by = c("chr" = "Sequence")) %>%
#   mutate(pos_cum = diminution_pos + cumstart)

# p <- p +
#   geom_vline(xintercept = custom_lines$pos_cum, color = "gray50", linewidth = 0.3, alpha = 0.5, linetype = "dashed", ) +
#   geom_hline(yintercept = custom_lines$pos_cum, color = "gray50", linewidth = 0.3, alpha = 0.5, linetype = "dashed")


# Add checkerboard shading (vertical and horizontal)
max_coord <- max(chrom_lengths$cumend)
p <- p +
  # Vertical shading (full height columns)
  geom_rect(data = ranges, 
            aes(xmin = start_cum, xmax = end_cum, ymin = 0, ymax = max_coord), 
            fill = "gray80", alpha = 0.2, inherit.aes = FALSE) +
  # Horizontal shading (full width rows)
  geom_rect(data = ranges, 
            aes(xmin = 0, xmax = max_coord, ymin = start_cum, ymax = end_cum), 
            fill = "gray80", alpha = 0.2, inherit.aes = FALSE)

# Save and display
ggsave("hic_like_plot.png", p, width = 10, height = 10, dpi = 500)

