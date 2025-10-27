library(rtracklayer)
library(tidyverse)

# Classification of repeats

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
  filter(!grepl("MT", Sequence)) %>%
  mutate(multispecies_sequence = paste0(assembly, "_", Sequence))




# GRS files
GRS_files <- list.files("analyses/genome_features/elim_coords/",
                        full.names = T, pattern = "nx.*.GRS.bed$")
names(GRS_files) <- make.names(sub(".+//(.+).GRS.bed", "\\1", GRS_files))
GRS <- map_df(GRS_files, read_tsv,
              col_names = c("Sequence", "start", "end", "has_GRS"),
              col_types = c("ciic"),
              .id = "assembly") %>%
  mutate(multispecies_sequence = paste0(assembly, "_", Sequence))



# Protein files
protfiles <- list.files("analyses/genes/GTFs/",
                        full.names = T, pattern = ".gz")
names(protfiles) <- make.names(sub(".+//(.+).gff.gz", "\\1", protfiles))
genes_df <- map_df(protfiles, read_tsv,
                   comment = "#", quote = "\"",
                   col_names = c("Sequence", "source", "type", "start", "end", "score",
                                 "strand", "phase", "ID"),
                   .id = "assembly") %>%
  filter(!grepl("MT", Sequence),
         !type %in% c("intron", "start_codon", "stop_codon")) %>%
  mutate(seqnames = paste0(assembly, "_", Sequence)) %>%
  select(-assembly, -Sequence)

# RModeler's classification of collapsed (TranposonPSI, gt_LTR and RepeatModeler + nema_repeats)
# repeats
reps_names_files <- list.files("analyses/genome_features/repeats/custom_tblastn_RMasker/",
                               full.names = T, pattern = "repts.dict.tsv")
names(reps_names_files) <- make.names(sub(".+//(.+).collapsed.classfd.repts.dict.tsv.gz", "\\1", reps_names_files))
reps_names <- map_df(reps_names_files, read_tsv,
                     col_names = c("cluster", "assignment"),
                     col_types = c("cc"),
                     .id = "assembly") %>%
  mutate(multispecies_cluster = paste0(assembly, "_", cluster)) %>%
  select(-assembly, -cluster)

# Rmasker of collapsed (TranposonPSI, gt_LTR and RepeatModeler + nema_repeats)
# repeats
rmaskfiles <- list.files("analyses/genome_features/repeats/custom_tblastn_RMasker/",
                         full.names = T, pattern = "gff.gz")
names(rmaskfiles) <- make.names(sub(".+//(.+).fasta.out.gff.gz", "\\1", rmaskfiles))
rmask_df <- map_df(rmaskfiles, read_tsv, comment = "#",
                   col_names = c("Sequence", "source", "type", "start", "end",  "score", "strand",
                                 "phase", "Target"), # score is perc_divergence (% identity from consus)
                   .id = "assembly") %>%
  mutate(seqnames = paste0(assembly, "_", Sequence),
         motif = sub(".+tif:(.+)\".+", "\\1", Target),
         multispecies_cluster = paste0(assembly, "_", motif)) %>%
  filter(!grepl("MT", Sequence)) %>%
  select(-assembly, -Sequence, -motif, -Target)


# TRF files
trffiles <- list.files("analyses/genome_features/repeats/TRF/",
                       full.names = T, pattern = "trf.tsv.gz")
names(trffiles) <- make.names(sub(".+//windowed100K_(.+).trf.tsv.gz", "\\1", trffiles))
trf_df <- map_df(trffiles, read_tsv,
                 .id = "assembly",
                 col_names = c("chrom", "start", "end", "unit_size",
                               "num_copies", "perc_ident", "perc_indels",
                               "score", "entropy", "sequence")) %>%
  mutate(subcoord = as.integer(sub(".+:(.+)-.+", "\\1", chrom)),
    start = start + subcoord, 
    end = end + subcoord, 
    chrom = sub(":.+", "", chrom)) %>%
  mutate(id = 1:n()) %>%
  mutate(seqnames = paste0(assembly, "_", chrom),
         source = "TRF",
         type = "tandem_repeat",
         assignment = "Tandem_repeat") %>%
  filter(!grepl("MT", chrom)) %>%
  select(-assembly, -chrom, -sequence)


# RNAmmer files
rrnafiles <- list.files("analyses/genome_features/repeats/rnammer/",
                        full.names = T, pattern = "gff3.gz")
names(rrnafiles) <- make.names(sub(".+//(.+).rnammer.gff3.gz", "\\1", rrnafiles))
rrna_df <- map_df(rrnafiles, read_tsv,
                  .id = "assembly", col_names = c("Sequence", "source", "type", "start", "end", "score",
                                                  "strand", "phase", "ID")) %>%
  mutate(seqnames = paste0(assembly, "_", Sequence)) %>%
  filter(!grepl("MT", Sequence)) %>%
  select(-assembly, -Sequence, -phase)

# tRNAscan files
trnafiles <- list.files("analyses/genome_features/repeats/tRNAscan/",
                        full.names = T, pattern = "gff.gz")
names(trnafiles) <- make.names(sub(".+//(.+).trnas.gff.gz", "\\1", trnafiles))
trna_df <- map_df(trnafiles, read_tsv,
                  .id = "assembly", col_names = c("Sequence", "source", "type", "start", "end", "score",
                                                  "strand", "phase", "ID")) %>% 
  mutate(seqnames = paste0(assembly, "_", Sequence)) %>%
  filter(!grepl("MT", Sequence)) %>%
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
  mutate(evalue = sub("evalue=([^;]+).+", "\\1", ID),
         description = sub(".+desc=([^;]+)", "\\1", ID),
         seqnames = paste0(assembly, "_", Sequence)) %>%
  filter(!grepl("MT", Sequence)) %>%
  select(-assembly, -Sequence, -phase)

# count(miscrna_df, type) %>% arrange(desc(n)) %>% View()

# Dictionary of sequence names
rmask_df_nmd <- left_join(rmask_df, reps_names, by = c("multispecies_cluster")) %>%
  mutate(assignment = ifelse(is.na(assignment), "Simple_repeat", assignment)) # The "Target" column shows the simple repeat

#### Separating distal from internal ####

unloc_seqs <- filter(seq_sizes, grepl("loc|scaffold_534", multispecies_sequence)) %>%
  mutate(start = 1, end = size-1)
seq_ends <- bind_rows(mutate(seq_sizes, start = 1, end = 1e4),
                      mutate(seq_sizes, start = size - 1e4, end = size)
) %>%
  filter(!grepl("loc|scaffold_534", multispecies_sequence)) %>%
  bind_rows(unloc_seqs)



#### GRanges creation ####
# Whole genome sizes
gnm_gr <- Seqinfo(seqnames = seq_sizes$multispecies_sequence, seqlengths=seq_sizes$size,
                  isCircular = rep(F, nrow(seq_sizes)), genome = "Auanema")

seq_ends_gr <- GRanges(seqnames = Rle(seq_ends$multispecies_sequence), 
                       ranges = IRanges(start = seq_ends$start,
                                        end = seq_ends$end), 
                       strand = "*",
                       seqinfo=gnm_gr)
# GRS regions

GRS_gr <- GRanges(seqnames = Rle(GRS$multispecies_sequence), 
                  ranges = IRanges(start = GRS$start,
                                   end = GRS$end), 
                  strand = "*")
GRS_distal_gr <- GRS_gr[GRS_gr%over% seq_ends_gr]
GRS_distal_gr$type = "distal"
GRS_intrnl_gr <- GRS_gr[!GRS_gr%over% seq_ends_gr]
GRS_intrnl_gr$type = "internal"
# Protein-coding genes
genes_gr <- makeGRangesFromDataFrame(genes_df, keep.extra.columns=T, seqinfo=gnm_gr)


# Rmasker
rmask_gr <- makeGRangesFromDataFrame(rmask_df_nmd, keep.extra.columns=T, seqinfo=gnm_gr)

# TRF
trf_gr <- makeGRangesFromDataFrame(trf_df, keep.extra.columns=T, seqinfo=gnm_gr)

# RNAmmer
rrna_gr <- makeGRangesFromDataFrame(rrna_df, keep.extra.columns=T, seqinfo=gnm_gr)

# tRNAscan
trna_gr <- makeGRangesFromDataFrame(trna_df, keep.extra.columns=T, seqinfo=gnm_gr)

# Infernal for miscelaneious ncRNAs
miscrna_gr <- makeGRangesFromDataFrame(miscrna_df, keep.extra.columns=T, seqinfo=gnm_gr)




#### Intersections ####

# Remove RMasker annots that overlap large tandem repeats
large_trf_gr <- trf_gr[width(trf_gr)>1000]
rmask_gr <- rmask_gr[!rmask_gr %over% large_trf_gr]

# Reclassify RMasker elements that overlap TRF as tandem repeats
rmask_gr$assignment <- ifelse(rmask_gr$assignment %in% c("Unknown", "Simple_repeat", "Satellite") &
                                rmask_gr %over% trf_gr, "Tandem_repeat", rmask_gr$assignment)
# Combine TRF with RMasker
# "reduce" TRF ranges because they are highly redundant
trf_no_rmask_gr <- GenomicRanges::reduce(trf_gr[!trf_gr %over% rmask_gr],
                                         min.gapwidth=150L)
trf_no_rmask_gr$assignment <- "Tandem_repeat"
trf_no_rmask_gr$type <- "tandem_repeat"
trf_no_rmask_gr$source <- "TRF"

# RMasker "Tandem_repeat" are also highly redundant
# Let's separate and reduce
rmask_tr_gr <- GenomicRanges::reduce(rmask_gr[rmask_gr$assignment == "Tandem_repeat"],
                                     min.gapwidth=150L)
rmask_tr_gr$assignment <- "Tandem_repeat"
rmask_tr_gr$type <- "tandem_repeat"
rmask_tr_gr$source <- "RepeatMasker"
rmask_no_tr_gr <- rmask_gr[!rmask_gr %over% rmask_tr_gr]

reps_gr <- c(rmask_no_tr_gr, rmask_tr_gr, trf_no_rmask_gr)

# Give priority to rRNA and tRNAs over TRF with RMasker and combine them
classic_ncRNAs_gr <- c(trna_gr, rrna_gr)
reps_no_trnas_no_rrnas <- GenomicRanges::subtract(reps_gr, classic_ncRNAs_gr,
                                                  ignore.strand=T) %>% 
  unlist()
all_reps_gr <- c(classic_ncRNAs_gr, reps_no_trnas_no_rrnas)

# miscrna_gr[miscrna_gr %over% all_reps_gr]
# all_reps_gr[all_reps_gr %over% miscrna_gr]
# There are 78 intersections between the repeats and misc RNAs predicted by Infernal
# 7I will override repeats because I would be more interested in looking at miscRNA than repeats
all_reps_gr_no_misc_gr <- GenomicRanges::subtract(all_reps_gr, miscrna_gr) %>%
  unlist()
all_nc_gr <- c(all_reps_gr_no_misc_gr, miscrna_gr)

# Check which CDS are classified as high confidence ncRNAs or repeats
cds_gr <- genes_gr[genes_gr$type == "CDS"]
# cds_gr[cds_gr %over% all_nc_gr] %>% as_tibble() %>% count(type) %>% arrange(desc(n))
# all_nc_gr[all_nc_gr %over% cds_gr] %>% as_tibble() %>% count(type) %>% arrange(desc(n))
# all_nc_gr[all_nc_gr %over% cds_gr] %>% as_tibble() %>% count(assignment) %>% arrange(desc(n))
# all_nc_gr[all_nc_gr %over% cds_gr] %>% as_tibble() %>% View()
# Becuase most of the intersection is with repeats, high confidence ranges will focus on repeats
# The few tRNAs and rRNAs that appear seem reliable as to replace gene models

# Intersect reps_nc against CDS
nc_inter_prots_gr <- all_nc_gr[all_nc_gr %over% cds_gr]
# Low confidence repeats are short and simple
unreliable_nc_inter_prots_gr <- nc_inter_prots_gr[width(nc_inter_prots_gr) < 100 & 
                                                    nc_inter_prots_gr$assignment %in% c("Unknown", "Simple_repeat", "Tandem_repeat")]
# Everything else is considered as reliable
reliable_nc_inter_prots_gr <- nc_inter_prots_gr[!nc_inter_prots_gr %over% unreliable_nc_inter_prots_gr]
# reliable_nc_inter_prots_gr %>% as_tibble() %>% View()


# Get the intersected CDS to then get the gene and then remove everything the gene contains
discard_cds_gr <- cds_gr[cds_gr %over% reliable_nc_inter_prots_gr]
genes_gene_gr <- genes_gr[genes_gr$type == "gene"]
discard_genes_gene_gr <- genes_gene_gr[genes_gene_gr %over% discard_cds_gr]
discard_genes_gr <- genes_gr[genes_gr %over% discard_genes_gene_gr]
accepted_genes_gr <- genes_gr[!genes_gr %over% discard_genes_gene_gr]


# How many genes are left per assembly at the end?
stats_accepted_genes <- as_tibble(accepted_genes_gr) %>%
  mutate(assembly = sub("([^_]+)_.+", "\\1", seqnames)) %>%
  filter(type == "gene") %>%
  group_by(assembly) %>%
  count() %>%
  ungroup() %>%
  arrange(n)
# Makes sense


#### GRS intersections ####
GRS_size <- GRS_gr %>% as_tibble() %>%
  mutate(assembly = sub("([^_]+)_.+", "\\1", seqnames)) %>%
  group_by(assembly) %>%
  summarise(n_GRS = n(), span_GRS = sum(width))

wgnm_gr <- GRanges(seqnames = Rle(seq_sizes$multispecies_sequence), 
                   ranges = IRanges(start = 1,
                                    end = seq_sizes$size), 
                   strand = "*",
                   seqinfo = gnm_gr)
core_gr <- unlist(subtract(wgnm_gr, GRS_gr))
core_gr$type <- "non_GRS"
GRS_by_type_size <- c(GRS_distal_gr, GRS_intrnl_gr, core_gr) %>% as_tibble() %>%
  mutate(assembly = sub("([^_]+)_.+", "\\1", seqnames)) %>%
  group_by(assembly, type) %>%
  summarise(span_GRS = sum(width), .groups = "drop") %>%
  pivot_wider(names_from = type, values_from = contains("span"), values_fill = 0, id_expand = T) %>%
  mutate(repeat_fam = "Total") %>%
  select(assembly, repeat_fam, distal_GRS = distal, internal_GRS = internal, non_GRS)


genome_size <- seq_sizes %>% as_tibble() %>%
  group_by(assembly) %>%
  summarise(span = sum(size))

all_nc_gr[all_nc_gr %over% GRS_gr] %>% as_tibble() %>%
  mutate(assembly = sub("([^_]+)_.+", "\\1", seqnames)) %>%
  group_by(assembly) %>%
  count(type) %>%
  ungroup() %>%
  arrange(assembly, desc(n))


frac_by_genome_type <- all_nc_gr[all_nc_gr %over% GRS_gr] %>% as_tibble() %>%
  mutate(assembly = sub("([^_]+)_.+", "\\1", seqnames)) %>%
  group_by(assembly, type) %>%
  summarise(n_type = n(), span_type = sum(width), .groups= "drop") %>%
  arrange(assembly, desc(span_type)) %>% 
  left_join(GRS_size, by = "assembly") %>%
  mutate(frac_rep = round((span_type/span_GRS) * 100, 2))


all_nc_gr %>% as_tibble() %>%
  #   filter(!duplicated(source)) %>%
  #   View 
  # filter(type %in% c("similarity", "tandem_repeat")) %>%
  filter(source %in% c("tRNAScan-SE", "RNAmmer", "RepeatMasker", "TRF")) %>%  
  mutate(assembly = sub("([^_]+)_.+", "\\1", seqnames)) %>%
  group_by(assembly) %>%
  summarise(span_rep = sum(width), .groups= "drop") %>%
  arrange(assembly, desc(span_rep)) %>% 
  left_join(genome_size, by = "assembly") %>%
  mutate(frac_rep = round((span_rep/span) * 100, 2)) %>% 
  # write_tsv("~/Downloads/tmp.tsv")
  View

as_tibble(all_nc_gr) %>%
  mutate(assembly = sub("([^_]+)_.+", "\\1", seqnames),
         repeat_fam = sub("([^/]+)/.*", "\\1", assignment),
         repeat_fam = ifelse(is.na(repeat_fam), source, repeat_fam),
         repeat_fam = ifelse(repeat_fam == "RNAmmer", "rRNA", repeat_fam),
         repeat_fam = ifelse(repeat_fam == "tRNAScan-SE", "tRNA", repeat_fam),
         repeat_fam = ifelse(repeat_fam == "SINE?", "SINE", repeat_fam)) %>%
  group_by(assembly, repeat_fam) %>%
  summarise(span = sum(width), .groups = "drop") %>%
  arrange(assembly, repeat_fam) %>%
  left_join(whole_gnm_sizes, by = "assembly") %>%
  mutate(perc = round((span/(total_span)) * 100, 2)#,
         # diff_perc = perc_GRS - perc_non_GRS
  ) %>%
  select(assembly, repeat_fam, perc) %>%
  pivot_wider(names_from = assembly, values_from = contains("perc"), values_fill = 0, id_expand = T) %>%
  # View()
  write_tsv("~/Downloads/tmp.tsv")



#### Export genes ####



# Export ncRNA and repeats to make a new prediction with softmasked genome
assmbly_col <- seqnames(all_nc_gr) %>% as.character() %>%
  sub("([^_]+)_.+", "\\1", .)
all_nc_gr$Sequence <- seqnames(all_nc_gr) %>% as.character() %>%
  sub("[^_]+_(.+)", "\\1", .)
for (assembly in unique(seq_sizes$assembly)) {
  tmp_gr <- all_nc_gr[assmbly_col == assembly]
  renm_tmp_gr <- GRanges(tmp_gr$Sequence,
                         ranges(tmp_gr, use.mcols=TRUE), strand = strand(tmp_gr))
  export(renm_tmp_gr, paste0("~/Downloads/tmp_annot/", assembly, ".nc_reps.gff2"))
  
}



#### Calculate repeat stats ####

whole_gnm_sizes <- group_by(seq_sizes, assembly) %>%
  summarise(total_span = sum(size), .groups = "drop")

distal_grs_gr <- GRS_gr[GRS_gr %over% GRS_distal_gr]

all_nc_gr$is_GRS <- ifelse(all_nc_gr %over% GRS_gr,
                           ifelse(all_nc_gr %over% GRS_distal_gr, "distal_GRS",
                                  "internal_GRS"), "non_GRS")

# Table of repeats by category

as_tibble(all_nc_gr) %>%
  mutate(assembly = sub("([^_]+)_.+", "\\1", seqnames),
         repeat_fam = sub("([^/]+)/.*", "\\1", assignment),
         repeat_fam = ifelse(is.na(repeat_fam), source, repeat_fam),
         repeat_fam = ifelse(repeat_fam == "RNAmmer", "rRNA", repeat_fam),
         repeat_fam = ifelse(repeat_fam == "tRNAScan-SE", "tRNA", repeat_fam),
         repeat_fam = ifelse(repeat_fam == "SINE?", "SINE", repeat_fam)) %>%
  group_by(assembly, is_GRS, repeat_fam) %>%
  summarise(span = sum(width), .groups = "drop") %>% 
  pivot_wider(names_from = is_GRS, values_from = contains("span"), values_fill = 0, id_expand = T) %>%
  arrange(assembly, repeat_fam) %>%
  bind_rows(GRS_by_type_size) %>%
  select(assembly, repeat_fam, distal_GRS, internal_GRS, non_GRS) %>%
  pivot_wider(names_from = assembly, values_from = contains("GRS"), values_fill = 0, id_expand = T) %>%
  arrange(non_GRS_nxOscDolc1.1) %>%
  filter(repeat_fam != "ARTEFACT") %>%
  write_tsv("~/Downloads/repeats_size_by_GRS.tsv")




as_tibble(all_nc_gr) %>%
  mutate(assembly = sub("([^_]+)_.+", "\\1", seqnames),
         repeat_fam = sub("([^/]+)/.*", "\\1", assignment),
         repeat_fam = ifelse(is.na(repeat_fam), source, repeat_fam),
         repeat_fam = ifelse(repeat_fam == "RNAmmer", "rRNA", repeat_fam),
         repeat_fam = ifelse(repeat_fam == "tRNAScan-SE", "tRNA", repeat_fam),
         repeat_fam = ifelse(repeat_fam == "SINE?", "SINE", repeat_fam)) %>%
  group_by(assembly, is_GRS, repeat_fam) %>%
  summarise(span = sum(width), .groups = "drop") %>% 
  pivot_wider(names_from = is_GRS, values_from = contains("span"), values_fill = 0, id_expand = T) %>%
  arrange(assembly, repeat_fam) %>% 
  left_join(GRS_size, by = "assembly") %>% 
  left_join(whole_gnm_sizes, by = "assembly") %>%
  select(-n_GRS) %>% 
  mutate(perc_GRS = round(((distal_GRS + internal_GRS)/span_GRS) * 100, 2),
         perc_distal_GRS = round((distal_GRS/span_GRS) * 100, 2),
         perc_internal_GRS = round((internal_GRS/span_GRS) * 100, 2),
         perc_non_GRS = round((non_GRS/(total_span - span_GRS)) * 100, 2)#,
         # diff_perc = perc_GRS - perc_non_GRS
  ) %>%
  select(assembly, repeat_fam, perc_GRS, perc_distal_GRS, perc_internal_GRS, perc_non_GRS) %>%
  pivot_wider(names_from = assembly, values_from = contains("perc"), values_fill = 0, id_expand = T) %>% 
  # View()
  write_tsv("~/Downloads/tmp.tsv")

as_tibble(all_nc_gr) %>%
  mutate(assembly = sub("([^_]+)_.+", "\\1", seqnames),
         repeat_fam = sub("([^/]+)/.*", "\\1", assignment),
         repeat_fam = ifelse(is.na(repeat_fam), source, repeat_fam),
         repeat_fam = ifelse(repeat_fam == "RNAmmer", "rRNA", repeat_fam)) %>%
  filter(repeat_fam=="Tandem_repeat", assembly == "nxOscOnir1.2") %>%
  View
all_nc_gr[all_nc_gr$assignment=="Tandem_repeat"]

# What about "SINE?"?
as_tibble(all_nc_gr) %>%
  mutate(assembly = sub("([^_]+)_.+", "\\1", seqnames),
         repeat_fam = sub("([^/]+)/.*", "\\1", assignment),
         repeat_fam = ifelse(is.na(repeat_fam), source, repeat_fam)) %>%
  filter(repeat_fam == "SINE?") %>%
  View


# Percent of repeats for Table 1
GenomicRanges::reduce(trf_gr[!trf_gr %over% rmask_gr]) %>%
  as_tibble() %>%
  mutate(assembly = sub("([^_]+)_.+", "\\1", seqnames)) %>%
  group_by(assembly) %>%
  summarise(tot_TRs = sum(width))

GenomicRanges::reduce(rmask_gr) %>%
  as_tibble() %>%
  mutate(assembly = sub("([^_]+)_.+", "\\1", seqnames)) %>%
  group_by(assembly) %>%
  summarise(tot_TEs = sum(width))







