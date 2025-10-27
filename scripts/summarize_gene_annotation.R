library(tidyverse)

aua_ipr <- read_tsv("analyses/genes/annotation/nxAuaRhod1_1.iprscan.guoying.tsv.gz",
                    comment = "#",
                    col_names = c("transcript_id", "md5", "slen", "annotator",
                                  "signature_id", "signature_info", "start", "end",
                                  "score", "status", "date", "interpro_id", "interpro_info",
                                  "GO_id", "pathway_id"),
                    col_types = 'cciccciiclccccc') %>%
  mutate(gene_id = sub(".t[0-9]+", "", transcript_id))


# eggnog annotation
aua_egg <- read_tsv("analyses/genes/annotation/nxAuaRhod1_1.eggnog_mapper.tsv.gz",
                    comment = "#") %>%
  mutate(gene_id = sub(".t[0-9]+", "", query))

# ipro transposon domains
transp_doms <- read_tsv("raw/intreproscan_transposon_related_domains.tsv.gz")


# Gene expression
gexp <- read_tsv("analyses/genes/star_stringtie/arhod_star_stringtie_gene_expression.tsv.gz",
                 col_names = c("ID", "Name", "Reference", "Strand", "Start",
                               "End", "Coverage", "FPKM", "TPM", "SRA"))

gexp_max <- group_by(gexp, ID) %>%
  slice_max(TPM, with_ties = F) %>% ungroup() %>%
  left_join(select(metadata, sample, name), by = c("SRA" = "sample"))
rm(gexp)


ortoCounts <- read_tsv("analyses/genes/orthofinder/Orthogroups/Orthogroups.GeneCount.tsv.gz")
unassigned <- read_tsv("analyses/genes/orthofinder/Orthogroups/Orthogroups_UnassignedGenes.tsv.gz",
                       guess_max = 1e6)
orthoGenes <- read_tsv("analyses/genes/orthofinder/Orthogroups/Orthogroups.tsv.gz") %>%
  bind_rows(unassigned)

orthoGenelist <- tibble(transcript_id = character(),
                        OG_id = character(),
                        assembly = character())
for (sp in colnames(orthoGenes)[-1]) {
  og <- orthoGenes %>%
    filter(!is.na(.data[[sp]]))
  orthoGs <- map2(pull(og, .data[[sp]]), og$Orthogroup, function(x, y){
    spt_gs <- as.character(unlist(str_split(x, ", ")))
    data.frame(transcript_id = spt_gs,
               OG_id = rep(y, length(spt_gs)),
               stringsAsFactors=FALSE)
  }) %>%
    bind_rows() %>%
    add_column(assembly = sp)
  orthoGenelist <- bind_rows(orthoGenelist, orthoGs)
}


# Filters ####
withExpr <- filter(gexp_max, TPM >= 1, grepl("file", ID)) #%>%
  # nrow()
annot_gIDs <- c(aua_egg$gene_id, aua_ipr$gene_id)
transp_gIDs <- filter(aua_ipr, signature_id %in% transp_doms$Accession) %>%
  pull(gene_id) %>% unique()

withOrts <- filter(ortoCounts, nxAuaRhod1.1 > 0) %>%
  # select(-c(Orthogroup, nxAuaRhod1.1, Total)) %>%
  select(-c(Orthogroup, nxAuaRhod1.1, arabidopsis_thaliana,
            saccharomyces_cerevisiae, Total)) %>%
  rowSums() > 1 

auaOGwithOrts <- filter(ortoCounts, nxAuaRhod1.1 > 0) %>%
  filter(withOrts) %>% pull(Orthogroup) 
  # summarise(Auawithorts = sum(nxAuaRhod1.1))
auaOGwithOrts <- filter(orthoGenelist, OG_id %in% auaOGwithOrts,
                        assembly == "nxAuaRhod1.1") %>%
  mutate(gene_id = sub(".t[0-9]+", "", transcript_id))

allAuaGenes <- tibble(gene_id = gexp_max$ID[grepl("file", gexp_max$ID)]) %>%
  mutate(withOrts = gene_id %in% auaOGwithOrts$gene_id,
         exprssed = gene_id %in% withExpr$ID,
         withAnnot = gene_id %in% annot_gIDs,
         fullEv = withOrts & exprssed & withAnnot,
         transp_rel = gene_id %in% transp_gIDs,
         expTransp = transp_rel & exprssed)

colSums(select(allAuaGenes, -gene_id))


# Size of chromosomal sequences
filter(seq_sizes, grepl("SUPER", Sequence)) %>%
  summarise(sum(size))
filter(seq_sizes, !grepl("unloc|scaf|MT", Sequence)) %>%
  summarise(sum(size))
filter(seq_sizes, grepl("MT", Sequence)) %>%
  summarise(sum(size))

## ncRNAs ####
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
  mutate(seqnames = Sequence,
         gene_id = as.integer(sub(".+Parent=[^0-9]+(.+)_gen.*", "\\1", ID))) %>%
  filter(!grepl("MT|scaf|unloc", Sequence)) %>%
  select(-assembly, -Sequence, -phase)

utrna_df <- filter(trna_df, !duplicated(gene_id))

# infernal files
miscrnafiles <- list.files("analyses/genome_features/repeats/infernal/",
                           full.names = T, pattern = "gff.gz")
names(miscrnafiles) <- make.names(sub(".+//(.+).infernal.gff.gz", "\\1", miscrnafiles))
miscrna_df <- map_df(miscrnafiles, read_tsv,
                     .id = "assembly", col_names = c("Sequence", "source", "type", "start", "end", "score",
                                                     "strand", "phase", "ID")) %>% 
  filter(!grepl("rRNA|tRNA|Protozoa", type),
         !type %in% c("K_chan_RES", "RNase_MRP", "Fluoride", "GlsR7", "RAGATH-21", "snosnR60_Z15",
                      "Histone3")) %>%
  mutate(seqnames = Sequence) %>%
  filter(!grepl("MT|scaf|unloc", Sequence))

nrow(miscrna_df)
count(miscrna_df, type) %>% arrange(desc(n)) %>% View()
nrow(rrna_df) + nrow(utrna_df) + nrow(miscrna_df)
filter(utrna_df, grepl("pseudo", ID)) %>% nrow

## Genome coverage #####
covFile <- "analyses/genome_features/nemaChromQC/nxAuaRhod1_1.regions.bed.gz"
cov <- read_tsv(covFile,
                col_names = c("seqnames", "start",
                              "end", "median_cov"))%>%
  filter(!grepl("MT", seqnames)) %>%
  mutate(start = start +1)

# Sequence sizes
fai_files <- list.files("analyses/genome_features/sequence_sizes/",
                        full.names = T, pattern = ".primary.fa.gz.fai")
names(fai_files) <- make.names(sub(".+//(.+).primary.fa.gz.fai", "\\1", fai_files))
seq_sizes <- map_df(fai_files, read_tsv,
                    col_names = c("Sequence", "size", "L1", "L2", "L3"),
                    col_types = c("ciiii"),
                    .id = "assembly") %>%
  select(assembly, Sequence, size) %>%
  filter(!grepl("MT", Sequence))

# GRS files
GRS_files <- list.files("analyses/genome_features/elim_coords/",
                        full.names = T, pattern = "nx.*.GRS.bed$")
names(GRS_files) <- make.names(sub(".+//(.+).GRS.bed", "\\1", GRS_files))
GRS <- map_df(GRS_files, read_tsv,
              col_names = c("Sequence", "start", "end", "has_GRS"),
              col_types = c("ciic"),
              .id = "assembly")


elim_seqs <- filter(seq_sizes, grepl("loc|scaffold", Sequence)) %>%
  mutate(start = 1, end = size-1) %>%
  select(-size) %>%
  rbind(GRS)


#### GRanges creation ####
# Whole genome sizes
gnm_gr <- Seqinfo(seqnames = seq_sizes$Sequence,
                  seqlengths=seq_sizes$size,
                  isCircular = rep(F, nrow(seq_sizes)), 
                  genome = "Auanema")

# GRS regions
GRS_gr <- GRanges(seqnames = Rle(elim_seqs$Sequence), 
                  ranges = IRanges(start = elim_seqs$start,
                                   end = elim_seqs$end), 
                  strand = "*")

cov_gr <- makeGRangesFromDataFrame(cov, keep.extra.columns=T,
                                     seqinfo=gnm_gr)

cov_gr[cov_gr %over% GRS_gr] %>%
  as_tibble() %>% summary()

cov_gr[!cov_gr %over% GRS_gr] %>%
  as_tibble() %>% summary()

cov_gr[!cov_gr %over% GRS_gr] %>%
  as_tibble() %>% filter(seqnames == "SUPER_X") %>%
  summary()

cov_gr[!cov_gr %over% GRS_gr] %>%
  as_tibble() %>% filter(seqnames != "SUPER_X") %>%
  summary()
  
cov_gr[cov_gr %over% GRS_gr] %>%
  as_tibble() %>% filter(seqnames == "SUPER_X") %>%
  summary()

