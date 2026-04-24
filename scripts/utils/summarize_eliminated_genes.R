library(rtracklayer)
library(tidyverse)

# Classification of repeats

### Core tibble creation ####


# eggnog annotation
aua_ipr <- read_tsv("analyses/genes/annotation/nxAuaRhod1_1.iprscan.guoying.tsv.gz",
                    comment = "#",
                    col_names = c("transcript_id", "md5", "slen", "annotator",
                                  "signature_id", "signature_info", "start", "end",
                                  "score", "status", "date", "interpro_id", "interpro_info",
                                  "GO_id", "pathway_id"),
                    col_types = 'cciccciiclccccc') %>%
  mutate(gene_id = sub(".t[0-9]+", "", transcript_id))

# ipro transposon domains
transp_doms <- read_tsv("raw/intreproscan_transposon_related_domains.tsv.gz") %>%
  filter(!is.na(Name))

transp_doms <- filter(aua_ipr, grepl("transpos", interpro_info))


# eggnog annotation
eggm <- read_tsv("analyses/genes/annotation/nxAuaRhod1_1.eggnog_mapper.tsv.gz",
                    comment = "#") %>%
  mutate(gene_id = sub(".t[0-9]+", "", query))

# Protein lengths
prot_lens <- read_tsv("analyses/genes/annotation/nxAuaRhod1_1.longest_isoform.lengths.tsv.gz",
                      col_names = c("tID", "length")) %>%
  mutate(gID = sub(".t[1-9]+", "", tID))

# RNA-seq metadata
metadata <- read_csv("analyses/genes/DEseq2/samplesheet.csv") %>%
  mutate(condition = sub("APS4_(.+)_[1-3]", "\\1", sample_title),
         LifeStage = ifelse(grepl("L2", condition), "L2",
                            ifelse(grepl("MA", condition), "adult",
                                   "mixed")),
         Sex = sub("L2_(.+)", "\\1", condition),
         simpleSex = sub("DA_(.+)", "\\1", Sex),
         name = sub("APS4_", "", sample_title))
# Gene expression
gexp <- read_tsv("analyses/genes/star_stringtie/arhod_star_stringtie_gene_expression.tsv.gz",
                 col_names = c("ID", "Name", "Reference", "Strand", "Start",
                               "End", "Coverage", "FPKM", "TPM", "SRA"))

gexp_max <- group_by(gexp, ID) %>%
  slice_max(TPM, with_ties = F) %>% ungroup() %>%
  left_join(select(metadata, sample, name), by = c("SRA" = "sample"))
rm(gexp)

filter(gexp_max, TPM >= 0.5, grepl("file", ID)) %>%
  nrow()
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

elim_seqs <- filter(seq_sizes, grepl("loc|scaffold", multispecies_sequence)) %>%
  mutate(start = 1, end = size-1) %>%
  select(-size) %>%
  rbind(GRS)


#### GRanges creation ####
# Whole genome sizes
gnm_gr <- Seqinfo(seqnames = seq_sizes$multispecies_sequence,
                  seqlengths=seq_sizes$size,
                  isCircular = rep(F, nrow(seq_sizes)), 
                  genome = "Auanema")

# GRS regions
GRS_gr <- GRanges(seqnames = Rle(elim_seqs$multispecies_sequence), 
                  ranges = IRanges(start = elim_seqs$start,
                                   end = elim_seqs$end), 
                  strand = "*")

genes_gr <- makeGRangesFromDataFrame(genes_df, keep.extra.columns=T,
                                     seqinfo=gnm_gr)

elim_genes_annot <- genes_gr[genes_gr %over% GRS_gr & genes_gr$type == "mRNA"] %>%
  as_tibble() %>% mutate(tID=sub("ID=([^;]+);.+", "\\1", ID),
                         Dbxref = ifelse(grepl("Dbxref=", ID),
                                         sub("ID=.+Dbxref=([^;]+).*", "\\1", ID),
                                         NA))

annot_gIDs <- c(eggm$gene_id, aua_ipr$gene_id)
transp_gIDs <- filter(aua_ipr, signature_id %in% transp_doms$Accession) %>%
  pull(gene_id) %>% unique()


elim_g_table <- left_join(elim_genes_annot, eggm, by = c("tID" = "query")) %>%
  mutate(seqnames = sub("nx[^_]+_1_", "", seqnames)) %>%
  left_join(prot_lens, by = c("tID")) %>%
  rename(aa_len = length) %>%
  left_join(select(gexp_max, gID = ID, TPM, name), by = c("gID")) %>%
  mutate(tID = sub("file_1_file_1_", "", tID),
         Preferred_name = ifelse(Preferred_name == "-", NA, Preferred_name),
         Description = ifelse(Description == "-", NA, Description),
         TPM = round(TPM, 1),
         transposon_rel_dom = ifelse(gene_id %in% transp_gIDs |
                                       grepl("Ribonuclease H|Reverse transcriptase|K02A2.6-like|transposition|ntegrase",
                                             Description), 
                                     "yes", "no")) %>%
  arrange(desc(TPM), seqnames, start) %>%
  select(seqnames:end, tID, aa_len, Preferred_name, Description, transposon_rel_dom, TPM, name, Dbxref) 

elim_g_table %>%
  write_tsv("report/tables/eliminated_genes.tsv")

nrow(elim_g_table)
count(elim_g_table, TPM)
count(elim_g_table, name)
unique(elim_g_table)

sp_tids <- filter(elim_g_table, grepl("Spa2", Description)) %>%
  pull(tID)

filter(aua_ipr, transcript_id %in% paste0("file_1_file_1_", sp_tids),
       signature_id %in% transp_doms$Accession) %>%
  View
  

# Tandem repeats lengths
tr_lens <- read_tsv()
