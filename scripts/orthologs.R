library(tidyverse)

aua_ipr <- read_tsv("analyses/genes/annotation/nxAuaRhod1_1.iprscan.guoying.tsv.gz",
                    comment = "#",
                    col_names = c("transcript_id", "md5", "slen", "annotator",
                                  "signature_id", "signature_info", "start", "end",
                                  "score", "status", "date", "interpro_id", "interpro_info",
                                  "GO_id", "pathway_id"),
                    col_types = 'cciccciiclccccc')
# eggnog annotation
aua_egg <- read_tsv("analyses/genes/annotation/nxAuaRhod1_1.eggnog_mapper.tsv.gz",
                    comment = "#")



annot_gIDs <- c(aua_egg$query, aua_ipr$transcript_id) %>%
  sub(".t[0-9]+", "", .)


length(unique(aua_ipr$transcript_id))
length(unique(annot_gIDs))

put_ipr_kin <- filter(aua_ipr, grepl("centromere|kinetochore|CENP", signature_info, ignore.case = T))
put_ipr_zinc <- filter(aua_ipr, grepl("Zinc finger|Zinc-finger", signature_info, ignore.case = T))



put_egg_kin <- filter(aua_egg, grepl("0000775", GOs))

any_protein_bind <- filter(aua_egg, grepl("0003676", GOs))


table(put_ipr_kin$transcript_id %in% put_egg_kin$query)

humankin <- read_tsv("raw/human_kinetochora_proteins.tsv.gz") %>%
  rename(unip_name = `Entry Name`)
celegkin <- read_tsv("raw/celeg_kinetochora_proteins.tsv.gz") %>%
  mutate(pID = sub(".+ ", "", `Gene Names`))



humanorth <- read_tsv("analyses/genes/orthofinder/orthologues/nxAuaRhod1.1__v__homo_sapiens.tsv") 
celegorth <- read_tsv("analyses/genes/orthofinder/orthologues/nxAuaRhod1.1__v__caenorhabditis_elegans.PRJNA13758.tsv")

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


withOrts <- filter(ortoCounts, nxAuaRhod1.1 > 0) %>%
  # select(-c(Orthogroup, nxAuaRhod1.1, Total)) %>%
  select(-c(Orthogroup, nxAuaRhod1.1, arabidopsis_thaliana,
            saccharomyces_cerevisiae, homo_sapiens, Total)) %>%
  rowSums() > 1 

filter(ortoCounts, nxAuaRhod1.1 > 0) %>%
  filter(withOrts) %>%
  summarise(Auawithorts = sum(nxAuaRhod1.1))


kineg_ortg <- filter(orthoGenelist,
                     transcript_id %in%
                       c(put_egg_kin$query, put_ipr_kin$transcript_id)) %>%
  pull(OG_id) %>% unique()

dnabind_ortg <- filter(orthoGenelist,
                     transcript_id %in%
                       c(any_protein_bind$query)) %>%
  pull(OG_id) %>% unique()

znf_ortg <- filter(orthoGenelist,
                       transcript_id %in%
                         c(put_ipr_zinc$transcript_id)) %>%
  pull(OG_id) %>% unique()



filter(orthoGenes, Orthogroup %in% znf_ortg) %>% View

#
filter(orthoGenes, Orthogroup %in% znf_ortg,
       is.na(caenorhabditis_elegans.PRJNA13758)) %>% View

# file_1_file_1_g7838 CTCFL
# file_1_file_1_g10230






orthoGenelist <- tibble(transcript_id = character(),
                        OG_id = character(),
                        assembly = character())
for (sp in c( "nxAuaRhod1.1")) { # "homo_sapiens" "caenorhabditis_elegans.PRJNA13758
  og <- celegorth %>%
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


table(put_egg_kin$query %in% orthoGenelist$transcript_id)
filter()

tmp <- mutate(orthoGenelist,
              pID = sub("\\.+[0-9]+", "", transcript_id)) %>%
  left_join(celegkin, by = "pID")

tmp <- mutate(orthoGenelist,
                        Entry = sub(".+\\|([^\\|]+)\\|.+", "\\1", transcript_id)) %>%
  left_join(humankin, by = "Entry")

table(is.na(tmp$Reviewed))


orthoGenelist <- mutate(orthoGenelist,
       unip_name = sub(".+\\|([^\\|]+)", "\\1", transcript_id)) %>%
  left_join(humankin, by = "unip_name")








  
