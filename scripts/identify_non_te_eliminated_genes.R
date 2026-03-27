#!/usr/bin/env Rscript

#' Identify Non-TE Eliminated Genes in Auanema rhodense
#' 
#' This script identifies the "legitimate" (non-TE-derived) protein-coding genes 
#' that are eliminated during Programmed DNA Elimination (PDE) in Auanema rhodense.
#' 
#' Logic:
#' 1. Identify genes in GRS regions or small scaffolds (excluding Mitochondria).
#' 2. Filter out genes with >50% physical overlap with TE models (Earlgrey).
#' 3. Filter out genes with TE-related domains (transposase, gag, pol, etc.) unless 
#'    they have confirmed cross-species orthology in OrthoFinder.

library(tidyverse)
library(rtracklayer)

# 1. Input Paths
TE_GFF <- "analyses/genome_features/repeats/earlgrey/Auanema_rhodensis.filteredRepeats.gff"
GRS_BED <- "analyses/genome_features/elim_coords/nxAuaRhod1_1.GRS.bed"
GENE_GFF <- "analyses/genes/GTFs/nxAuaRhod1_1.gff.gz"
EPP_ANNOT <- "analyses/genes/annotation/nxAuaRhod1_1.eggnog_mapper.tsv.gz"
IPR_ANNOT <- "analyses/genes/annotation/nxAuaRhod1_1.iprscan.guoying.tsv.gz"
ORTHOGROUPS <- "analyses/genes/orthofinder/Orthogroups/Orthogroups.GeneCount.tsv.gz"
ORTHO_GENELIST <- "analyses/genes/orthofinder/Orthogroups/Orthogroups.tsv.gz"
EXPRESSION <- "analyses/genes/expression/max_tpm_per_gene.tsv"

# 2. Source the base annotation info (reusing existing workspace logic)
if (!exists("genes_gr")) {
    source("scripts/summarize_gene_annotation.R")
}

# 3. Define Elimination Criteria
GRS_gr <- import(GRS_BED)
elim_genes_gr <- genes_gr[(genes_gr %over% GRS_gr | grepl("unloc|scaf", seqnames(genes_gr))) & 
                          !grepl("MT", seqnames(genes_gr))]
elim_gene_ids <- elim_genes_gr$gene_id

# 4. Filter for Non-TE Eliminated Genes
# TE keyword search (using word boundaries to avoid false positives like polypeptide)
te_keywords <- c("\\btransposon\\b", "\\btransposase\\b", "\\bretrotransposon\\b", "\\breverse transcriptase\\b", 
                 "\\bintegrase\\b", "\\bgag\\b", "\\bpol\\b", "\\bviropepsin\\b", "\\brvt\\b", "\\brnase h\\b", 
                 "\\baspartic peptidase\\b", "\\bpao\\b", "\\bgypsy\\b", "\\bcopia\\b", "\\bbel\\b", 
                 "\\bdirs\\b", "\\btc1\\b", "\\bmariner\\b")

# Find genes with matching domains in IPR
te_ipr_gids <- aua_ipr %>%
    filter(grepl(paste(te_keywords, collapse="|"), signature_info, ignore.case = TRUE) |
           grepl(paste(te_keywords, collapse="|"), interpro_info, ignore.case = TRUE)) %>%
    pull(gene_id) %>% unique()

# Find genes with matching names in EggNOG
te_egg_gids <- aua_egg %>%
    filter(grepl(paste(te_keywords, collapse="|"), Preferred_name, ignore.case = TRUE) |
           grepl(paste(te_keywords, collapse="|"), Description, ignore.case = TRUE)) %>%
    pull(gene_id) %>% unique()

te_domain_gids <- unique(c(te_ipr_gids, te_egg_gids))

# 5. Integrate OrthoFinder data
ortoCounts <- read_tsv(ORTHOGROUPS)
multi_species_ogs <- ortoCounts %>%
    filter(Total > nxAuaRhod1.1) %>%
    pull(Orthogroup)

# Map Aua genes to OrthoFinder OGs using orthoGenelist
aua_to_og <- orthoGenelist %>%
    filter(assembly == "nxAuaRhod1.1") %>%
    mutate(gene_id = sub(".t[0-9]+", "", transcript_id)) %>%
    select(gene_id, OG_id) %>%
    distinct()

# Final Identification Logic
# We keep ALL eliminated genes but tag them for potential TE origin
elim_all_annot <- allAuaGenes %>%
    filter(gene_id %in% elim_gene_ids) %>%
    left_join(aua_to_og, by = "gene_id") %>%
    mutate(
        is_legit_homolog = !is.na(OG_id) & OG_id %in% multi_species_ogs,
        has_te_domains = gene_id %in% te_domain_gids,
        is_non_te_candidate = !te_contained & (is_legit_homolog | !has_te_domains)
    )

# 6. Final Annotation Join
orthoGroups <- read_tsv(ORTHO_GENELIST) %>%
    select(Orthogroup, celegans = caenorhabditis_elegans.PRJNA13758)

elim_all_annot <- elim_all_annot %>%
    left_join(select(gexp_max, ID, TPM), by = c("gene_id" = "ID")) %>%
    left_join(select(aua_egg, gene_id, Preferred_name, Description), by = "gene_id") %>%
    left_join(orthoGroups, by = c("OG_id" = "Orthogroup")) %>%
    select(gene_id, Preferred_name, Description, OG_id, celegans, TPM, 
           te_contained, has_te_domains, is_legit_homolog, is_non_te_candidate) %>%
    arrange(desc(is_non_te_candidate), desc(TPM))

# 7. Write Results
dir.create("analyses/genes/PDE_characterization", showWarnings = FALSE)
write_tsv(elim_all_annot, "analyses/genes/PDE_characterization/all_eliminated_genes_annotated.tsv")
# High-confidence non-TE subset
write_tsv(filter(elim_all_annot, is_legit_homolog, !has_te_domains), 
          "analyses/genes/PDE_characterization/legitimate_non_te_elim_candidates.tsv")

cat("Total eliminated genes (nuclear):", nrow(elim_all_annot), "\n")
cat("Identified", sum(elim_all_annot$is_non_te_candidate), "non-TE candidates.\n")
