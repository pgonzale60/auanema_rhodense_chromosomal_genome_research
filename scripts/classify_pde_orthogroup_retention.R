#!/usr/bin/env Rscript

#' Classify PDE Orthogroups by Retention Status
#' 
#' This script classifies eliminated genes into categories:
#' - multi copy_in_retained: Paralog exists in the somatic genome.
#' - multi eliminated_only: All gene copies are eliminated from the soma.
#' - single: Unique gene (no paralogs) that is eliminated.
#' 
#' This logic helps determine if PDE removes a unique function or just reduces dosage.

library(tidyverse)

# 1. Inputs
# Result from identify_non_te_eliminated_genes.R
PDE_ALL <- "analyses/genes/PDE_characterization/all_eliminated_genes_annotated.tsv"
ORTHO_GENECOUNT <- "analyses/genes/orthofinder/Orthogroups/Orthogroups.GeneCount.tsv.gz"

# 2. Identification of All Auanema Genes and their Elimination Status
if (!exists("allAuaGenes")) {
    source("scripts/summarize_gene_annotation.R")
}

# 3. Get eliminated genes list
pde_genes <- read_tsv(PDE_ALL)

# 4. Count Eliminated vs Retained per Orthogroup
# First, map all Aua genes to OGs using orthoGenelist (which is available from summarize_gene_annotation.R)
aua_to_og <- orthoGenelist %>%
    filter(assembly == "nxAuaRhod1.1") %>%
    mutate(gene_id = sub(".t[0-9]+", "", transcript_id)) %>%
    select(gene_id, OG_id) %>%
    distinct()

aua_og_map <- allAuaGenes %>%
    left_join(aua_to_og, by = "gene_id") %>%
    filter(!is.na(OG_id))

# Get eliminated IDs for quick lookup
elim_ids <- pde_genes$gene_id

og_stat <- aua_og_map %>%
    group_by(OG_id) %>%
    summarise(
        total_aua = n(),
        elim_aua = sum(gene_id %in% elim_ids),
        retained_aua = total_aua - elim_aua
    ) %>%
    mutate(
        copy_status = ifelse(total_aua == 1, "single", "multi"),
        retention_status = case_when(
            copy_status == "single" ~ "-",
            retained_aua == 0 ~ "eliminated_only",
            retained_aua > 0 ~ "copy_in_retained"
        )
    )

# 5. Apply Classification to PDE Genes
pde_classified <- pde_genes %>%
    left_join(og_stat, by = c("OG_id" = "OG_id")) %>%
    select(gene_id, Preferred_name, Description, OG_id, celegans, TPM, 
           copy_status, retention_status, somatic_copies_retained = retained_aua)

# 6. Write Final Summary
write_tsv(pde_classified, "analyses/genes/PDE_characterization/pde_genes_classified_retention.tsv")

# Print Summary Table
summary_tab <- pde_classified %>%
    filter(!is.na(copy_status)) %>%
    count(copy_status, retention_status)

print("Classification Summary of Eliminated Genes:")
print(summary_tab)
