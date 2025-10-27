library(DESeq2)
library(tidyverse)



metadata <- read_csv("analyses/genes/DEseq2/samplesheet.csv") %>%
  mutate(condition = sub("APS4_(.+)_[1-3]", "\\1", sample_title),
         LifeStage = ifelse(grepl("L2", condition), "L2",
                            ifelse(grepl("MA", condition), "adult",
                                   "mixed")),
         Sex = sub("L2_(.+)", "\\1", condition),
         simpleSex = sub("DA_(.+)", "\\1", Sex))

# Load the deseq2 object
load("analyses/genes/DEseq2/deseq2.dds.RData")
# Convert Sex and LifeStage columns to factors
metadata$Sex <- factor(metadata$Sex)
metadata$LifeStage <- factor(metadata$LifeStage)
# Prepare the sample information for DESeq2
coldata <- metadata[, c("sample", "Sex", "LifeStage")]
rownames(coldata) <- coldata$sample
# Create the DESeqDataSet object
mdds <- DESeqDataSetFromMatrix(countData = round(counts(dds, normalized = TRUE)),
                               colData = coldata,
                               design = ~ Sex)
# Estimate size factors and dispersions
mdds <- DESeq(mdds)

# Test for differential expression
results <- results(mdds, contrast = c("Sex","DA_FEM","FEM"))
write.csv(as.data.frame(results),
          file = "analyses/genes/DEseq2/gene_expression.csv.gz")


# Filter for significant differentially expressed genes
sig_results <- subset(results, padj < 0.05 & abs(log2FoldChange) > 1)

# Write the results to a CSV file
write.csv(as.data.frame(sig_results),
          file = "analyses/genes/DEseq2/DE_induced_hermaphrodites.csv.gz")

sig_results

# Explore gene names
gAnnot <- read_tsv("analyses/genes/nxAuaRhod1_1.annie.tsv",
                   col_names = c("ID", "type", "annot")) %>%
  mutate(ID = sub(".t1", "", ID))

wgAnnot <- filter(gAnnot, type %in% c("name", "product")) %>%
  pivot_wider(names_from = type, values_from = annot)
tdif <- as_tibble(sig_results)
tdif$ID <- rownames(sig_results)
left_join(tdif, wgAnnot, by = "ID") %>% View




