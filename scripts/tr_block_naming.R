#!/usr/bin/env Rscript
# tr_block_naming.R
#
# Identifies the dominant tandem repeat (TR) family in each of the five
# genomic blocks per chromosome of Auanema rhodensis:
#
#   Block 1 – Left eliminated    (GRS 1, terminal)
#   Block 2 – Left somatic       (between GRS 1 and GRS 2)
#   Block 3 – Central eliminated (GRS 2, internal)
#   Block 4 – Right somatic      (between GRS 2 and GRS 3)
#   Block 5 – Right eliminated   (GRS 3, terminal)
#
# Output: TSV with top TR families per block, with systematic names applied.
#
# Usage (from project root):
#   Rscript scripts/tr_block_naming.R

library(tidyverse)

# ── Parameters ─────────────────────────────────────────────────────────────────
FAI_FILE     <- "analyses/genome_features/sequence_sizes/nxAuaRhod1.1.primary.fa.gz.fai"
GRS_FILE     <- "analyses/genome_features/elim_coords/nxAuaRhod1_1.GRS.bed"
MEMBERS_FILE <- "analyses/genome_features/repeats/TRF/trf_family_members.tsv"
OUTPUT_TSV   <- "analyses/genome_features/repeats/TRF/tr_block_naming.tsv"
TOP_N        <- 3   # number of top families to report per block

# ── TR naming scheme ───────────────────────────────────────────────────────────
# Maps internal family IDs to systematic names.
# Autosomal terminal repeats: TR1 (dominant), TR2 (paired companion)
# Autosomal central repeats:  AR-I to AR-VI (chromosome-specific, ordered by
#                             chromosome number)
# X chromosome eliminated:    TR25 (X-left terminal), TR26 (X-right terminal),
#                             XR-central (X central GRS, minor)
# X chromosome somatic:       TR30 (left somatic domain), TR30b (right somatic
#                             domain)
TR_NAME_MAP <- c(
  # Autosomal terminal
  "P351_F001" = "TR1",
  "P176_F001" = "TR2",
  # Autosomal central (one dominant family per autosome)
  "P348_F001" = "AR-I",    # Chr 1
  "P399_F001" = "AR-II",   # Chr 2
  "P167_F002" = "AR-III",  # Chr 3
  "P332_F002" = "AR-IV",   # Chr 4
  "P231_F004" = "AR-V",    # Chr 5
  "P412_F002" = "AR-VI",   # Chr 6
  # X chromosome
  "P347_F001" = "TR25",       # X left terminal (eliminated)
  "P342_F001" = "TR26",       # X right terminal (eliminated)
  "P248_F001" = "XR-central", # X central GRS (minor, eliminated)
  "P348_F003" = "TR30",       # X left somatic domain (retained)
  "P291_F001" = "TR30b"       # X right somatic domain (retained)
)

# ── Load data ──────────────────────────────────────────────────────────────────
message("Loading genome info...")
seq_sizes <- read_tsv(FAI_FILE,
    col_names = c("Chrom", "size", "L1", "L2", "L3"), col_types = "ciiii") %>%
  select(Chrom, size) %>%
  filter(!grepl("MT|unloc|scaffold", Chrom))

message("Loading GRS regions...")
grs <- read_tsv(GRS_FILE, col_names = c("Chrom", "Start", "End"),
                col_types = "cii") %>%
  filter(Chrom %in% seq_sizes$Chrom)

message("Loading TR family members...")
members <- read_tsv(MEMBERS_FILE, col_types = cols(.default = "c",
    Start = "i", End = "i", Length = "i")) %>%
  filter(Chrom %in% seq_sizes$Chrom) %>%
  select(Family_ID, Chrom, Start, End)

# ── Build 5-block structure per chromosome ─────────────────────────────────────
# Each chromosome has exactly 3 GRS regions ordered left to right.
# Blocks 1, 3, 5 are eliminated (GRS); blocks 2 and 4 are somatic (gaps).
message("Building 5-block structure...")

grs_split <- grs %>%
  group_by(Chrom) %>%
  mutate(grs_idx = row_number()) %>%
  ungroup()

grs1 <- grs_split %>% filter(grs_idx == 1) %>% rename(g1s = Start, g1e = End) %>% select(-grs_idx)
grs2 <- grs_split %>% filter(grs_idx == 2) %>% rename(g2s = Start, g2e = End) %>% select(-grs_idx)
grs3 <- grs_split %>% filter(grs_idx == 3) %>% rename(g3s = Start, g3e = End) %>% select(-grs_idx)

grs_wide <- grs1 %>%
  left_join(grs2, by = "Chrom") %>%
  left_join(grs3, by = "Chrom")

blocks <- bind_rows(
  grs_wide %>% transmute(Chrom, Block = 1L, Block_Type = "eliminated", Block_Start = g1s, Block_End = g1e),
  grs_wide %>% transmute(Chrom, Block = 2L, Block_Type = "somatic",    Block_Start = g1e + 1L, Block_End = g2s - 1L),
  grs_wide %>% transmute(Chrom, Block = 3L, Block_Type = "eliminated", Block_Start = g2s, Block_End = g2e),
  grs_wide %>% transmute(Chrom, Block = 4L, Block_Type = "somatic",    Block_Start = g2e + 1L, Block_End = g3s - 1L),
  grs_wide %>% transmute(Chrom, Block = 5L, Block_Type = "eliminated", Block_Start = g3s, Block_End = g3e)
) %>%
  mutate(Block_Size = Block_End - Block_Start + 1L) %>%
  arrange(Chrom, Block)

# ── Compute TR coverage per block ──────────────────────────────────────────────
message("Computing TR coverage per block (may take a moment)...")

# Cross-join blocks with members, then filter to overlapping pairs
block_coverage <- inner_join(blocks, members, by = "Chrom",
                             relationship = "many-to-many") %>%
  filter(Start < Block_End, End > Block_Start) %>%
  mutate(Overlap = pmin(End, Block_End) - pmax(Start, Block_Start)) %>%
  group_by(Chrom, Block, Block_Type, Block_Start, Block_End, Block_Size, Family_ID) %>%
  summarise(Coverage = sum(Overlap), .groups = "drop") %>%
  mutate(Coverage_Fraction = Coverage / Block_Size)

# ── Select top N families per block and apply naming scheme ───────────────────
message("Applying TR naming scheme...")

top_families <- block_coverage %>%
  group_by(Chrom, Block) %>%
  slice_max(Coverage, n = TOP_N, with_ties = FALSE) %>%
  mutate(Rank = row_number()) %>%
  ungroup() %>%
  mutate(
    Chr_Name = sub("SUPER_", "Chr ", Chrom),
    TR_Name  = dplyr::recode(Family_ID, !!!TR_NAME_MAP, .default = Family_ID)
  ) %>%
  select(Chrom, Chr_Name, Block, Block_Type, Block_Start, Block_End,
         Block_Size, Rank, Family_ID, TR_Name, Coverage, Coverage_Fraction) %>%
  arrange(Chrom, Block, Rank)

# ── Write output ───────────────────────────────────────────────────────────────
write_tsv(top_families, OUTPUT_TSV)
message("Written: ", OUTPUT_TSV)

# ── Print summary (rank-1 only) ────────────────────────────────────────────────
cat("\n── Dominant TR per block ─────────────────────────────────────────────────\n")
top_families %>%
  filter(Rank == 1) %>%
  mutate(
    Coords   = sprintf("%s:%d-%d", Chrom, Block_Start, Block_End),
    Coverage = sprintf("%.1f Mb (%.0f%%)", Coverage / 1e6, Coverage_Fraction * 100)
  ) %>%
  select(Chr_Name, Block, Block_Type, Coords, TR_Name, Family_ID, Coverage) %>%
  print(n = Inf, width = 120)
