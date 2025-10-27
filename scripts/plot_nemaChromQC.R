#!/usr/bin/env Rscript

library(ggpubr)
library(gtools)
library(scales)
library(tidyverse)


assemName               <- "nxAuaRhod1_1"
nigonDictFile           <- "~/Programs/vis_ALG/data/gene2Nigon_busco20200927.tsv.gz"
buscoFile               <- "analyses/genome_features/nemaChromQC/busco/nxAuaRhod1_1.curated_primary_nematoda_odb10_full_table.tsv"
teloMappedFile          <- "analyses/genome_features/nemaChromQC/teloMaps/nxAuaRhod1_1.curated_primary.teloMapped.paf.gz"
teloRepsFile            <- "analyses/genome_features/nemaChromQC/teloRepeatCounts/nxAuaRhod1_1.curated_primary_teloRepeatCounts.tsv.gz"
allRepsFile             <- "analyses/genome_features/nemaChromQC/red/nxAuaRhod1_1.curated_primary.red.bed.gz"
gcFile                  <- "analyses/genome_features/nemaChromQC/gc/nxAuaRhod1_1.curated_primary.gc.bed.gz"
covFile                 <- "analyses/genome_features/nemaChromQC/nxAuaRhod1_1.regions.bed.gz"
minimumGenesPerSequence <- 30
minNigonFrac            <- 0.2
minFracAlignedTeloReads <- 0.1
windwSize               <- 5e5

# Parameters used when testing
# assemName <- "DF5120.canu.purged.hic_scaff"
# nigonDictFile <- "gene2Nigon_busco20200927.tsv.gz"
# buscoFile <- "DF5120.canu.purged.hic_scaff_nematoda_odb10_full_table.tsv"
# teloMappedFile <- "DF5120.canu.purged.hic_scaff.teloMapped.paf.gz"
# teloRepsFile <- "DF5120.canu.purged.hic_scaff_teloRepeatCounts.tsv.gz"
# allRepsFile <- "DF5120.canu.purged.hic_scaff.red.bed.gz"
# gcFile <- "DF5120.canu.purged.hic_scaff.gc.bed.gz"
# minimumGenesPerSequence <- 15
# minNigonFrac <- .90
# minFracAlignedTeloReads <- 0.1
# windwSize <- 5e5

cols <- c("A" = "#af0e2b", "B" = "#e4501e",
          "C" = "#4caae5", "D" = "#f3ac2c",
          "E" = "#57b741", "N" = "#8880be",
          "X" = "#81008b", "-" = "#aaaaaa")


# functions to read PAF copied from github.com/thackl/thacklr/R/read.R
read_paf <- function (file, max_tags = 20){
  col_names <- c("query_name", "query_length", "query_start", 
                 "query_end", "strand", "target_name", "target_length", 
                 "target_start", "target_end", "map_match", "map_length", 
                 "map_quality")
  col_types <- "ciiicciiiiin"
  
  if(max_tags > 0){
    col_names <- c(col_names, paste0("tag_", seq_len(max_tags)))
    col_types <- paste0(col_types, paste(rep("?", max_tags), collapse=""))
  }
  
  read_tsv(file, col_names = col_names, col_types = col_types) %>%
    tidy_paf_tags
}

tidy_paf_tags <- function(.data){
  tag_df <- tibble(.rows=nrow(.data))
  tag_types <- c()
  seen_empty_tag_col <- FALSE
  
  for (x in select(.data, starts_with("tag_"))){
    tag_mx <- str_split(x, ":", 3, simplify=T)
    tag_mx_nr <- na.omit(unique(tag_mx[,1:2]))
    if(nrow(tag_mx_nr) == 0){
      seen_empty_tag_col <- TRUE
      break; # empty col -> seen all tags
    }
    tags <- tag_mx_nr[,1]
    tag_type <- tag_mx_nr[,2]
    names(tag_type) <- tags
    # add to global tag_type vec
    tag_types <- c(tag_types, tag_type)
    tag_types <- tag_types[unique(names(tag_types))]
    # sort tag values into tidy tag tibble
    for (tag in tags){
      if(!has_name(tag_df, tag)){ # init tag
        tag_df[[tag]] <- NA
      }
      tag_idx <- tag_mx[,1] %in% tag
      tag_df[[tag]][tag_idx] <- tag_mx[tag_idx,3]
    }
  }
  
  tag_df <- tag_df %>%
    mutate_at(names(tag_types)[tag_types == "i"], as.integer) %>%
    mutate_at(names(tag_types)[tag_types == "f"], as.numeric)
  
  if(!seen_empty_tag_col)
    rlang::warn("Found tags in max_tags column, you should increase max_tags to ensure all tags for all entries got included")
  
  rlang::inform(
    str_glue("Read and tidied up a .paf file with {n_tags} optional tag fields:\n{s_tags}", s_tags = toString(names(tag_types)), n_tags = length(tag_types)))
  
  bind_cols(select(.data, -starts_with("tag_")), tag_df)
}

# Function to get mode 
getMode <- function(x) {
  keys <- unique(x)
  keys[which.max(tabulate(match(x, keys)))]
}

# Function to group reads in blocks
block_mappings <- function(teloMappings){
  
  hqTeloMappings <- filter(teloMappings,
                           tp == "P",
                           map_length > query_length * 0.8,
                           target_length > windwSize * 2)
  duplicateIds <- filter(hqTeloMappings, duplicated(query_name)) %>%
    pull(query_name) %>% unique
  unduplicated <- filter(hqTeloMappings, query_name %in% duplicateIds) %>%
    group_by(query_name) %>%
    mutate(onlyInFrag = max(target_length) < windwSize * 2) %>%
    filter((target_length < windwSize * 2 & onlyInFrag) | target_length >= windwSize * 2) %>%
    slice_min(query_start) %>%
    slice_max(map_quality) %>%
    slice_min(target_start) %>%
    ungroup()
  onlyOneMap <- filter(hqTeloMappings, !query_name %in% duplicateIds) %>%
    bind_rows(unduplicated)
  
  
  
  teloBlocks <- mutate(onlyOneMap, rstart = ifelse(strand == "+",
                                                   target_start, target_end)) %>%
    arrange(target_name, rstart) %>%
    group_by(target_name, strand) %>%
    mutate(block = cumsum(c(1, diff(rstart) > 100))) %>% #filter(target_name == "ptg000011l") %>% select(target_name, strand, target_start, target_end, block) %>% View
    group_by(target_name, strand, block) %>%
    summarise(strand = unique(strand),
              regionStart = ifelse(strand == "-", getMode(target_start),
                                   getMode(target_start)),
              regionEnd = ifelse(strand == "+", getMode(target_end),
                                 getMode(target_end)),
              teloPos = ifelse(strand == "+", regionStart, regionEnd) ,
              target_length = unique(target_length),
              regSupport = n(),
              .groups = "drop") %>%
    arrange(target_name, teloPos) %>%
    select(target_name, target_length, strand, block, teloPos, regSupport)
  return(teloBlocks)
}




# Load data
nigonDict <- read_tsv(nigonDictFile,
                      col_types = c(col_character(), col_character()))
busco <- suppressWarnings(read_tsv(buscoFile,
                                   col_names = c("Busco_id", "Status", "Sequence",
                                                 "start", "end", "strand", "Score", "Length",
                                                 "OrthoDB_url", "Description"),
                                   col_types = c("ccciicdicc"),
                                   comment = "#"))

teloMappings <- read_paf(teloMappedFile) %>%
  mutate(target_name = paste0("chr", sub("SUPER_", "", target_name)))
  
telomWind <- read_tsv(teloRepsFile,
                      col_names = c("contig", "wStart",
                                    "wEnd", "value",
                                    "feat")) %>%
  mutate(contig = paste0("chr", sub("SUPER_", "", contig)))
repsWind <- read_tsv(allRepsFile,
                     col_names = c("contig", "wStart",
                                   "wEnd", "value",
                                   "feat")) %>%
  mutate(contig = paste0("chr", sub("SUPER_", "", contig)))

gcWind <- read_tsv(gcFile,
                   col_names = c("contig", "wStart",
                                 "wEnd", "value")) %>%
  mutate(contig = paste0("chr", sub("SUPER_", "", contig)))

cov <- read_tsv(covFile,
                col_names = c("contig", "wStart",
                              "wEnd", "value")) %>%
  mutate(contig = paste0("chr", sub("SUPER_", "", contig)))


# Three kinds of panels: 
# 1) multifactor feature in big windows (Nigons)
# 2) feature in sparse ranges (telomeric mappings)
# 3) single feature in small windows (GC, repeats, telomeric repeats)

# The tricky part is to include the same panels for telomeric mappings
# even if they are have no telomere. I will add 1 count at the very end
# of each contig to keep axis and levels.


# Filter data
fbusco <- filter(busco, !Status %in% c("Missing")) %>%
  left_join(nigonDict, by = c("Busco_id" = "Orthogroup")) %>%
  mutate(nigon = ifelse(is.na(nigon), "-", nigon),
         stPos = start) %>%
  filter(nigon != "-")


consUsco <- group_by(fbusco, Sequence) %>%
  mutate(nGenes = n(),
         mxGpos = max(stPos)) %>%
  ungroup() %>%
  filter(nGenes > minimumGenesPerSequence,
         mxGpos > windwSize * 2) %>%
  mutate(Sequence = paste0("chr", sub("SUPER_", "", Sequence)))

teloRepsForPlot <- filter(telomWind, contig %in% consUsco$Sequence)
allRepsForPlot <- filter(repsWind, contig %in% consUsco$Sequence)
gcForPlot <- filter(gcWind, contig %in% consUsco$Sequence,
                    value < quantile(value, 0.98),
                    value > quantile(value, 0.02))

covForPlot <- filter(cov, contig %in% consUsco$Sequence,
                     value < quantile(value, 0.98))

seqSizes <- group_by(gcForPlot, contig) %>%
  slice_max(wEnd) %>%
  ungroup() %>%
  select(Sequence = contig, size = wEnd)


if(nrow(teloMappings) > 0){
  # get telomere mapping coordinates into blocks
  teloBlocks <- block_mappings(teloMappings)
  
  longSeqTeloMappings <- filter(teloMappings,
                                tp == "P",
                                #map_quality >= 50,
                                map_length > query_length * 0.8,
                                target_length > windwSize * 2)
  
  mappedTelo <- mutate(longSeqTeloMappings,
                       frac_target_start = (target_start / target_length)) %>%
    group_by(target_name) %>%
    mutate(tReads = n()) %>%
    ungroup() %>%
    filter(tReads > minFracAlignedTeloReads * max(tReads)) %>%
    select(-tReads) %>%
    select(target_name, target_start, target_end, strand, frac_target_start)
  
  # Add missing levels
  fLev_mappedTelo <- filter(teloRepsForPlot, 
                            !contig %in% mappedTelo$target_name) %>%
    group_by(contig) %>%
    arrange(desc(wEnd)) %>% slice(1) %>%
    ungroup() %>% mutate(strand = "+", frac_target_start = 1) %>%
    select(target_name = contig, target_start = wStart,
           target_end = wEnd, strand, frac_target_start) %>%
    bind_rows(mappedTelo)
  
  # Ensure the longest sequence has at least one count
  if(max(fLev_mappedTelo$target_start) != max(seqSizes$size)){
    fLev_mappedTelo <- slice_max(seqSizes, size) %>%
      mutate(strand = "+", frac_target_start = 1,
             target_end = size-1) %>%
      select(target_name = Sequence, target_start = size,
             target_end, frac_target_start) %>%
      bind_rows(fLev_mappedTelo)
  }
  
  # Ensure there is at least one count at
  if(min(fLev_mappedTelo$target_start) != 0){
    fLev_mappedTelo <- slice(seqSizes, 1) %>%
      mutate(frac_target_start = 0,
             size = 0, target_end = size+1) %>%
      select(target_name = Sequence, target_start = size,
             target_end, frac_target_start) %>%
      bind_rows(fLev_mappedTelo)
  }
  
} else {
  teloBlocks <- tibble(target_name = character(), strand = character())
  mappedTelo <- tibble(target_name = character(), strand = character())
  
  fLev_mappedTelo <- group_by(teloRepsForPlot, contig) %>%
    arrange(desc(wEnd)) %>% slice(c(1, n())) %>%
    ungroup() %>% mutate(strand = "+") %>%
    select(target_name = contig, target_start = wStart,
           target_end = wEnd, strand)
}

fLev_mappedTelo <- filter(fLev_mappedTelo, !is.na(strand))



### Plot ####
if(nrow(consUsco) > 0){
  plNigon <- group_by(consUsco, Sequence) %>%
    mutate(ints = as.numeric(as.character(cut(stPos,
                                              breaks = seq(0, max(stPos), windwSize),
                                              labels = seq(windwSize, max(stPos), windwSize)))),
           ints = ifelse(is.na(ints), max(ints, na.rm = T) + windwSize, ints)) %>%
    count(ints, nigon) %>%
    ungroup() %>%
    mutate(Sequence = factor(Sequence,
                             levels = mixedsort(unique(Sequence)))) %>%
    ggplot(aes(fill=nigon, y=n, x=ints-(windwSize/2))) + 
    facet_grid(. ~ Sequence) +
    geom_bar(position="stack", stat="identity") +
    theme_bw() +
    scale_y_continuous("Nigon loci", 
                       breaks = scales::pretty_breaks(4),
                       position = "left") +
    scale_x_continuous(labels = label_number(scale_cut = cut_short_scale()),
                       limits = c(0, max(gcForPlot$wEnd))) +
    scale_fill_manual(values = cols) +
    guides(fill = guide_legend(nrow = 1,
                               title = "Nigon")) +
    ggtitle("") +
    theme(axis.title.x=element_blank())
  
  plTeloCov <- mutate(fLev_mappedTelo, target_name = factor(target_name,
                                                            levels = mixedsort(unique(target_name)))) %>%
    ggplot(aes(x = target_start, color = strand)) +
    facet_grid(. ~ target_name) +
    geom_histogram(bins = 100, alpha=1, position="identity",
                   fill = "white") +
    theme_bw() +
    scale_y_continuous("Telomeric reads", position = "left",
                       expand = c(0, 0.1)) +
    scale_x_continuous(labels = label_number(scale_cut = cut_short_scale()),
                       limits = c(0, max(gcForPlot$wEnd))) +
    ggtitle("") +
    guides(color = guide_legend(nrow = 1,
                               title = "Telomere strand")) +
    theme(axis.title.x=element_blank())
  
  plTelo <- mutate(teloRepsForPlot, contig = factor(contig,
                                                    levels = mixedsort(unique(contig)))) %>%
    ggplot(aes(x=wStart, y=value)) + 
    facet_grid(. ~ contig) +
    geom_line() +
    theme_bw() +
    scale_y_continuous("Telomeric repeat", position = "left") +
    scale_x_continuous(labels = label_number(scale_cut = cut_short_scale()),
                       limits = c(0, max(gcForPlot$wEnd))) +
    ggtitle("") +
    theme(axis.title.x=element_blank())
  
  plGC <- mutate(gcForPlot, contig = factor(contig,
                                            levels = mixedsort(unique(contig)))) %>%
    ggplot(aes(x=wStart, y=value)) + 
    facet_grid(. ~ contig) +
    geom_point(alpha=0.1, color = "gray") +
    theme_bw() +
    scale_y_continuous("GC", labels = scales::percent, position = "left") +
    scale_x_continuous(labels = label_number(scale_cut = cut_short_scale()),
                       limits = c(0, max(gcForPlot$wEnd))) +
    ggtitle("") +
    theme(axis.title.x=element_blank())
  
  plCov <- mutate(covForPlot, contig = factor(contig,
                                            levels = mixedsort(unique(contig)))) %>%
    ggplot(aes(x=wStart, y=value)) + 
    facet_grid(. ~ contig) +
    geom_point(alpha=0.1, color = "gray") +
    theme_bw() +
    scale_y_continuous("Coverage", position = "left") +
    scale_x_continuous(labels = label_number(scale_cut = cut_short_scale()),
                       limits = c(0, max(gcForPlot$wEnd))) +
    ggtitle("") +
    theme(axis.title.x=element_blank())
  
  plReps <- mutate(allRepsForPlot, contig = factor(contig,
                                                   levels = mixedsort(unique(contig)))) %>%
    ggplot(aes(x=wStart, y=value)) + 
    facet_grid(. ~ contig) +
    geom_point(alpha=0.1, color = "gray") +
    theme_bw() +
    scale_y_continuous("Repeat density",
                       labels = scales::percent,
                       position = "left",
                       limits = c(0, max(allRepsForPlot$value))) +
    scale_x_continuous(labels = label_number(scale_cut = cut_short_scale()),
                       limits = c(0, max(gcForPlot$wEnd))) +
    ggtitle("") +
    theme(axis.title.x=element_blank())
  
  pExpTelo <- ggarrange(plNigon, plReps, plGC, 
                        plTeloCov, plCov, 
                        nrow = 5, align = c("v"),
                        common.legend = TRUE, legend="bottom",
                        labels = LETTERS[1:5])
    
  
} else {
  plNigon <- ggplot() + theme_void()
  plTeloCov <- ggplot() + theme_void()
  plTelo <- ggplot() + theme_void()
  pExpTelo <- ggplot() + theme_void()
}



# QC metrics

# identify erroneous Nigon fusions. Assuming O. tipulae karyotype
if(any(grepl("\\+", mappedTelo$strand)) & any(grepl("-", mappedTelo$strand))) {
  teloFracByStrand <- mutate(mappedTelo, cstrand = ifelse(strand == "+", "pos", "neg")) %>%
    group_by(target_name, cstrand) %>%
    summarise(avgFrac = mean(frac_target_start),
              .groups = "drop") %>%
    pivot_wider(names_from = cstrand,
                id_cols = target_name,
                values_from = avgFrac)
  seqsWithInternalTelomere <- filter(teloFracByStrand, pos > .02, neg < .98) %>%
    pull(target_name)
} else {
  seqsWithInternalTelomere <- character()
}


if(nrow(consUsco > 0)){
  nigonFused <- filter(consUsco, nigon != "-") %>%
    count(Sequence, nigon) %>%
    group_by(Sequence) %>%
    mutate(nSeqNigons = sum(n),
           fracSeqNigon = n/nSeqNigons,
           maxFracSeqNigon = max(fracSeqNigon)) %>%
    filter(nSeqNigons > minimumGenesPerSequence,
           fracSeqNigon > 0.2 * maxFracSeqNigon) %>%
    mutate(ndiffNigons = n()) %>%
    filter(ndiffNigons > 1) %>%
    arrange(nigon) %>%
    summarise(diffNigons = paste(nigon, collapse = ","),
              nSeqNigons = unique(nSeqNigons),
              .groups = "drop") %>%
    filter(diffNigons != "E,X")
} else {
  nigonFused <- tibble(Sequence = "",
                       diffNigons = "",
                       nSeqNigons = 0)
}

nigonFusedAndInternalTelomere <- filter(nigonFused, Sequence %in%
                                          seqsWithInternalTelomere)


# I add an artificicial duplicate in case in doesn't exist
# to be able to create a "Duplicated" column in the busco_score step ahead.
if(!"Duplicated" %in% busco$Status){
  busco <- tibble(Busco_id = "fakeBusco1", Status = "Duplicated",
                  Sequence = "NonExistent", start = 0, end = 1008100,
                  strand = "+", Score = 100, Length = 1000, OrthoDB_url = "url",
                  Description = "none") %>%
    bind_rows(busco)
  print("not Duplicated")
}

if(!"Fragmented" %in% busco$Status){
  busco <- tibble(Busco_id = "fakeBusco2", Status = "Fragmented",
                  Sequence = "NonExistent", start = 0, end = 1008100,
                  strand = "+", Score = 100, Length = 1000, OrthoDB_url = "url",
                  Description = "none") %>%
    bind_rows(busco)
  print("not Fragmented")
}

if(!"Missing" %in% busco$Status){
  busco <- tibble(Busco_id = "fakeBusco3", Status = "Missing",
                  Sequence = "NonExistent", start = 0, end = 1008100,
                  strand = "+", Score = 100, Length = 1000, OrthoDB_url = "url",
                  Description = "none") %>%
    bind_rows(busco)
  print("not Missing")
}

# calculate busco string because BUSCO sometimes fails in rounding
busco_score <- filter(busco, !duplicated(Busco_id)) %>%
  mutate(Status = sub("Complete", "Single", Status)) %>%
  count(Status) %>%
  mutate(Total = sum(n)) %>%
  pivot_wider(names_from = Status, values_from = n) %>%
  mutate(Complete = Single + Duplicated)  %>%
  pivot_longer(cols = c(Duplicated, Fragmented,
                        Missing, Single, Complete)) %>%
  mutate(frac = round((value/Total) * 100, 1)) %>%
  select(-value) %>%
  pivot_wider(names_from = name, values_from = frac)


busco_string <- paste("C:", busco_score$Complete, "%",
                      "[S:", busco_score$Single, "%,",
                      "D:", busco_score$Duplicated, "%],",
                      "F:", busco_score$Fragmented, "%,",
                      "M:", busco_score$Missing, "%,",
                      "n:", busco_score$Total,
                      sep = "")



if(any(grepl("\\+", mappedTelo$strand)) & any(grepl("-", mappedTelo$strand))) {
  teloCompleteSeqs <- filter(teloFracByStrand, pos < .05, neg > .95) %>%
    pull(target_name)
} else {
  teloCompleteSeqs <- character()
}




# This assumes low duplication rate
nigonCompleteSeqs <- count(fbusco, nigon, Sequence) %>%
  group_by(nigon) %>%
  mutate(fracTot = n/sum(n)) %>%
  ungroup() %>%
  filter(fracTot > minNigonFrac,
         nigon != "-") %>%
  select(Sequence, nigon)

t2t <- nigonCompleteSeqs$Sequence[
  nigonCompleteSeqs$Sequence %in% teloCompleteSeqs]


chromQC <- tibble(assemblyName = assemName,
                  nT2t = length(t2t),
                  nCompleteNigonSeqs = nrow(nigonCompleteSeqs),
                  nInternalTelomere = length(seqsWithInternalTelomere),
                  nNigonFusedAndInternalTelomere = nrow(nigonFusedAndInternalTelomere),
                  nNigonFused = nrow(nigonFused),
                  completeUscos = busco_score$Complete,
                  singleUscos = busco_score$Single,
                  duplicatedUscos = busco_score$Duplicated,
                  fragmentedUscos = busco_score$Fragmented,
                  missingUscos = busco_score$Missing,
                  t2tSeqs = paste(t2t, collapse = ","),
                  completeNigonSeqs = paste(nigonCompleteSeqs$Sequence, collapse = ","),
                  completeNigons = paste(nigonCompleteSeqs$nigon, collapse = ","))

# Export plot and result
write_tsv(chromQC, paste(assemName,
                         ".chromQC.tsv",
                         sep = ""))
write_tsv(teloBlocks, paste(assemName,
                            ".teloMappedBlocks.tsv",
                            sep = ""))

write(busco_string, paste(assemName,
                          ".buscoString.txt",
                          sep = ""))

plotWidth <- min(1 +
                   length(unique(consUsco$Sequence)) * 1.5,
                 50)

ggsave(paste(assemName, ".pdf", sep = ""),
       pExpTelo,
       width = plotWidth, height = 9)

ggsave(paste(assemName, ".jpg", sep = ""),
       pExpTelo, units = "in", 
       width = 8.5*1.1,
       height = 11*0.7)
