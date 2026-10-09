# Nematode rDNA Phylogenetic Analysis (13 Taxa)

This directory contains the dataset preparation, sequence alignments, and phylogenetic tree inference supporting Panel A of **Figure 3** in the manuscript.

## Overview

The molecular cladogram features 13 nematode taxa across Rhabditina, Ascaridomorpha, and the outgroup *Trichinella spiralis* to contextualize Programmed DNA Elimination (PDE) phenotypes:
- **Precise PDE**: *Auanema rhodense*, *Oscheius tipulae*, *Mesorhabditis belari*, *Caenorhabditis auriculariae*, *Caenorhabditis monodelphis*, *Caenorhabditis parvicauda*
- **No PDE**: *Caenorhabditis elegans*, *Caenorhabditis briggsae*
- **Imprecise PDE**: *Toxocara canis*, *Parascaris univalens*, *Ascaris lumbricoides*, *Ascaris suum*
- **Outgroup**: *Trichinella spiralis*

## Workflow

1. **Sequence Retrieval**:
   - `build_dataset.py`: Fetches GenBank 18S and 28S rRNA sequences for 12 species, incorporating *C. auriculariae* rDNA directly from long-read genome assembly annotations (`c_auriculariae_18S.fa`, `c_auriculariae_28S.fa`).
   - Produces `data/all_18S_13taxa.fasta` and `data/all_28S_13taxa.fasta`.

2. **Multiple Sequence Alignment**:
   ```bash
   mafft --auto data/all_18S_13taxa.fasta > data/all_18S_13taxa.aligned.fa
   mafft --auto data/all_28S_13taxa.fasta > data/all_28S_13taxa.aligned.fa
   ```

3. **Concatenation & Tree Inference**:
   - Alignments concatenated into `data/concatenated_13taxa.fasta` (6,297 bp).
   - Maximum Likelihood tree inferred using IQ-TREE with GTR+F+I+R4 model:
     ```bash
     iqtree3 -s data/concatenated_13taxa.fasta -o Trichinella_spiralis -pre data/iqtree_13taxa
     ```

4. **Cladogram Ordering**:
   - The tree is ordered to place the focal species *Auanema rhodense* at the top and outgroup *Trichinella spiralis* at the bottom: `data/ordered_cladogram_13taxa.nwk`.

5. **Figure Generation**:
   - `data/ordered_cladogram_13taxa.nwk` is directly ingested by `scripts/manuscript/fig_3_motif_analysis.R` to render Figure 3 Panel A.
