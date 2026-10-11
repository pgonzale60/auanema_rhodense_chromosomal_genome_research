# Nematode rDNA Phylogenetic Analysis (15 Taxa)

This directory contains the dataset preparation, sequence alignments, and phylogenetic tree inference supporting Panel A of **Figure 3** in the manuscript.

## Overview

The molecular cladogram features 15 nematode taxa across Rhabditina, Ascaridomorpha, Spiruromorpha, and the outgroup *Trichinella spiralis* to contextualize Programmed DNA Elimination (PDE) phenotypes:
- **Precise PDE**: *Auanema rhodense*, *Oscheius tipulae*, *Caenorhabditis auriculariae*, *Caenorhabditis monodelphis*, *Caenorhabditis parvicauda*, *Mesorhabditis belari*
- **No PDE**: *Caenorhabditis briggsae*, *Caenorhabditis elegans*, *Pristionchus pacificus*, *Brugia malayi*
- **Imprecise PDE**: *Toxocara canis*, *Parascaris univalens*, *Ascaris lumbricoides*, *Ascaris suum*
- **Outgroup**: *Trichinella spiralis*

## Workflow

1. **Sequence Retrieval**:
   - `build_dataset.py`: Fetches GenBank 18S and 28S rRNA sequences, incorporating *C. auriculariae* rDNA directly from long-read genome assembly annotations (`c_auriculariae_18S.fa`, `c_auriculariae_28S.fa`).
   - Produces `data/all_18S.fasta` and `data/all_28S.fasta`.

2. **Multiple Sequence Alignment**:
   ```bash
   mafft --auto data/all_18S.fasta > data/all_18S.aligned.fa
   mafft --auto data/all_28S.fasta > data/all_28S.aligned.fa
   ```

3. **Concatenation & Tree Inference**:
   - Alignments concatenated into `data/concatenated.fasta`.
   - Maximum Likelihood tree inferred using IQ-TREE:
     ```bash
     iqtree3 -s data/concatenated.fasta -o Trichinella_spiralis -pre data/iqtree
     ```

4. **Cladogram Ordering**:
   - The tree is ordered to place the focal species *Auanema rhodense* at the top and outgroup *Trichinella spiralis* at the bottom: `data/ordered_cladogram.nwk`.

5. **Figure Generation**:
   - `data/ordered_cladogram.nwk` is directly ingested by `scripts/manuscript/fig_3_motif_analysis.R` to render Figure 3 Panel A.
