# 🧬 Origin, Duplication, and Diversification of Opsins in Vertebrates

## 📖 Overview

**This project was my final year dissertation project as part of my MSc Bioinformatics course.**

The project focuses on the evolutionary history and functional diversification of opsin genes in vertebrates. 
Opsins, as key light-sensitive proteins in the G protein-coupled receptor (GPCR) family, play essential roles in the visual systems of many organisms, enabling a variety of light detection functions. 
By studying opsins, we gain insights into the molecular evolution that allows species to adapt to different light environments and visual requirements.

The study examines opsin gene duplications and losses, and their impact on gene expression patterns across vertebrate lineages.
We used bioinformatics tools, including phylogenomic analyses and single-cell RNA sequencing to reconstruct the evolutionary history of opsins.
The analysis includes data from 364 vertebrate species and over 500 known opsin sequences, covering a broad taxonomic range.

Our methodology integrates multiple approaches:

    Data Collection: A total of 1115 whole vertebrate species proteomes were downloaded from the two databases, from which 364 species were manually selected based on their BUSCO
    (Benchmarking Universal Single-copy Orthologs) completeness score, ensuring taxonomical diversity and a total of 525 available known opsin sequences were downloaded 
    from the UniProtKB database (367 query sequences), and the GPCRdb database (158 query sequences).
    We searched the literature for Single-cell transcriptomic datasets of the eye/retina and collected 23 datasets from previous studies. We also obtained the 
    corresponding proteomes from the genomes the datasets were mapped to.
    
    Phylogenetic Analysis: To explore the relationships between opsins, we constructed phylogenetic trees, helping to identify conserved and divergent evolutionary paths among these genes.
    
    Single-Cell Analysis: Single-cell RNA sequencing data was used to analyse opsin gene expression in specific tissues, such as the retina, across different species. 
    
    This approach allowed us to speculate about cell-type-specific expression patterns and examine how opsins function in various environmental contexts.

🌟 Key Findings

    Significant gene duplications in specific opsin subfamilies.

    Identification of conserved opsins expressed across a wide range of vertebrates.

    Insights into functional divergence and cell-type-specific roles of opsins in the vertebrate retina.

These results enhance our understanding of visual system evolution and highlight the genomic and transcriptomic complexity underlying light detection mechanisms.

---

## Phylogenetic Analysis

- Data mining and pre-processing
- Homologous identification
- Generation of Phylogenetic trees

Conducted multiple sequence alignment, trimming, and constructed phylogenetic trees to infer evolutionary relationships.

---

## Single-cell RNA seq. Analysis

- Processed and analyzed scRNA-seq data using the Seurat R package.
- Investigated cell-type-specific expression of opsins in chick and zebrafish retina.
  
---

## 🛠️ Tools & Environment

    Linux Mint 21.1 “Vera” (based on Ubuntu 22.04)

    Bash (Bourne Again Shell)

    ALICE High-Performance Computing (University of Leicester)

    Tools: Seurat, MAFFT, TrimAl, IQ-TREE, ggtree, rotl, diamond

## 🧭 Project Structure

01_data_collection/

Scripts for downloading proteomes and known opsin sequences.
02_homolog_search/

DIAMOND BLASTP searches, FASTA header renaming, and homolog filtering.
03_phylogenetic_analysis/

Multiple sequence alignment, trimming, and opsin phylogenetic tree generation.
04_single_cell_analysis/

Seurat-based pipelines for analyzing chick and zebrafish retinal scRNA-seq datasets.
