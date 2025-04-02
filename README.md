# Origin, Duplication, and diversification of Opsins in Vertebrates

## Overview

**This project was my final year dissertation project as part of my MSc Bioinformatics course.**

The project focuses on the evolutionary history and functional diversification of opsin genes in vertebrates. 
Opsins, as key light-sensitive proteins in the G protein-coupled receptor (GPCR) family, play essential roles in the visual systems of many organisms, enabling a variety of light detection functions. 
By studying opsins, we gain insights into the molecular evolution that allows species to adapt to different light environments and visual requirements.

The study examines opsin gene duplications and losses, and their impact on gene expression patterns across vertebrate lineages.
We utilised a combination of bioinformatics tools, including phylogenomic analyses and single-cell RNA sequencing to reconstruct the evolutionary history of opsins.
The analysis includes data from 364 vertebrate species and over 500 known opsin sequences, covering a wide taxonomic range.

Our methodology integrates multiple approaches:

    Data Collection: A total of 1115 whole vertebrate species proteomes were downloaded from the two databases, from which 364 species were manually selected based on their BUSCO
    (Benchmarking Universal Single-copy Orthologs) completeness score, ensuring taxonomical diversity and a total of 525 available known opsin sequences were downloaded 
    from the UniProtKB database (367 query sequences), and the GPCRdb database (158 query sequences).
    We searched the literature for Single-cell transcriptomic datasets of the eye/retina and collected 23 datasets from previous studies. We also obtained the 
    corresponding proteomes from the genomes the datasets were mapped to.
    
    Phylogenetic Analysis: To explore the relationships between opsins, we constructed phylogenetic trees, helping to identify conserved and divergent evolutionary paths among these genes.
    
    Single-Cell Analysis: Single-cell RNA sequencing data was used to analyse opsin gene expression in specific tissues, such as the retina, across different species. 
    
    This approach allowed us to speculate about cell-type-specific expression patterns and examine how opsins function in various environmental contexts.

Key findings include:

    Significant gene duplications in certain opsin families.
    Identification of conserved opsins that are expressed across a broad range of vertebrates.

This project provides a detailed exploration of opsin evolution, offering valuable insights into the molecular mechanisms that drive visual diversity in vertebrates.
Our findings contribute to a better understanding of the genetic basis for light perception and adaptation, with potential applications in vision science and evolutionary biology.

---

## Tools

The analysis was performed on a Linux laptop (intel Core i7, Linux Mint 21.1 “Vera”
operating system, based on Ubuntu Linux 22.04, with bash (Bourne again shell) as the
default shell), and we used the ALICE High-Performance Computing facility at the
University of Leicester.

---

## Phylogenetic Analysis

- Data mining and pre-processing
- Homologous identification
- Generation of Phylogenetic trees

---

## Single-cell RNA seq. Analysis

- Analysis using the Seurat R package pipeline

---

## Contents

**Data Collection**: Scripts to download the necessary data

**Phylogenetic Analysis**: Scripts to process the data and generate the phylogenetic tree

- Preparing the files
- Sorting the genomes
