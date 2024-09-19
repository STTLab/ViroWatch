# ViroWatch

ViroWatch is a bioinformatics software pipeline designed to facilitate real-time molecular surveillance of viral genomes. This lightweight, graph-based system efficiently handles and aggregates genomic sequencing data, particularly focusing on HIV-1 genomes. It is optimized for low-resource settings and can process data from multiple sequencing platforms, providing rapid and accessible insights into viral transmission dynamics and potential drug resistance.

> [!IMPORTANT]  
> This project is currently under active development, no release have been made.

## Introduction

ViroWatch aims to enhance the management of genomic sequencing data using graph-based approaches, improving the scalability and flexibility of viral surveillance systems. It is designed to analyze HIV-1 genomes from clinical samples and public databases, integrating molecular and clinical data for comprehensive viral monitoring.

### Key objectives

- Integrating molecular surveillance with clinical data and existing virological databases.
- Developing a scalable and portable system for rapid deployment, even in low-resource settings.
- ViroWatch enables public health authorities and researchers to:
  - Track the emergence of drug-resistant HIV-1 strains.
  - Monitor potential viral transmission clusters.

### Features

- **Graph-Based Data Management:** Utilizes graph databases for flexible and scalable storage of genomic data.
- **Compatibility with Multiple Sequencing Platforms:** Supports long-read sequencing technologies (e.g., Oxford Nanopore, PacBio).
- **Real-Time Surveillance:** Enables rapid analysis and visualization of transmission clusters and genomic diversity.
- **Portable and Lightweight:** Runs on standard computer systems, making it accessible to labs in low-resource settings.
- **Data Integration:** Combines clinical data with existing virological and molecular epidemiology datasets.

## Pipeline

| ![Figure 1: An illustration of entities relationship pattern for managing bacterial whole genome sequencing data and all relevant information by NosoGraph.](./README/Images/figure_1.png?raw=true "Figure 1")
|:--
| *Figure 1:* The sequencing data from Oxford Nanopore (ONT) sequencer will produce either signal file in POD5 format or raw reads FASTQ file based on the settings of the machine, alternatively, users can input raw read FASTQ file from other long read sequencing platform such as PacBio as well, the raw read will subjected to quality control, the first process is to identify the composition of organisms within the read set using Kraken2 to determine the identity of each read, with this data, we can identify contaminants or artifact sequences and remove them along with low quality reads in the QC process prior to assembly. For the assembly process user can choose to do a de novo assembly or reference-based assembly or both, Trycycler will be used to merge contigs from both type of assemblers together, however, this process requires user to manually curate each contig themselves before putting it into Trycycler. The final assembled contigs will be identify with BLASTn to retrieve complete or near complete HIV-1 genome, these genome will be assigned a subtype based on the subtype of the closest reference genome, additionally, the genomes are sent to Stanford HIV drug resistant database for identification of drug resistance mutation which are cross-verify with locally compute variant calling and variant effect prediction. The final clinical genomic report will comprised of a strain composition found within a sample, a list of drug resistant mutation, variants and their significance along with patient’s clinical data such as historical viral load or CD4 count.

## Tools

- QC Tools (selectable)
  - [FiltLong](https://github.com/rrwick/Filtlong) (v0.2.1)
  - [LongQC](https://github.com/yfukasawa/LongQC) (v1.2.1)
- Taxonomic sequence classification
  - [Kraken2](https://github.com/DerrickWood/kraken2) (v2.1.3)
- Assemblers (selectable)
  - [Canu](https://github.com/marbl/canu) (v2.2)
  - [Flye\/MetaFlye](https://github.com/mikolmogorov/Flye) (v2.9.5)
  - [HaploDMF](https://github.com/dhcai21/HaploDMF) (latest)
- Variant caller
  - [NanoCaller](https://github.com/WGLab/NanoCaller) (v3.6.0) - Variant calling from reads.
  - [Snippy](https://github.com/tseemann/snippy) (latest) - Variant calling from assembled contigs.

## Data Sources

- HIV-1 sequences from Los Alamos HIV databases.
- Stanford HIV Drug Resistance Database.
- Publicly available sequences from NCBI.

## Knowledge graph design

| ![Figure 2: An illustration of entities relationship pattern for managing bacterial whole genome sequencing data and all relevant information by NosoGraph.](./README/Images/figure_2.png?raw=true "Figure 2")
|:--
| *Figure 2:* A knowledge graph is a structured representation of knowledge, typically consisting of entities, their attributes, and the relationships between them. This structure enables powerful data modeling, exploration, and analysis, making complex relationships within data more accessible and actionable. Knowledge graphs are derived from organizing unstructured information into a graph format, where nodes represent entities (such as people, places, or concepts) and edges represent the relationships between these entities. The knowledge graph for this study is illustrated in Figure 1. Starting from patients participate in sequencing experiment which follows certain protocol to produce sequencing data which then subjected to alignment with reference genome and/or de novo assembly and strain identification. Variants called from the sample may contain ARV resistant mutation which can be linked back to ARV medication which the patient may taking. This structure provides robust data transparency and discoverability.
