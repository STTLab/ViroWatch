# ViroWatch

ViroWatch is a bioinformatics pipeline and graph database blueprint designed to facilitate real-time molecular surveillance of HIV-1 genomes. It is optimized for low-resource settings and can provide rapid and accessible insights into viral transmission dynamics and potential drug resistance.

> [!IMPORTANT]
> Use a command line implementation of the pipeline in [HIV-64148](https://github.com/STTLab/HIV-64148) assembly pipeline.

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
| *Figure 1:* The sequencing data from Oxford Nanopore (ONT) sequencer will produce either signal file in POD5 format or raw reads FASTQ file based on the settings of the machine, the raw read will subjected to quality control, the assembly process a de novo assembly. The final assembled contigs will be identify with BLASTn to retrieve complete or near complete HIV-1 genome, these genome will be assigned a subtype based on the subtype of the closest reference genome, additionally, the genomes are sent to Stanford HIV drug resistant database for identification of drug resistance mutation which are cross-verify with locally compute variant calling and variant effect prediction. The final clinical genomic report will comprised of a strain composition found within a sample, a list of drug resistant mutation, variants and their significance along with patient’s clinical data such as historical viral load or CD4 count.

## Tools

- QC Tools (selectable)
  - [FiltLong](https://github.com/rrwick/Filtlong) (v0.2.1)
  - [LongQC](https://github.com/yfukasawa/LongQC) (v1.2.1)
- Assemblers
  - [Flye\/MetaFlye](https://github.com/mikolmogorov/Flye) (v2.9.5)
- Variant caller
  - [Medaka](https://github.com/nanoporetech/medaka) (v2.2.0) - Variant calling from reads.
  - [Snippy](https://github.com/tseemann/snippy) (latest) - Variant calling from assembled contigs.

## Data Sources

- HIV-1 sequences from Los Alamos HIV databases.
- Stanford HIV Drug Resistance Database.
- Publicly available sequences from NCBI.

## Knowledge graph design

| ![Figure 2: An illustration of entities relationship pattern for managing bacterial whole genome sequencing data and all relevant information by NosoGraph.](./README/Images/figure_2.png?raw=true "Figure 2")
|:--
| *Figure 2:* The structure of the Knowledge Graph for ViroWatch

The structure of the Knowledge Graph developed in the current system is illustrated in Figure 2. It is designed to be divided into three main domains based on interconnected clinical and biological knowledge, linked through defined relationships:

1. Clinical terminology domain represents standardized clinical concepts using SNOMED Clinical Terms (Systematized Nomenclature of Medicine—Clinical Terms) as a controlled vocabulary. This domain covers entities such as disorders, clinical findings, health-related situations, and morphologic abnormalities. The use of such standards ensures consistency in patient data recording, enables systematic disease classification, and facilitates interoperability with external clinical databases.
2. Patient and clinical metadata is used to store patient-related information, including key entities such as patients, specimens, and laboratory test results. These include measurements such as HIV viral load and CD4+ cell counts. This domain captures the clinical context by addressing “who” the patient is, “what” samples were collected, and “what” tests were performed along with their results, forming a critical foundation for downstream analysis.
3. Microbiology and genomics domain represents biological data and analytical results derived from patient samples. It includes entities such as isolates, organisms, genome assemblies, and genetic variants generated from sequencing pipelines. This domain links molecular-level information with clinical data, enabling integrated analysis of relationships between pathogen genomics and clinical or epidemiological characteristics.

## Usage

This repository provides:

- A conceptual schema defining node labels, relationship types, and data domains
- Example CSV files for data import
- Cypher queries demonstrating common operations and analytical use cases
- Guidance for setting up Neo4j as a working environment

Users can adopt the schema as a starting point, extend it to fit their specific use cases, and integrate it with custom pipelines or applications as needed.

It is important to note that, **NosoGraph is not a database management system (DBMS)** and does not provide a complete software platform for data ingestion, storage, or analysis. Instead, it defines a blueprint outlining structured conceptual model that guides how clinical, microbiological, and genomic data should be organized and linked within a graph database. The implementation of the underlying infrastructure (e.g., data pipelines, deployment environment, access control, and application interfaces) is intentionally out of scope of this repository. Users are expected to adapt the schema to their own systems and integrate it with existing workflows or tools.

We recommend using Neo4j as the platform offers an intuitive desktop interface, providing ease-of-use for general users and a mature ecosystem for graph-based development.

> [info]
> **Disclaimer:** This project is not affiliated with, endorsed by, or sponsored by Neo4j, Inc. “Neo4j” and related trademarks are the property of Neo4j, Inc. All references to Neo4j within this repository are for informational and implementation purposes only.

### Quick Start (Neo4j Desktop)

#### 1. Install Neo4j Desktop

Download and install Neo4j Desktop from:

[https://neo4j.com/download/](https://neo4j.com/download/)

Follow instructions to download, install, and launch the application.

#### 2. Create a New Database

1. Choose "Local instances" on the sidebar menu
2. Click "Create instance"
3. Fill instance details according to instructions.
4. Set a database name (e.g., nosograph-db)
5. Set a password and store it securely
6. Click “Create”.
7. Connect to the instance through "Query" or "Explore" menu

#### 3. Prepare Data Import

To import data into Neo4j instance, if using CSV files, the file must be put into an import directory within an instance path. The path can be looked up in instances list in the connection screen `Path: C:\Users\<username>\.Neo4jDesktop2\Data\dbmss\dbms-<instance-id>\import`

```Cyppher
LOAD CSV WITH HEADERS FROM 'file:///<file_name>.csv' AS row
RETURN row;
```

#### 4. Explore the Graph

From Query menu after connected to an instance you may use Neo4j Browser to:

- Visualize relationships interactively
- Expand nodes (double-click)
- Run example queries from this repository
