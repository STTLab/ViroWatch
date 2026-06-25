# ViroWatch

ViroWatch is a Nextflow pipeline for HIV-1 genome surveillance from Oxford Nanopore reads. It takes per-sample FASTQ files through quality control, de novo assembly, consensus polishing, drug resistance analysis, and optional BLAST-based subtyping.

Designed for low-resource settings — portable, single conda environment, resume-capable.

## Pipeline overview

```
FASTQ → dedup → filter → NanoStat
              → minimap2 → qualimap
              → Flye (de novo assembly)
              → Racon ×3 (polishing)
              → Medaka (consensus)
              → QUAST (assembly QC)
              → SierraPy (Stanford HIVDB drug resistance)
              → BLAST vs LosAlamos (optional subtyping)
              → MultiQC (aggregated QC report)
```

## Requirements

- [Nextflow](https://www.nextflow.io/) ≥ 23.10.0
- Environment manager (pick one)
  - [Conda](https://docs.conda.io/)
  - [Mamba](https://mamba.readthedocs.io/) (recommended for environment resolution)
  - [Micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html) (light-weight Mamba equivalent)

All tools are installed automatically into a single conda environment (`envs/virowatch.yaml`) on first run.

## Quick start

```bash
# Test run with bundled data (CRF01_AE reference + test FASTQ)
nextflow run . -profile test

# Standard run
nextflow run . --input samplesheet.csv --outdir ./results

# With optional LosAlamos BLAST subtyping
nextflow run . --input samplesheet.csv --blast_db /path/to/LosAlamos_db

# Resume a failed run
nextflow run . --input samplesheet.csv -resume
```

## Samplesheet format

Provide a CSV file with `sample_id` and `fastq` columns:

```csv
sample_id,fastq
sample_01,/path/to/sample_01.fq.gz
sample_02,/path/to/sample_02.fq.gz
```

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--input` | *(required)* | Path to samplesheet CSV |
| `--outdir` | `./results` | Output directory |
| `--ref_fa` | `assets/refs/CRF01_AE.fa` | Reference FASTA for mapping QC |
| `--ref_gff` | `assets/refs/CRF01_AE.gff` | GFF for qualimap/QUAST |
| `--medaka_model` | `r1041_e82_400bps_sup_v5.2.0` | Medaka model — must match basecalling model |
| `--blast_db` | `null` (disabled) | Path to pre-built BLAST DB |
| `--chopper_q` | `10` | Minimum read quality score |
| `--chopper_minlen` | `2000` | Minimum read length (bp) |
| `--chopper_maxlen` | `6000` | Maximum read length (bp) |

Bundled references: `CRF01_AE` (default) and `HXB2` — both in `assets/refs/`.

## Output structure

Each sample produces `results/<sample_id>/`:

```
nanostat/             Read QC stats
aln.bam               Reference-mapped reads
qualimap/             Mapping QC
flye/                 De novo assembly
racon_iter_*.fa       Polishing intermediates
medaka_consensus/     Final consensus FASTA
quast/                Assembly QC vs reference
sierrapy.json         Stanford HIVDB drug resistance result
blast/                BLAST results (only if --blast_db provided)
multiqc/              Aggregated QC report
```

## LosAlamos BLAST database setup

The bundled `assets/blast/LosAlamos_db.gz` is a compressed FASTA. Build the BLAST database before use:

```bash
gunzip -c assets/blast/LosAlamos_db.gz > LosAlamos_db.fa
makeblastdb -in LosAlamos_db.fa -dbtype nucl -out LosAlamos_db
```

Then pass `--blast_db /path/to/LosAlamos_db` when running the pipeline.

## Tools

| Step | Tool | Version |
|---|---|---|
| Deduplication | seqkit rmdup | 2.10.0 |
| Read filtering | chopper | 0.10.0 |
| Read QC | NanoStat | 1.6.0 |
| Reference mapping | minimap2 | 2.29 |
| Mapping QC | qualimap | 2.3 |
| De novo assembly | Flye (--meta) | 2.9.6 |
| Polishing | Racon | 1.5.0 |
| Consensus | Medaka | 2.1.1 |
| Assembly QC | QUAST | 5.3.0 |
| Drug resistance | SierraPy | 0.4.3 |
| Subtyping | BLAST+ | 2.16.0 |
| Aggregated QC | MultiQC | 1.28 |

## Data sources

- Drug resistance: [Stanford HIV Drug Resistance Database](https://hivdb.stanford.edu/) (via SierraPy)
- Subtyping: [Los Alamos HIV Sequence Database](https://www.hiv.lanl.gov/)

---

## Knowledge graph (Neo4j)

> **Note:** The graph database component is under active development. Documentation below describes the planned schema.

| ![Figure 2: An illustration of entities relationship pattern for managing bacterial whole genome sequencing data and all relevant information by NosoGraph.](./README/Images/figure_2.png?raw=true "Figure 2")
|:--
| *Figure 2:* The structure of the Knowledge Graph for ViroWatch

The Knowledge Graph is structured across three interconnected domains:

1. **Clinical terminology** — standardized concepts using SNOMED CT (disorders, clinical findings, morphologic abnormalities).
2. **Patient and clinical metadata** — patient records, specimens, and lab results (viral load, CD4+ counts).
3. **Microbiology and genomics** — isolates, assemblies, and genetic variants linked to clinical data.

### Quick Start (Neo4j Desktop)

#### 1. Install Neo4j Desktop

Download from [https://neo4j.com/download/](https://neo4j.com/download/) and follow installation instructions.

#### 2. Create a New Database

1. Choose "Local instances" on the sidebar menu
2. Click "Create instance" and fill in instance details
3. Set a database name (e.g., `virowatch-db`) and a password
4. Click "Create" and connect via the "Query" or "Explore" menu

#### 3. Import Data

Place CSV files in the Neo4j import directory (`Path: C:\Users\<username>\.Neo4jDesktop2\Data\dbmss\dbms-<instance-id>\import`) then load with Cypher:

```cypher
LOAD CSV WITH HEADERS FROM 'file:///<file_name>.csv' AS row
RETURN row;
```

#### 4. Explore the Graph

Use Neo4j Browser to visualize relationships, expand nodes (double-click), and run analytical queries.

> **Disclaimer:** This project is not affiliated with, endorsed by, or sponsored by Neo4j, Inc. "Neo4j" and related trademarks are the property of Neo4j, Inc.
