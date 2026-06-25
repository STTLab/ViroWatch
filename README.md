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
              → BLAST vs core_nt (optional; NCBI taxonomy)
              → MultiQC (aggregated QC report)
              → HTML HIV sequence analysis report (per sample)
```

## Requirements

- [Nextflow](https://www.nextflow.io/) ≥ 23.10.0
- [Micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html) or [Mamba](https://mamba.readthedocs.io/) (recommended; plain conda is slow on medaka's dependency tree)

Most tools are installed automatically into `envs/virowatch.yaml` on first run. Medaka runs in a separate isolated environment (`envs/medaka.yaml`) because its PyTorch/CUDA dependencies conflict with the main environment.

### Hardware note — AVX2

Newer builds of Flye (≥ 2.9.6) and Racon (1.5.0) are compiled with AVX2 instructions and will crash with `Illegal instruction (SIGILL)` on CPUs that pre-date Haswell (Intel Xeon E5 v1/v2, some Broadwell Xeons). The environment file already pins `flye=2.9.5` to avoid this. Racon 1.5.0 has no non-AVX2 conda build; a workaround is to manually replace the binary after environment creation:

```bash
# Check whether your CPU supports AVX2
grep -c avx2 /proc/cpuinfo     # 0 = no AVX2; >0 = fine, no workaround needed

# Workaround: replace racon binary with an AVX2-free build (Linux x86_64)
wget https://conda.anaconda.org/bioconda/linux-64/racon-1.4.20-hd03093a_2.tar.bz2
tar xjf racon-1.4.20-hd03093a_2.tar.bz2 -C $(conda env list | awk '/virowatch/{print $NF}') bin/racon bin/rampler
```

On AVX2-capable hardware (Haswell and later) no action is required — remove the `flye=2.9.5` pin in `envs/virowatch.yaml` to get the latest build.

## Quick start

```bash
# Test run with bundled data (CRF01_AE reference + test FASTQ)
nextflow run . -profile test

# Standard run
nextflow run . --input samplesheet.csv --outdir ./results

# With optional LosAlamos BLAST subtyping
nextflow run . --input samplesheet.csv --blast_db /path/to/LosAlamos_db

# With both BLAST DBs and clinical data
nextflow run . --input samplesheet.csv \
  --blast_db /path/to/LosAlamos_db \
  --core_nt_db /path/to/core_nt/core_nt \
  --vl_csv vl.csv --cd4_csv cd4.csv

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
| `--blast_db` | `null` | Path to pre-built LosAlamos BLAST DB (disabled if null) |
| `--core_nt_db` | `null` | Path to NCBI core_nt DB; requires `BLASTDB` env var pointing to taxdb |
| `--vl_csv` | `null` | Viral load history CSV (`sample_id`, `date`, `vl` columns) |
| `--cd4_csv` | `null` | CD4 count history CSV (`sample_id`, `date`, `cd4_pct`, `cd4_count` columns) |
| `--chopper_q` | `10` | Minimum read quality score |
| `--chopper_minlen` | `2000` | Minimum read length (bp) |
| `--chopper_maxlen` | `6000` | Maximum read length (bp) |

Bundled references: `CRF01_AE` (default) and `HXB2` — both in `assets/refs/`.

### Configuration layout

| File | Committed | Purpose |
|---|---|---|
| `nextflow.config` | yes | Pipeline defaults and profile stubs |
| `conf/base.config` | yes | Process resource caps and retry strategy |
| `conf/test.config` | yes | Test profile — bundled data, relaxed filters |
| `conf/site.config.template` | yes | Template for site-specific settings |
| `conf/site.config` | **no** (gitignored) | Filled-in site config — BLAST paths, workDir, etc. |

Machine-specific settings (BLAST DB paths, `workDir`, executor) go in `conf/site.config`, which is gitignored. Copy the template and fill it in:

```bash
cp conf/site.config.template conf/site.config
# edit conf/site.config, then:
nextflow run . --input samplesheet.csv -c conf/site.config
```

Alternatively, add site params to `~/.nextflow/config` — Nextflow auto-loads it on every run with no `-c` flag needed:

```groovy
// ~/.nextflow/config
params {
    blast_db   = '/mnt/central/BLAST/LosAlamos_db'
    core_nt_db = '/mnt/central/BLAST/core_nt/core_nt'
}
workDir = '/scratch/nextflow_work'
```

## Output structure

Each sample produces `results/<sample_id>/`:

```
nanostat/                    Read QC stats
aln.bam                      Reference-mapped reads
qualimap/                    Mapping QC
flye/                        De novo assembly
racon_iter_*.fa              Polishing intermediates
medaka_consensus/            Final consensus FASTA
quast/                       Assembly QC vs reference
sierrapy.json                Stanford HIVDB drug resistance result
blast/los_alamos.blast.json  LosAlamos BLAST results (if --blast_db)
blast/core_nt.blast.json     core_nt BLAST results (if --core_nt_db)
multiqc/                     Aggregated QC report
<sample_id>_report.html      Per-sample HIV sequence analysis report
```

## LosAlamos BLAST database setup

Two assets are bundled in `assets/blast/`:

**Option A — pre-built DB (recommended):** `LosAlamos_db.tar.gz` contains the BLAST index files and taxdb (66 MB). Extract to a persistent location and point `--blast_db` at it:

```bash
tar -xzf assets/blast/LosAlamos_db.tar.gz -C /path/to/blast_dbs/
nextflow run . --input samplesheet.csv --blast_db /path/to/blast_dbs/LosAlamos_db
```

**Option B — rebuild from FASTA:** `LosAlamos_db.gz` is the raw sequence file (20 MB compressed). Use this if you need to rebuild with custom taxid mappings:

```bash
gunzip -c assets/blast/LosAlamos_db.gz > LosAlamos_db.fa
makeblastdb -in LosAlamos_db.fa -dbtype nucl -out LosAlamos_db \
    -taxid_map sequence_to_taxid.txt -parse_seqids
```

### Database contents

15,471 HIV-1 sequences from the [Los Alamos HIV Sequence Database](https://www.hiv.lanl.gov/), with taxids assigned per subtype:

| Count | TaxID | Subtype |
|------:|------:|---------|
| 10,096 | 505185 | HIV-1 M:B |
| 2,402 | 505186 | HIV-1 M:C |
| 2,124 | 1345266 | HIV-1 M:CRF01_AE |
| 232 | 1287874 | HIV-1 M:CRF02_AG |
| 201 | 505226 | HIV-1 M:D |
| 183 | 1385609 | HIV-1 M:CRF07_BC |
| 101 | 505228 | HIV-1 M:G |
| 115 | 11676 | HIV-1 (root/unclassified) |
| 14 | 1392219 | HIV-1 M:F2 |
| 3 | 505184 | HIV-1 M:A |

CRF01_AE (13.7%) is well represented, relevant for Southeast Asian surveillance. Subtype B dominates (65.3%) reflecting the historical composition of the Los Alamos database.

## Tools

| Step | Tool | Version |
|---|---|---|
| Deduplication | seqkit rmdup | 2.10.0 |
| Read filtering | chopper | 0.10.0 |
| Read QC | NanoStat | 1.6.0 |
| Reference mapping | minimap2 | latest |
| Mapping QC | qualimap | 2.3 |
| De novo assembly | Flye (--meta) | 2.9.5 † |
| Polishing | Racon | 1.5.0 ‡ |
| Consensus | Medaka | 2.1.1 (separate env) |
| Assembly QC | QUAST | 5.3.0 |
| Drug resistance | SierraPy | 0.4.3 |
| Subtyping | BLAST+ | 2.16.0 |
| Aggregated QC | MultiQC | 1.28 |
| Report | Jinja2 / Python | 3.1.4 / 3.x |

† Pinned at 2.9.5 — the last build without AVX2; safe to upgrade on AVX2-capable hardware.  
‡ bioconda 1.5.0 binary requires AVX2; see [Hardware note](#hardware-note--avx2) for the workaround on older CPUs.

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
