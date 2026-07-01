#!/usr/bin/env bash
# Per-sample ViroWatch pipeline.
# Called by Nextflow once per sample; all tools are available via the
# activated virowatch conda environment.
set -euo pipefail

SAMPLE=""
FASTQ=""
REF_FA=""
REF_GFF=""
THREADS=4
MODEL="r1041_e82_400bps_sup_v5.2.0"
MINLEN=2000
MAXLEN=6000
MINQUAL=10
BLAST_DB=""
CORE_NT_DB=""
REPORT_DIR=""
KRAKEN2_DB=""
KRAKEN2_CONFIDENCE=0.0
KRAKEN2_Z_MIN=-1.0
KRAKEN2_MIN_TAXA=3
SNPEFF_DB=""
SNPEFF_CONFIG=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --sample)      SAMPLE="$2";      shift 2 ;;
        --fastq)       FASTQ="$2";       shift 2 ;;
        --ref_fa)      REF_FA="$2";      shift 2 ;;
        --ref_gff)     REF_GFF="$2";     shift 2 ;;
        --threads)     THREADS="$2";     shift 2 ;;
        --model)       MODEL="$2";       shift 2 ;;
        --minlen)      MINLEN="$2";      shift 2 ;;
        --maxlen)      MAXLEN="$2";      shift 2 ;;
        --minqual)     MINQUAL="$2";     shift 2 ;;
        --blast_db)    BLAST_DB="$2";    shift 2 ;;
        --core_nt_db)  CORE_NT_DB="$2";  shift 2 ;;
        --report_dir)  REPORT_DIR="$2";  shift 2 ;;
        --kraken2_db)         KRAKEN2_DB="$2";         shift 2 ;;
        --kraken2_confidence) KRAKEN2_CONFIDENCE="$2"; shift 2 ;;
        --kraken2_z_min)      KRAKEN2_Z_MIN="$2";      shift 2 ;;
        --kraken2_min_taxa)   KRAKEN2_MIN_TAXA="$2";   shift 2 ;;
        --snpeff_db)          SNPEFF_DB="$2";          shift 2 ;;
        --snpeff_config)      SNPEFF_CONFIG="$2";      shift 2 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

OUT="${SAMPLE}"
mkdir -p "${OUT}"

# ── 1. Read QC ────────────────────────────────────────────────────────────────
echo "[$(date)] ${SAMPLE}: dedup + filter"
seqkit rmdup -o "${OUT}/${SAMPLE}.dedup.fq.gz" "${FASTQ}"
chopper -q "${MINQUAL}" --minlength "${MINLEN}" --maxlength "${MAXLEN}" \
    -i "${OUT}/${SAMPLE}.dedup.fq.gz" > "${OUT}/${SAMPLE}.filtered.fq"
FILTERED="${OUT}/${SAMPLE}.filtered.fq"

echo "[$(date)] ${SAMPLE}: NanoStat"
NanoStat -n "${SAMPLE}" -t "${THREADS}" \
    --outdir "${OUT}/nanostat" --fastq "${FILTERED}"

# ── 1b. Kraken2 read-set QC (optional) ────────────────────────────────────────
# Runs before assembly so taxonomic QC is produced even if Flye later yields no
# assembly. Skipped unless a Kraken2 DB is provided (mirrors the optional BLAST steps).
if [[ -n "${KRAKEN2_DB}" ]]; then
    echo "[$(date)] ${SAMPLE}: Kraken2 classification"
    mkdir -p "${OUT}/kraken2"
    kraken2 --db "${KRAKEN2_DB}" --threads "${THREADS}" \
        --confidence "${KRAKEN2_CONFIDENCE}" \
        --report "${OUT}/kraken2/${SAMPLE}.kraken2.report.txt" \
        --output "${OUT}/kraken2/${SAMPLE}.kraken2.output.txt" \
        "${FILTERED}"
    if [[ -n "${REPORT_DIR}" ]]; then
        python "${REPORT_DIR}/meta_kg_export.py" \
            --sample         "${SAMPLE}" \
            --kraken2-report "${OUT}/kraken2/${SAMPLE}.kraken2.report.txt" \
            --reads          "${FASTQ}" \
            --outdir         "${OUT}" \
            --z-min          "${KRAKEN2_Z_MIN}" \
            --min-taxa       "${KRAKEN2_MIN_TAXA}"
    fi
fi

# ── 2. Map QC ─────────────────────────────────────────────────────────────────
echo "[$(date)] ${SAMPLE}: minimap2 + qualimap"
minimap2 -t "${THREADS}" -ax map-ont "${REF_FA}" "${FILTERED}" \
    | samtools sort -@ "${THREADS}" -o "${OUT}/aln.bam"
samtools index "${OUT}/aln.bam"
qualimap bamqc -bam "${OUT}/aln.bam" -gff "${REF_GFF}" \
    -nt "${THREADS}" -outdir "${OUT}/qualimap" --java-mem-size=4G

# ── 3. Assembly ───────────────────────────────────────────────────────────────
echo "[$(date)] ${SAMPLE}: Flye assembly"
flye --nano-hq "${FILTERED}" -o "${OUT}/flye" -t "${THREADS}" --meta

if [[ ! -f "${OUT}/flye/assembly.fasta" ]]; then
    echo "[WARN] ${SAMPLE}: Flye produced no assembly — skipping downstream steps"
    exit 0
fi

# ── 4. Racon polish (3×) ──────────────────────────────────────────────────────
echo "[$(date)] ${SAMPLE}: Racon polish"
CONTIGS="${OUT}/flye/assembly.fasta"
for i in 1 2 3; do
    minimap2 -x map-ont -o "${OUT}/racon_overlap_${i}.paf" "${CONTIGS}" "${FILTERED}"
    racon -t "${THREADS}" "${FILTERED}" "${OUT}/racon_overlap_${i}.paf" "${CONTIGS}" \
        > "${OUT}/racon_iter_${i}.fa"
    CONTIGS="${OUT}/racon_iter_${i}.fa"
done

# ── 5. Medaka consensus ───────────────────────────────────────────────────────
echo "[$(date)] ${SAMPLE}: Medaka consensus"
micromamba run -n medaka medaka_consensus -i "${FILTERED}" -d "${CONTIGS}" \
    -o "${OUT}/medaka_consensus" -t "${THREADS}" -m "${MODEL}"
CONSENSUS="${OUT}/medaka_consensus/consensus.fasta"

# ── 6. Variant calling & annotation ───────────────────────────────────────────
# Calls variants from the filtered reads against the reference (medaka_variant),
# reheaders + fills tags, then — only when a --snpeff_db is supplied — annotates
# with SnpEff and injects INFO/GENE (mirrors the reference VARIANT subworkflow:
# MEDAKA_VARIANT → BCFTOOLS_PREP → SNPEFF → ANNOTATE_GENES). medaka/bcftools/htslib
# all live in the medaka env; snpEff runs from the virowatch env.
echo "[$(date)] ${SAMPLE}: Medaka variant calling"
mkdir -p "${OUT}/variants"
micromamba run -n medaka medaka_variant \
    -i "${FILTERED}" -r "${REF_FA}" \
    -o "${OUT}/variants/medaka_var_out" -t "${THREADS}"
cp "${OUT}/variants/medaka_var_out/medaka.annotated.vcf" "${OUT}/variants/medaka.raw.vcf"
cp "${OUT}/variants/medaka_var_out/calls_to_ref.bam"     "${OUT}/variants/calls_to_ref.bam"
cp "${OUT}/variants/medaka_var_out/calls_to_ref.bam.bai" "${OUT}/variants/calls_to_ref.bam.bai"

# Reheader sample column to the sample ID and fill INFO/END + INFO/TYPE tags
printf '%s\n' "${SAMPLE}" > "${OUT}/variants/sample_names.txt"
micromamba run -n medaka bcftools reheader \
    -s "${OUT}/variants/sample_names.txt" \
    -o "${OUT}/variants/reheadered.vcf" "${OUT}/variants/medaka.raw.vcf"
micromamba run -n medaka bcftools +fill-tags \
    -Ov -o "${OUT}/variants/tagged.vcf" "${OUT}/variants/reheadered.vcf" \
    -- -t INFO/END,INFO/TYPE

if [[ -n "${SNPEFF_DB}" ]]; then
    echo "[$(date)] ${SAMPLE}: SnpEff annotation (${SNPEFF_DB})"
    SNPEFF_CFG_FLAG=""
    [[ -n "${SNPEFF_CONFIG}" ]] && SNPEFF_CFG_FLAG="-c ${SNPEFF_CONFIG}"
    snpEff ann ${SNPEFF_CFG_FLAG} "${SNPEFF_DB}" \
        -csvStats  "${OUT}/variants/snpeff.stats.csv" \
        -htmlStats "${OUT}/variants/snpEff_summary.html" \
        "${OUT}/variants/tagged.vcf" > "${OUT}/variants/snpeff.annotated.vcf"

    # Inject INFO/GENE extracted from the SnpEff ANN field (field 4 = gene name)
    micromamba run -n medaka bcftools query \
        -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\n' "${OUT}/variants/snpeff.annotated.vcf" \
    | awk -F'\t' '{
        split($5, ann, ",");
        genes = "";
        for (i in ann) {
            split(ann[i], f, "|");
            if (f[4] != "") { genes = (genes == "") ? f[4] : genes "," f[4]; }
        }
        print $1 "\t" $2 "\t" $3 "\t" $4 "\t" genes;
      }' > "${OUT}/variants/genes.tsv"
    micromamba run -n medaka bgzip -f "${OUT}/variants/genes.tsv"
    micromamba run -n medaka tabix -s1 -b2 -e2 "${OUT}/variants/genes.tsv.gz"

    printf '##INFO=<ID=GENE,Number=.,Type=String,Description="Gene names extracted from SnpEff ANN">\n' \
        > "${OUT}/variants/gene_header.txt"
    micromamba run -n medaka bcftools annotate \
        -a "${OUT}/variants/genes.tsv.gz" \
        -c CHROM,POS,REF,ALT,INFO/GENE \
        -h "${OUT}/variants/gene_header.txt" \
        "${OUT}/variants/snpeff.annotated.vcf" \
        -Oz -o "${OUT}/variants/medaka.annotated.vcf.gz"
    # Random-access index (.gzi) for the already-bgzipped output
    micromamba run -n medaka bgzip -r "${OUT}/variants/medaka.annotated.vcf.gz"
else
    echo "[$(date)] ${SAMPLE}: no --snpeff_db, skipping SnpEff annotation"
fi

# ── 7. QUAST assembly QC ──────────────────────────────────────────────────────
echo "[$(date)] ${SAMPLE}: QUAST"
quast -o "${OUT}/quast" -t "${THREADS}" --nanopore "${FILTERED}" \
    -g "${REF_GFF}" -r "${REF_FA}" "${CONSENSUS}"

# ── 8. SierraPy drug resistance ───────────────────────────────────────────────
echo "[$(date)] ${SAMPLE}: SierraPy"
sierrapy fasta "${CONSENSUS}" --no-sharding -o "${OUT}/sierrapy.json"

# ── 9. BLAST (optional) ───────────────────────────────────────────────────────
if [[ -n "${BLAST_DB}" || -n "${CORE_NT_DB}" ]]; then
    mkdir -p "${OUT}/blast"
fi

if [[ -n "${CORE_NT_DB}" ]]; then
    echo "[$(date)] ${SAMPLE}: BLAST vs core_nt"
    blastn -task megablast \
        -db "${CORE_NT_DB}" \
        -query "${CONSENSUS}" \
        -outfmt 15 \
        -out "${OUT}/blast/core_nt.blast.json" \
        -num_threads "${THREADS}" -mt_mode 1 \
        -max_target_seqs 25 -evalue 1e-20 -perc_identity 85
fi

if [[ -n "${BLAST_DB}" ]]; then
    echo "[$(date)] ${SAMPLE}: BLAST vs LosAlamos"
    blastn -task megablast \
        -db "${BLAST_DB}" \
        -query "${CONSENSUS}" \
        -outfmt 15 \
        -out "${OUT}/blast/los_alamos.blast.json" \
        -num_threads "${THREADS}" -mt_mode 1 \
        -max_target_seqs 25 -evalue 1e-20 -perc_identity 85
fi

# ── 10. MultiQC ───────────────────────────────────────────────────────────────
echo "[$(date)] ${SAMPLE}: MultiQC"
multiqc --force --outdir "${OUT}/multiqc" "${OUT}"

# ── 11. HTML report ───────────────────────────────────────────────────────────
if [[ -n "${REPORT_DIR}" ]]; then
    echo "[$(date)] ${SAMPLE}: generating report"
    VL_FLAG=""
    CD4_FLAG=""
    [[ -n "${VL_CSV:-}"  ]] && VL_FLAG="--vl_csv ${VL_CSV}"
    [[ -n "${CD4_CSV:-}" ]] && CD4_FLAG="--cd4_csv ${CD4_CSV}"
    python "${REPORT_DIR}/render.py" \
        --sample     "${SAMPLE}" \
        --sample_dir "${OUT}" \
        --report_dir "${REPORT_DIR}" \
        --output     "${OUT}/${SAMPLE}_report.html" \
        ${VL_FLAG} ${CD4_FLAG}
fi

# ── 12. KG CSV export ─────────────────────────────────────────────────────────
if [[ -n "${REPORT_DIR}" ]]; then
    echo "[$(date)] ${SAMPLE}: KG CSV export"
    python "${REPORT_DIR}/kg_export.py" \
        --sample     "${SAMPLE}" \
        --sample_dir "${OUT}" \
        --report_dir "${REPORT_DIR}" \
        --fastq      "${FASTQ}" \
        --ref_fa     "${REF_FA}"
fi

echo "[$(date)] ${SAMPLE}: done"
