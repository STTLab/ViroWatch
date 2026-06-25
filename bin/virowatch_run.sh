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

# ── 6. QUAST assembly QC ──────────────────────────────────────────────────────
echo "[$(date)] ${SAMPLE}: QUAST"
quast -o "${OUT}/quast" -t "${THREADS}" --nanopore "${FILTERED}" \
    -g "${REF_GFF}" -r "${REF_FA}" "${CONSENSUS}"

# ── 7. SierraPy drug resistance ───────────────────────────────────────────────
echo "[$(date)] ${SAMPLE}: SierraPy"
sierrapy fasta "${CONSENSUS}" --no-sharding -o "${OUT}/sierrapy.json"

# ── 8. BLAST (optional) ───────────────────────────────────────────────────────
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

# ── 9. MultiQC ────────────────────────────────────────────────────────────────
echo "[$(date)] ${SAMPLE}: MultiQC"
multiqc --force --outdir "${OUT}/multiqc" "${OUT}"

# ── 10. HTML report ───────────────────────────────────────────────────────────
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

echo "[$(date)] ${SAMPLE}: done"
