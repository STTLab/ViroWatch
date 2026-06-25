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

while [[ $# -gt 0 ]]; do
    case $1 in
        --sample)   SAMPLE="$2";   shift 2 ;;
        --fastq)    FASTQ="$2";    shift 2 ;;
        --ref_fa)   REF_FA="$2";   shift 2 ;;
        --ref_gff)  REF_GFF="$2";  shift 2 ;;
        --threads)  THREADS="$2";  shift 2 ;;
        --model)    MODEL="$2";    shift 2 ;;
        --minlen)   MINLEN="$2";   shift 2 ;;
        --maxlen)   MAXLEN="$2";   shift 2 ;;
        --minqual)  MINQUAL="$2";  shift 2 ;;
        --blast_db) BLAST_DB="$2"; shift 2 ;;
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
medaka_consensus -i "${FILTERED}" -d "${CONTIGS}" \
    -o "${OUT}/medaka_consensus" -t "${THREADS}" -m "${MODEL}"
CONSENSUS="${OUT}/medaka_consensus/consensus.fasta"

# ── 6. QUAST assembly QC ──────────────────────────────────────────────────────
echo "[$(date)] ${SAMPLE}: QUAST"
quast -o "${OUT}/quast" -t "${THREADS}" --nanopore "${FILTERED}" \
    -g "${REF_GFF}" -r "${REF_FA}" "${CONSENSUS}"

# ── 7. SierraPy drug resistance ───────────────────────────────────────────────
echo "[$(date)] ${SAMPLE}: SierraPy"
sierrapy fasta "${CONSENSUS}" -o "${OUT}/sierrapy.json"

# ── 8. BLAST — LosAlamos (optional) ──────────────────────────────────────────
if [[ -n "${BLAST_DB}" ]]; then
    echo "[$(date)] ${SAMPLE}: BLAST vs LosAlamos"
    mkdir -p "${OUT}/blast"
    blastn -task megablast \
        -db "${BLAST_DB}" \
        -query "${CONSENSUS}" \
        -outfmt 11 \
        -out "${OUT}/blast/result.asn" \
        -num_threads "${THREADS}" -mt_mode 1 \
        -max_target_seqs 25 -evalue 1e-20 -perc_identity 85
    blast_formatter \
        -archive "${OUT}/blast/result.asn" \
        -outfmt "7 qacc sacc staxid sscinames stitle qstart qend sstart send qcovs pident evalue bitscore" \
        -out "${OUT}/blast/los_alamos.blast.tsv"
fi

# ── 9. MultiQC ────────────────────────────────────────────────────────────────
echo "[$(date)] ${SAMPLE}: MultiQC"
multiqc --force --outdir "${OUT}/multiqc" "${OUT}"

echo "[$(date)] ${SAMPLE}: done"
