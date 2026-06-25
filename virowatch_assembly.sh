#!/usr/bin/env bash
set -euo pipefail

# Check input argument
if [ $# -ne 1 ]; then
    echo "Usage: $0 sample_id_list.txt"
    exit 1
fi

SAMPLE_ID_FILE="$1"

# Check if file exists
if [ ! -f "$SAMPLE_ID_FILE" ]; then
    echo "Error: File '$SAMPLE_ID_FILE' not found!"
    exit 1
fi

function NanoStatQC {
    local prefix=$1
    local input_fastq=$2
    local outdir=$3
    echo ""
    echo "[ $(date) ] - Running NanoStat"
    echo "Input : ${input_fastq}"
    echo "Output: ${outdir}"
    echo ""
    micromamba run -n qc_tools NanoStat \
        -n "${prefix}" -t 8 \
        --outdir "${outdir}" \
        --fastq "${input_fastq}"
}

function flye_assembly {
    local prefix=$1
    local input_fastq=$2
    local outdir=$3

    echo ""
    echo "[ $(date) ] - Running Flye"
    echo "Input : ${input_fastq}"
    echo "Output: ${outdir}"
    echo ""

    micromamba run -n assemblers flye --nano-hq "${input_fastq}" -o "${outdir}" -t 16 --meta
}


RACON_ITER=3
function racon_polish {
    local contigs_fasta=$1
    local input_fastq=$2
    local outdir=$3

    mkdir -p ${outdir}
    echo ""
    echo "[ $(date) ] - Running Racon"
    echo "Input : ${input_fastq}"
    echo "Output: ${outdir}"
    echo ""
    for ((i=0; i<=RACON_ITER; i++)); do
        if [[ "$i" -gt 0 ]]; then
            overlap_file="${outdir}/overlap_iter_${i}.paf"
            micromamba run -n assemblers minimap2 \
                -x map-ont \
                -o ${overlap_file} \
                ${contigs_fasta} ${input_fastq}
            exit_code=0 #$?
            if [[ "$exit_code" -ne 0 ]]; then
                break
            fi
            polished_sequence="${outdir}/racon_iter_${i}.fa"
            micromamba run -n assemblers racon -t 16 \
                ${input_fastq} \
                ${overlap_file} \
                ${contigs_fasta} > ${polished_sequence}
            exit_code=0 #$?
            if [[ "$exit_code" -ne 0 ]]; then
                echo "ERROR polishing with racon."
                break
            else
                contigs_fasta="${polished_sequence}"
            fi
        fi
    done
    if [[ -f "${outdir}/racon_iter_${RACON_ITER}.fa" ]]; then
        ln -fs "$(pwd)/${outdir}/racon_iter_${RACON_ITER}.fa" "$(pwd)/${outdir}/final_polished.fa"
    fi
}

# Load SRA toolkit environment if needed (optional)
# export PATH=$PATH:/path/to/sratoolkit/bin

# Process each sample_id
while IFS= read -r sample_id || [[ -n "$sample_id" ]]; do
    if [[ -z "$sample_id" ]]; then
        continue  # skip empty lines
    fi
    if [[ -f "./${sample_id}.fq.gz" && ! -d "./${sample_id}" ]]; then
        mkdir "${sample_id}"
        mv "./${sample_id}.fq.gz" "./${sample_id}/${sample_id}.fq.gz"
    fi
    if [[ -f "./${sample_id}/${sample_id}.fq.gz" ]]; then
        micromamba run -n qc_tools seqkit rmdup -o "./${sample_id}/${sample_id}.dedup.fq.gz" "./${sample_id}/${sample_id}.fq.gz" 
        filtered_fq="./${sample_id}/${sample_id}.filtered.fq"
        micromamba run -n qc_tools chopper -q 10 --minlength 2000 --maxlength 6000 -i "./${sample_id}/${sample_id}.dedup.fq.gz" > $filtered_fq
        NanoStatQC ${sample_id} $filtered_fq "./${sample_id}/nanostat"
        micromamba run -n assemblers minimap2 -t 16 -ax map-ont /mnt/sdb/REPOSITORY/REF_SEQUENCES/HIV1/CRF01_AE/AF164485.1.fa $filtered_fq > "./${sample_id}/aln.sam"
        samtools sort "./${sample_id}/aln.sam" -o "./${sample_id}/aln.bam"
        rm "./${sample_id}/aln.sam"
        samtools index "./${sample_id}/aln.bam"
        samtools stats "./${sample_id}/aln.bam" > "./${sample_id}/aln.stats"
        micromamba run -n qc_tools qualimap bamqc -bam "./${sample_id}/aln.bam" -gff /mnt/sdb/REPOSITORY/REF_SEQUENCES/HIV1/CRF01_AE/AF164485.1.gff -nt 16 -outdir "./${sample_id}/qualimap"
        micromamba run -n qc_tools qualimap comp-counts -bam "./${sample_id}/aln.bam" -gtf /mnt/sdb/REPOSITORY/REF_SEQUENCES/HIV1/CRF01_AE/AF164485.1.gff -out "./${sample_id}/qualimap/comp-counts.counts"
        flye_assembly ${sample_id} $filtered_fq "./${sample_id}/flye"
        contigs_fa="./${sample_id}/flye/assembly.fasta"
        if [[ ! -f "./${sample_id}/flye/assembly.fasta" ]]; then
            continue
        fi
        racon_polish $contigs_fa \
            $filtered_fq \
            "./${sample_id}/racon"
        racon_polished_fa="./${sample_id}/racon/final_polished.fa"
        micromamba run -n medaka medaka_consensus -i $filtered_fq -d $racon_polished_fa -o "./${sample_id}/medaka_consensus" -t 8 -m r1041_e82_400bps_sup_v5.2.0
        medaka_polished_fa="./${sample_id}/medaka_consensus/consensus.fasta"
        micromamba run -n qc_tools quast \
            -o "./${sample_id}/quast" -t 8 \
            --glimmer --nanopore $filtered_fq \
            -g /mnt/sdb/REPOSITORY/REF_SEQUENCES/HIV1/CRF01_AE/AF164485.1.gff \
            -r /mnt/sdb/REPOSITORY/REF_SEQUENCES/HIV1/CRF01_AE/AF164485.1.fa \
            $medaka_polished_fa
        micromamba run -n sierrapy sierrapy fasta $medaka_polished_fa -o "./${sample_id}/sierrapy_result.json"

        export BLASTDB="/mnt/sdb/REPOSITORY/BLAST/taxdb"
        BLAST_CORE_NT_DB="/mnt/sdb/REPOSITORY/BLAST/core_nt/core_nt"
        BLAST_LOS_ALAMOS_DB="/mnt/sdb/REPOSITORY/BLAST/HIV64148_LosAlamos/LosAlamos_db"
        mkdir "./${sample_id}/blast"
        micromamba run -n blast blastn -task megablast -db $BLAST_CORE_NT_DB -query $medaka_polished_fa \
            -outfmt "7 qacc sacc staxid sscinames sblastnames qstart qend sstart send qcovs pident evalue bitscore" \
            -num_threads 8 -mt_mode 1 \
            -max_target_seqs 25 \
            -evalue 1e-20 \
            -perc_identity 85 \
            -out "./${sample_id}/blast/medaka_consensus.fa.core_nt.blast.tsv"
        micromamba run -n blast blastn -task megablast -db $BLAST_LOS_ALAMOS_DB -query $medaka_polished_fa \
            -outfmt 11 \
            -num_threads 8 -mt_mode 1 \
            -max_target_seqs 25 \
            -evalue 1e-20 \
            -perc_identity 85 \
            -out "./${sample_id}/blast/medaka_consensus.fa.los_alamos.blast.asn"
        micromamba run -n blast blast_formatter \
            -archive "./${sample_id}/blast/medaka_consensus.fa.los_alamos.blast.asn" \
            -outfmt "7 qacc sacc staxid sscinames stitle qstart qend sstart send qcovs pident evalue bitscore" \
            -out "./${sample_id}/blast/medaka_consensus.fa.los_alamos.blast.tsv"


        micromamba run -n qc_tools multiqc --force --outdir "./${sample_id}/multiqc" "./${sample_id}"
    else
        echo $sample_id
    fi
done < "$SAMPLE_ID_FILE"
