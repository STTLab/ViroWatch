#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process VIROWATCH_SAMPLE {
    tag "$meta.id"
    publishDir "${params.outdir}", mode: 'copy'
    conda "${projectDir}/envs/virowatch.yaml"

    input:
    tuple val(meta), path(fastq)
    path ref_fa
    path ref_gff

    output:
    path "${meta.id}/", emit: results

    script:
    def blast_flag    = params.blast_db   ? "--blast_db ${params.blast_db}"      : ""
    def core_nt_flag  = params.core_nt_db ? "--core_nt_db ${params.core_nt_db}"  : ""
    def vl_flag       = params.vl_csv     ? "VL_CSV=${params.vl_csv}"            : ""
    def cd4_flag      = params.cd4_csv    ? "CD4_CSV=${params.cd4_csv}"          : ""
    """
    ${vl_flag} ${cd4_flag} virowatch_run.sh \\
        --sample     ${meta.id} \\
        --fastq      ${fastq} \\
        --ref_fa     ${ref_fa} \\
        --ref_gff    ${ref_gff} \\
        --threads    ${task.cpus} \\
        --model      ${params.medaka_model} \\
        --minlen     ${params.chopper_minlen} \\
        --maxlen     ${params.chopper_maxlen} \\
        --minqual    ${params.chopper_q} \\
        --report_dir ${projectDir}/report \\
        ${blast_flag} \\
        ${core_nt_flag}
    """
}

workflow {
    Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> [ [id: row.sample_id], file(row.fastq, checkIfExists: true) ] }
        .set { ch_reads }

    ref_fa  = file(params.ref_fa,  checkIfExists: true)
    ref_gff = file(params.ref_gff, checkIfExists: true)

    VIROWATCH_SAMPLE(ch_reads, ref_fa, ref_gff)
}
