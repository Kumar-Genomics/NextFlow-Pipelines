#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rnaseq-nf: RNA-seq Analysis Pipeline
    FASTQ → QC → Trimming → Alignment → Quantification → Count Matrix
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Author  : Dr. Suhirthakumar Puvanendran
    GitHub  : https://github.com/YOUR_USERNAME/rnaseq-nextflow
    License : MIT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Print pipeline banner
log.info """\
    ╔══════════════════════════════════════════════════════╗
    ║          R N A - S E Q   P I P E L I N E            ║
    ║     FASTQ → QC → Trim → Align → Quantify → Count    ║
    ╚══════════════════════════════════════════════════════╝

    Pipeline Parameters
    ───────────────────────────────────────────────
    reads        : ${params.reads}
    genome       : ${params.genome}
    gtf          : ${params.gtf}
    star_index   : ${params.star_index}
    outdir       : ${params.outdir}
    aligner      : ${params.aligner}
    strandedness : ${params.strandedness}
    ───────────────────────────────────────────────
    """.stripIndent()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RNASEQ } from './workflows/rnaseq'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    RNASEQ ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    log.info """\
        ╔══════════════════════════════════════════════════════╗
        ║              P I P E L I N E   D O N E              ║
        ╚══════════════════════════════════════════════════════╝
        Status     : ${workflow.success ? '✅ SUCCESS' : '❌ FAILED'}
        Duration   : ${workflow.duration}
        Output dir : ${params.outdir}
        Work dir   : ${workflow.workDir}
        """.stripIndent()
}
