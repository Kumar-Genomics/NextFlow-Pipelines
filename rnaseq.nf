/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RNASEQ WORKFLOW
    Orchestrates: Input validation → QC → Trimming → Indexing → Alignment → Quantification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

// ── Import modules ─────────────────────────────────────────────────────────────────

include { FASTQC              } from '../modules/local/fastqc'
include { TRIMMOMATIC         } from '../modules/local/trimmomatic'
include { FASTP               } from '../modules/local/fastp'
include { STAR_GENOMEGENERATE } from '../modules/local/star_genomegenerate'
include { STAR_ALIGN          } from '../modules/local/star_align'
include { HISAT2_BUILD        } from '../modules/local/hisat2_build'
include { HISAT2_ALIGN        } from '../modules/local/hisat2_align'
include { SALMON_INDEX        } from '../modules/local/salmon_index'
include { SALMON_QUANT        } from '../modules/local/salmon_quant'
include { SAMTOOLS_SORT       } from '../modules/local/samtools_sort'
include { SAMTOOLS_INDEX      } from '../modules/local/samtools_index'
include { SAMTOOLS_FLAGSTAT   } from '../modules/local/samtools_flagstat'
include { FEATURECOUNTS       } from '../modules/local/featurecounts'
include { HTSEQ_COUNT         } from '../modules/local/htseq_count'
include { MERGE_COUNTS        } from '../modules/local/merge_counts'
include { INFER_STRANDEDNESS  } from '../modules/local/infer_strandedness'
include { MULTIQC             } from '../modules/local/multiqc'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW: RNASEQ
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNASEQ {

    // ── Channel: paired-end reads from samplesheet or glob ────────────────────────
    if (params.reads) {
        if (params.single_end) {
            Channel
                .fromPath(params.reads)
                .map { file -> tuple(file.simpleName, file) }
                .set { ch_reads }
        } else {
            Channel
                .fromFilePairs(params.reads, checkIfExists: true)
                .set { ch_reads }
        }
    } else {
        error "No reads specified. Use --reads 'path/to/*_{1,2}.fastq.gz'"
    }

    // ── Channel: Reference files ──────────────────────────────────────────────────
    ch_genome     = params.genome     ? Channel.fromPath(params.genome,     checkIfExists: true) : Channel.empty()
    ch_gtf        = params.gtf        ? Channel.fromPath(params.gtf,        checkIfExists: true) : Channel.empty()
    ch_star_index = params.star_index ? Channel.fromPath(params.star_index, checkIfExists: true) : Channel.empty()

    // Accumulate QC reports for MultiQC
    ch_multiqc_files = Channel.empty()

    // ── STEP 1: FastQC (raw reads) ────────────────────────────────────────────────
    if (!params.skip_fastqc) {
        FASTQC(ch_reads)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect())
    }

    // ── STEP 2: Read Trimming ─────────────────────────────────────────────────────
    if (!params.skip_trimming) {
        if (params.trimmer == 'fastp') {
            FASTP(ch_reads)
            ch_trimmed       = FASTP.out.reads
            ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect())
        } else {
            TRIMMOMATIC(ch_reads)
            ch_trimmed       = TRIMMOMATIC.out.reads
            ch_multiqc_files = ch_multiqc_files.mix(TRIMMOMATIC.out.log.collect())
        }
    } else {
        ch_trimmed = ch_reads
    }

    // ── STEP 3: Post-trim FastQC ──────────────────────────────────────────────────
    if (!params.skip_fastqc && !params.skip_trimming) {
        FASTQC(ch_trimmed.map { meta, reads -> tuple("${meta}_trimmed", reads) })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect())
    }

    // ── STEP 4: Genome Indexing ───────────────────────────────────────────────────
    if (params.aligner == 'star') {
        if (!params.star_index) {
            STAR_GENOMEGENERATE(ch_genome, ch_gtf)
            ch_star_index = STAR_GENOMEGENERATE.out.index
        }
    } else if (params.aligner == 'hisat2') {
        if (!params.hisat2_index) {
            HISAT2_BUILD(ch_genome)
            ch_hisat2_index = HISAT2_BUILD.out.index
        } else {
            ch_hisat2_index = Channel.fromPath(params.hisat2_index, checkIfExists: true)
        }
    } else if (params.aligner == 'salmon') {
        if (!params.salmon_index) {
            SALMON_INDEX(ch_genome, ch_gtf)
            ch_salmon_index = SALMON_INDEX.out.index
        } else {
            ch_salmon_index = Channel.fromPath(params.salmon_index, checkIfExists: true)
        }
    }

    // ── STEP 5: Alignment / Quasi-mapping ─────────────────────────────────────────
    if (params.aligner == 'star') {
        STAR_ALIGN(ch_trimmed, ch_star_index, ch_gtf)
        ch_bam           = STAR_ALIGN.out.bam
        ch_multiqc_files = ch_multiqc_files.mix(STAR_ALIGN.out.log_final.collect())

    } else if (params.aligner == 'hisat2') {
        HISAT2_ALIGN(ch_trimmed, ch_hisat2_index, ch_gtf)
        ch_bam           = HISAT2_ALIGN.out.bam
        ch_multiqc_files = ch_multiqc_files.mix(HISAT2_ALIGN.out.summary.collect())

    } else if (params.aligner == 'salmon') {
        SALMON_QUANT(ch_trimmed, ch_salmon_index, ch_gtf)
        ch_multiqc_files = ch_multiqc_files.mix(SALMON_QUANT.out.results.collect())
    }

    // ── STEP 6: BAM Sorting, Indexing, Flagstat (genome aligners only) ────────────
    if (params.aligner != 'salmon') {
        SAMTOOLS_SORT(ch_bam)
        SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)
        SAMTOOLS_FLAGSTAT(SAMTOOLS_SORT.out.bam)
        ch_sorted_bam    = SAMTOOLS_SORT.out.bam
        ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_FLAGSTAT.out.stats.collect())

        // ── STEP 7: Strandedness Inference (optional auto mode) ──────────────────
        if (params.strandedness == 'auto') {
            INFER_STRANDEDNESS(ch_sorted_bam, ch_gtf)
            ch_strandedness = INFER_STRANDEDNESS.out.strandedness
        } else {
            ch_strandedness = Channel.value(params.strandedness)
        }

        // ── STEP 8: Read Quantification ──────────────────────────────────────────
        if (params.quantifier == 'featurecounts') {
            FEATURECOUNTS(ch_sorted_bam, ch_gtf)
            ch_counts        = FEATURECOUNTS.out.counts
            ch_multiqc_files = ch_multiqc_files.mix(FEATURECOUNTS.out.summary.collect())

        } else if (params.quantifier == 'htseq') {
            HTSEQ_COUNT(ch_sorted_bam, ch_gtf)
            ch_counts = HTSEQ_COUNT.out.counts
        }

        // ── STEP 9: Merge count files into a single matrix ───────────────────────
        MERGE_COUNTS(ch_counts.collect())
    }

    // ── STEP 10: MultiQC Report ───────────────────────────────────────────────────
    if (!params.skip_multiqc) {
        MULTIQC(ch_multiqc_files.collect())
    }
}
