#!/usr/bin/env nextflow

import java.time.LocalDateTime

nextflow.enable.dsl = 2

include { fastp } from './modules/qc.nf'
include { fastp_json_to_csv } from './modules/qc.nf'
include { errorcorrect } from './modules/qc.nf'
include { cutadapter } from './modules/qc.nf'
include { normalize } from './modules/qc.nf'
include { bbdukadapter } from './modules/qc.nf'
include { genotype } from './modules/hcv.nf'
include { makeconsensus } from './modules/hcv.nf'
include { maprawreads } from './modules/hcv.nf'
include {mapreadstoref} from './modules/hcv.nf'
include { mixscan } from './modules/mix.nf'
include { findamplicon } from './modules/findamplicon.nf'
include { QualiMap} from './modules/QualiMap.nf'
include { parseQMresults} from './modules/parseQMresults.nf'
include { segcov} from './modules/segcov.nf'
include { mafftraxmltree } from './modules/mafftraxmltree.nf'
include { report } from './modules/report.nf'
// prints to the screen and to the log
        log.info """

                 FluViewer Pipeline
                 ===================================
                 projectDir    : ${projectDir}
                 launchDir     : ${launchDir}
                 mode          : ${params.mode}
                 database      : ${params.db}
                 fastqInputDir : ${params.fastq_input}
                 outdir        : ${params.outdir}
                 """
                 .stripIndent()

workflow{
    ch_adapters = Channel.fromPath(params.adapters)
    ch_db = Channel.fromPath(params.db)
    ch_ref = Channel.fromPath(params.refhcv)
    ch_fastq_input = Channel.fromFilePairs(params.fastq_input + '/*_{R1,R2}*.fastq.gz', flat: true ).map{ it -> [it[0].split('_')[0], it[1], it[2]] }.unique{ it -> it[0] }
    ch_basic_qc = Channel.fromPath(params.basic_qc)
    ch_abundance_top_n = Channel.fromPath(params.abundance_top_n)

    main:

    fastp( ch_fastq_input )
    cutadapter(fastp.out.trimmed_reads.combine(ch_adapters))
    

    maprawreads(cutadapter.out.out_reads.combine(ch_db))
    mapreadstoref(cutadapter.out.out_reads.combine(ch_ref))
    ch_mix = mixscan(mapreadstoref.out.alignment.combine(ch_ref))
    genotype(cutadapter.out.out_reads)
    findamplicon(genotype.out.filtered_contigs)
    ch_consensus = makeconsensus(cutadapter.out.out_reads.combine(findamplicon.out.ref_seqs_mapping, by : 0))
    mafftraxmltree(makeconsensus.out.consensus_seqs)
    QualiMap(makeconsensus.out.alignment)
    ch_qc = parseQMresults(QualiMap.out.genome_results)
    segcov(makeconsensus.out.alignment)

    ch_combined_genotype = ch_consensus.genotyperesult
        .collectFile(it -> it[1], name: "combined_genotype_calls.csv", storeDir: params.outdir, keepHeader: true, skip: 1)

    ch_combined_consensus = ch_consensus.consensus_seqs_report
        .collectFile(it -> it[1], name: "combined_consensus_seqs_report.tsv", storeDir: params.outdir, keepHeader: true, skip: 1)

    ch_combined_demix = ch_mix.demix_report
        .collectFile(it -> it[1], name: "combined_demix_report.csv", storeDir: params.outdir, keepHeader: true, skip: 1)

    ch_combined_qc = ch_qc
        .collectFile(it -> it[1], name: "combined_qc_stats.csv", storeDir: params.outdir, keepHeader: true, skip: 1)
        
    ch_fastqlist = ch_fastq_input
        .collectFile(it -> it[0], name: "fastqlist", storeDir: params.outdir,newLine: true)

    report(ch_fastqlist.combine(ch_combined_genotype).combine(ch_combined_consensus).combine(ch_combined_demix).combine(ch_combined_qc).combine(ch_basic_qc).combine(ch_abundance_top_n))

    
}