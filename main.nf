#!/usr/bin/env nextflow

import java.time.LocalDateTime

nextflow.enable.dsl = 2

include { fastp } from './modules/qc.nf'
include { fastp_json_to_csv } from './modules/qc.nf'
include { errorcorrect } from './modules/qc.nf'
include { cutadapter } from './modules/qc.nf'
include { bbdukadapter } from './modules/qc.nf'
include { genotype } from './modules/hcv.nf'
include { makeconsensus } from './modules/hcv.nf'
include { findamplicon } from './modules/findamplicon.nf'
include { QualiMap} from './modules/QualiMap.nf'
include { parseQMresults} from './modules/parseQMresults.nf'
include { segcov} from './modules/segcov.nf'
include { mafftraxmltree } from './modules/mafftraxmltree.nf'
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
    ch_fastq_input = Channel.fromFilePairs(params.fastq_input + '/*_{R1,R2}*.fastq.gz', flat: true ).map{ it -> [it[0].split('_')[0], it[1], it[2]] }.unique{ it -> it[0] }
    
    //ch_fastq_input.view()

    main:

    //fastp(ch_fastq_input.combine(ch_adapters))
    //fastp_json_to_csv(fastp.out.json)
    cutadapter(ch_fastq_input)
    bbdukadapter(cutadapter.out.out_reads)
  
    genotype(bbdukadapter.out.cleaned_reads)
    findamplicon(genotype.out.filtered_contigs)
    //ch_fastq_input.combine(findamplicon.out.ref_seqs_mapping, by : 0).view()
    makeconsensus(bbdukadapter.out.cleaned_reads.combine(findamplicon.out.ref_seqs_mapping, by : 0))
    //mafftraxmltree(genotype.out.consensus_seqs)
    QualiMap(makeconsensus.out.alignment)
    parseQMresults(QualiMap.out.genome_results)
    segcov(makeconsensus.out.alignment)
}