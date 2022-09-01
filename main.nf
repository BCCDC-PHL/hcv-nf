#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { fastp } from './modules/fastp.nf'
include { fastp_json_to_csv } from './modules/fastp.nf'
include { genotype } from './modules/hcv.nf'
include { QualiMap} from './modules/QualiMap.nf'
include { parseQMresults} from './modules/parseQMresults.nf'
include { segcov} from './modules/segcov.nf'
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

    fastp( ch_fastq_input.combine(ch_adapters))
    genotype(fastp.out.trimmed_reads)
    QualiMap(genotype.out.alignment)
    parseQMresults(QualiMap.out.genome_results)
    segcov(genotype.out.alignment)
}