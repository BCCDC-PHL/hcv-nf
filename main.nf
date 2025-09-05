#!/usr/bin/env nextflow

import java.time.LocalDateTime

nextflow.enable.dsl = 2

include { hash_files }          from './modules/hash_files.nf'
include { pipeline_provenance } from './modules/provenance.nf'
include { collect_provenance }  from './modules/provenance.nf'
include { fastp }               from './modules/qc.nf'
include { blastn_and_filter }   from './modules/hcv.nf'
include { assemble }            from './modules/hcv.nf'
include { cutadapter }          from './modules/qc.nf'
//include { genotype }            from './modules/hcv.nf'
include { makeconsensus }       from './modules/hcv.nf'
include { blastconsensus }      from './modules/hcv.nf'
include { maprawreads }         from './modules/qc.nf'
include { plotdepthdb }         from './modules/qc.nf'
include { mapreadstoref }       from './modules/mix.nf'
include { mixscan }             from './modules/mix.nf'
include { findamplicon }        from './modules/findamplicon.nf'
include { QualiMap }            from './modules/QualiMap.nf'
include { parseQMresults}       from './modules/parseQMresults.nf'
include { segcov}               from './modules/segcov.nf'
include { mafftraxmltree }      from './modules/mafftraxmltree.nf'
include { plot_tree }           from './modules/mafftraxmltree.nf'
include { report }              from './modules/report.nf'

// prints to the screen and to the log
log.info """
         BCCDC-PHL/hcv-nf Pipeline
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
    ch_start_time = Channel.of(LocalDateTime.now())
    ch_pipeline_name = Channel.of(workflow.manifest.name)
    ch_pipeline_version = Channel.of(workflow.manifest.version)

    ch_pipeline_provenance = pipeline_provenance(ch_pipeline_name.combine(ch_pipeline_version).combine(ch_start_time))
    
    ch_adapters = Channel.fromPath(params.adapters)
    //ch_artifacts = Channel.fromPath(params.artifacts)
    ch_db = Channel.fromPath(params.db)
    ch_ref_core = Channel.fromPath(params.ref_core)
    ch_ref_ns5b = Channel.fromPath(params.ref_ns5b)

    ch_ref = Channel.fromPath(params.refhcv)
    ch_fastq_input = Channel.fromFilePairs( params.fastq_search_path, flat: true ).map{ it -> [it[0].split('_')[0], it[1], it[2]] }.unique{ it -> it[0] }
    ch_nt = Channel.fromPath(params.nt_dir)
    ch_db_name = Channel.of(params.db_name)


    main:

    hash_files(ch_fastq_input.map{ it -> [it[0], [it[1], it[2]]] }.combine(Channel.of("fastq_input")))
    //pre_fastqc(ch_fastq_input)
    fastp( ch_fastq_input )
    
    cutadapter(ch_fastq_input.combine(ch_adapters))
    //post_fastqc(cutadapter.out.out_reads)
    //bbdukclean(cutadapter.out.out_reads.combine(ch_artifacts))

    ch_maprawreads = maprawreads(cutadapter.out.out_reads.combine(ch_db)) //mapping raw reads to all 237 HCV references for debugging purposes, checking how reads mapped to core/ns5b regions
                                                        
    mapreadstoref(cutadapter.out.out_reads.combine(ch_ref)) //mapping raw reads to ref 1_AJ851228 for mix variant scan purpose

    plotdepthdb(ch_maprawreads.dbdepth)
    ch_mix = mixscan(mapreadstoref.out.alignment.combine(ch_ref))
    assemble(cutadapter.out.out_reads)  
    blastn_and_filter(assemble.out.contigs.combine(ch_db))
    //genotype(ch_fastq_input.combine(ch_db))

    ch_contigs = blastn_and_filter.out.filtered_contigs.combine(ch_ref_core).combine(ch_ref_ns5b)
    findamplicon(ch_contigs)
    ch_consensus = makeconsensus(cutadapter.out.out_reads.combine(findamplicon.out.ref_seqs_mapping, by : 0).combine(ch_db))
    ch_nt_calls = blastconsensus(ch_consensus.consensus_seqs.combine(ch_nt).combine(ch_db_name))
    mafftraxmltree(makeconsensus.out.consensus_seqs.combine(ch_ref_core).combine(ch_ref_ns5b))
    plot_tree_input = mafftraxmltree.out.core_besttree.mix(mafftraxmltree.out.ns5b_besttree)
    plot_tree(plot_tree_input)
    QualiMap(makeconsensus.out.alignment)
    ch_qc = parseQMresults(QualiMap.out.genome_results)
    segcov(makeconsensus.out.alignment)

    ch_combined_genotype = ch_nt_calls.genotyperesult
        .collectFile(it -> it[1], name: "combined_genotype_calls.csv", storeDir: params.outdir, keepHeader: true, skip: 1)

    ch_combined_consensus = ch_consensus.consensus_seqs_report
        .collectFile(it -> it[1], name: "combined_consensus_seqs_report.tsv", storeDir: params.outdir, keepHeader: true, skip: 1)

    ch_combined_demix = ch_mix.demix_report
        .collectFile(it -> it[1], name: "combined_demix_report.csv", storeDir: params.outdir, keepHeader: true, skip: 1)

    ch_combined_qc = ch_qc
        .collectFile(it -> it[1], name: "combined_qc_stats.csv", storeDir: params.outdir, keepHeader: true, skip: 1)
      
    ch_fastqlist = ch_fastq_input
        .collectFile(it -> it[0], name: "fastqlist", storeDir: params.outdir,newLine: true)
    ch_count_mapped_reads = ch_maprawreads.mappedreads
        .collectFile(it -> it[1], name: "combined_amplicon_mapped_reads_counts.csv", storeDir:params.outdir, keepHeader: true, skip: 1)

    report(ch_fastqlist.combine(ch_combined_consensus).combine(ch_combined_demix).combine(ch_combined_qc).combine(ch_count_mapped_reads).combine(ch_combined_genotype))

    
    ch_provenance = assemble.out.provenance
    ch_provenance = ch_provenance.join(fastp.out.provenance).map{ it -> [it[0], [it[1] , it[2]]] }
    ch_provenance = ch_provenance.join(hash_files.out.provenance).map{ it -> [it[0], it[1]<< it[2]] }
    ch_provenance = ch_provenance.join(blastn_and_filter.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(ch_nt_calls.provenance).map{ it -> [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(ch_mix.provenance).map{ it -> [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(ch_fastq_input.map{ it -> it[0] }.combine(ch_pipeline_provenance)).map{ it -> [it[0], it[1] << it[2]] }
    collect_provenance(ch_provenance)
}