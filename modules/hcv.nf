
process genotype {
    errorStrategy 'ignore'

    memory { 55.GB }
    tag {sample_id}

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}/${sample_id}*", mode:'copy', saveAs: { filename -> filename.split("/").last() }
    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}/logs", mode:'copy', saveAs: { filename -> "f-f logs" }


    input:
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    //tuple val(sample_id), path("${sample_id}/${sample_id}*.bam*"), emit: alignment, optional: true
    //tuple val(sample_id), path("${sample_id}/${sample_id}*_report.tsv"), emit: reports, optional: true
    tuple val(sample_id), path("${sample_id}/${sample_id}_blast_results_prefilter.csv"), emit: allblastresult, optional: true
    tuple val(sample_id), path("${sample_id}/${sample_id}_filtered_blast_results.csv"), emit: blastreport, optional: true
    //tuple val(sample_id), path("${sample_id}/${sample_id}_genotype_calls.csv"), emit: genotyperesult, optional: true
    tuple val(sample_id), path("${sample_id}/${sample_id}_contigs.fa"), emit: contigs, optional: true
    tuple val(sample_id), path("${sample_id}/${sample_id}_filtered_contigs.fa"), emit: filtered_contigs, optional: true
    

    //tuple val(sample_id), path("${sample_id}/${sample_id}*_ref_seqs_for_mapping.fa"), emit: refseq, optional: true
    //tuple val(sample_id), path("${sample_id}/${sample_id}_max_bitscores_per_subtype.csv"), emit: mean_bitscores_per_subtype, optional: true
    //tuple val(sample_id), path("${sample_id}/${sample_id}_consensus_blast_results.tsv"), emit: consensusmapping, optional: true
    //tuple val(sample_id), path("${sample_id}/${sample_id}_consensus_seqs_report.tsv"), emit: consensus_seqs_report, optional: true
    //tuple val(sample_id), path("${sample_id}/${sample_id}_consensus_seqs.fa"), emit: consensus_seqs, optional: true
    tuple val(sample_id), path("${sample_id}/logs"), emit: fluviewer_logs

    """
    genotype_v1.py -f ${reads_1} -r ${reads_2} -o ${sample_id} -d ${params.db} -m ${params.mode}

    """

}

process makeconsensus {

    memory { 55.GB }
    tag {sample_id}

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}/${sample_id}*", mode:'copy', saveAs: { filename -> filename.split("/").last() }
    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}/logs", mode:'copy', saveAs: { filename -> "f-f logs" }


    input:
    tuple val(sample_id), path(reads_1), path(reads_2), path(ref_seq_map)

    output:
    tuple val(sample_id), path("${sample_id}/${sample_id}*.bam*"), emit: alignment, optional: true
    tuple val(sample_id), path("${sample_id}/${sample_id}*_report.tsv"), emit: reports, optional: true
    //tuple val(sample_id), path("${sample_id}/${sample_id}_blast_results_prefilter.csv"), emit: allblastresult, optional: true
    //tuple val(sample_id), path("${sample_id}/${sample_id}_filtered_blast_results.csv"), emit: blastreport, optional: true
    tuple val(sample_id), path("${sample_id}/${sample_id}_genotype_calls.csv"), emit: genotyperesult, optional: true
    //tuple val(sample_id), path("${sample_id}/${sample_id}_contigs.fa"), emit: contigs, optional: true
    //tuple val(sample_id), path("${sample_id}/${sample_id}_filtered_contigs.fa"), emit: filtered_contigs, optional: true
    

    //tuple val(sample_id), path("${sample_id}/${sample_id}*_ref_seqs_for_mapping.fa"), emit: refseq, optional: true
    tuple val(sample_id), path("${sample_id}/${sample_id}_max_bitscores_per_subtype.csv"), emit: mean_bitscores_per_subtype, optional: true
    tuple val(sample_id), path("${sample_id}/${sample_id}_consensus_blast_results.tsv"), emit: consensusmapping, optional: true
    tuple val(sample_id), path("${sample_id}/${sample_id}_consensus_seqs_report.tsv"), emit: consensus_seqs_report, optional: true
    tuple val(sample_id), path("${sample_id}/${sample_id}_consensus_seqs.fa"), emit: consensus_seqs, optional: true
    //tuple val(sample_id), path("${sample_id}/${sample_id}_consensus_seqs_noprimers.fa"), emit: consensus_seqsnoprimers, optional: true
    tuple val(sample_id), path("${sample_id}/logs"), emit: fluviewer_logs

    """
    genotype_v1_2.py -f ${reads_1} -r ${reads_2} -s ${ref_seq_map} -o ${sample_id} -d ${params.db} -m ${params.mode}

    """

}
