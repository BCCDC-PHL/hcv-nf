process mixscan {
  errorStrategy 'ignore'
    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}*", mode:'copy'

    input:
    tuple val(sample_id), path(alignment_bam), path(ref)

    output:
    tuple val(sample_id), path("${sample_id}_variants.tsv"), emit: snp_freyja, optional: true
    tuple val(sample_id), path("${sample_id}_variants_depths"), emit: snp_depth, optional: true
    tuple val(sample_id), path("${sample_id}_demixing_results.tsv"), emit: demix_results, optional: true
    tuple val(sample_id), path("${sample_id}_demix.csv"), emit: demix_report, optional: true

    script:
    """      
    freyja variants ${alignment_bam} --variants ${sample_id}_variants --depths ${sample_id}_variants_depths --ref ${ref}
    demix_freyja_adapted.py ${sample_id}_variants.tsv ${sample_id}_variants_depths --output ${sample_id}_demixing_results.tsv --sample ${sample_id}
    """
}