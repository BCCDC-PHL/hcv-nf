process mapreadstoref {
    

    //publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_mapped_to_ref.bam*", mode:'copy'

    input:
    tuple val(sample_id), path(reads_1), path(reads_2), path(ref)

    output:
    tuple val(sample_id), path("${sample_id}_mapped_to_ref.bam"), emit: alignment, optional: true

    """
    
    bwa index ${ref}
    bwa mem ${ref} ${reads_1} ${reads_2} > ${sample_id}_align.sam
    samtools view -f 1 -F 2316 -q 30 -h ${sample_id}_align.sam | samtools sort -o ${sample_id}_mapped_to_ref.bam
    samtools index ${sample_id}_mapped_to_ref.bam

    """

}

process mixscan {
  errorStrategy 'ignore'
    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}/demix", pattern: "${sample_id}*", mode:'copy'

    input:
    tuple val(sample_id), path(alignment_bam), path(ref)

    output:
    tuple val(sample_id), path("${sample_id}_variants.tsv"), emit: snp_freyja, optional: true
    tuple val(sample_id), path("${sample_id}_variants_depths"), emit: snp_depth, optional: true
    tuple val(sample_id), path("${sample_id}_demixing_results.tsv"), emit: demix_results, optional: true
    tuple val(sample_id), path("${sample_id}_demix.csv"), emit: demix_report, optional: true
    tuple val(sample_id), path("${sample_id}_mixscan_provenance.yml"), emit: provenance

    script:
    """      
    printf -- "- process_name: mixscan\\n" > ${sample_id}_mixscan_provenance.yml
    printf -- "  tool_name: Freyja\\n  tool_version: \$(freyja --version | cut -d' ' -f3)\\n" >> ${sample_id}_mixscan_provenance.yml
    printf -- "  reference used: ${ref}\\n" >> ${sample_id}_mixscan_provenance.yml

    freyja variants ${alignment_bam} --variants ${sample_id}_variants --depths ${sample_id}_variants_depths --ref ${ref}
    demix_freyja_adapted.py ${sample_id}_variants.tsv ${sample_id}_variants_depths --output ${sample_id}_demixing_results.tsv --sample ${sample_id}
    """
}