process maprawreads {
    
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${sample_id}/debug", pattern: "${sample_id}_mapped_to_db.bam*", mode:'copy'
    publishDir "${params.outdir}/${sample_id}/debug", pattern: "${sample_id}_core_ns5b_mapped_reads.csv", mode:'copy'
    publishDir "${params.outdir}/${sample_id}/debug", pattern: "${sample_id}_mapped_to_db.depth", mode:'copy'

    input:
    tuple val(sample_id), path(reads_1), path(reads_2), path(ref)

    output:
    tuple val(sample_id), path("${sample_id}_mapped_to_db.bam*"), emit: readsbam, optional: true
    tuple val(sample_id), path("${sample_id}_core_ns5b_mapped_reads.csv"), emit: mappedreads, optional: true
    tuple val(sample_id), path("${sample_id}_mapped_to_db.depth"), emit: dbdepth, optional: true

    """
    
    bwa index ${ref}
    bwa mem ${ref} ${reads_1} ${reads_2} > ${sample_id}_align.sam
    samtools view -F 2828 -h ${sample_id}_align.sam | samtools sort -o ${sample_id}_mapped_to_db.bam
    samtools index ${sample_id}_mapped_to_db.bam

    samtools depth ${sample_id}_mapped_to_db.bam > ${sample_id}_mapped_to_db.depth

    samtools view -c -L ${params.corebed} ${sample_id}_mapped_to_db.bam | awk '{print "${sample_id},"\$0}' > core_mapped_reads.csv
    samtools view -c -L ${params.ns5bbed} ${sample_id}_mapped_to_db.bam | awk '{print "${sample_id},"\$0}' > ns5b_mapped_reads.csv

    join -t, -1 1 -2 1 core_mapped_reads.csv ns5b_mapped_reads.csv > ${sample_id}_core_ns5b_mapped_reads.csv
    sed -i -e '1isample_id,core_mapped_reads,ns5b_mapped_reads' ${sample_id}_core_ns5b_mapped_reads.csv
    """

}

process plotdepthdb {
    errorStrategy 'ignore'
    
    publishDir "${params.outdir}/${sample_id}/debug", pattern: "${sample_id}*_db_depth_plots.png", mode:'copy'

    input:
    tuple val(sample_id), path(mapped_db_depth)

    output:
    tuple val(sample_id), path("${sample_id}_*_db_depth_plots.png"), emit: dbdepthplot, optional: true


    """
    
    Rscript ${projectDir}/bin/plotting_db.R ${mapped_db_depth}

    """

}