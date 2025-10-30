process fastp {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}*.trim.fastq.gz", mode:'copy'
    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_fastp.*", mode:'copy'

    input:
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    tuple val(sample_id), path("${sample_id}_trim_R1.fastq.gz"), path("${sample_id}_trim_R2.fastq.gz"), emit: trimmed_reads
    //tuple val(sample_id), path("${sample_id}_fastp.json"), emit: json
    tuple val(sample_id), path("${sample_id}_fastp_provenance.yml"), emit: provenance
    tuple val(sample_id), path("${sample_id}_fastp.csv")            , emit: metrics
    tuple val(sample_id), path("${sample_id}_fastp.json")           , emit: report_json
    tuple val(sample_id), path("${sample_id}_fastp.html")           , emit: report_html

    script:
    """      
    printf -- "- process_name: fastp\\n" > ${sample_id}_fastp_provenance.yml
    printf -- "  tool_name: fastp\\n  tool_version: \$(fastp --version 2>&1 | cut -d ' ' -f 2)\\n" >> ${sample_id}_fastp_provenance.yml

    fastp \
      -t ${task.cpus} \
      -i ${reads_1} \
      -I ${reads_2} \
      -o ${sample_id}_trim_R1.fastq.gz \
      -O ${sample_id}_trim_R2.fastq.gz \
      --detect_adapter_for_pe \
      --trim_poly_g \
      --overrepresentation_analysis \
      --report_title "fastp report: ${sample_id}" \
      --json ${sample_id}_fastp.json \
      --html ${sample_id}_fastp.html


    fastp_json_to_csv.py -s ${sample_id} ${sample_id}_fastp.json > ${sample_id}_fastp.csv
    """
}

process pre_fastqc {
    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}*_fastqc.html", mode:'copy'

    input:
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    tuple val(sample_id), path("${sample_id}*_fastqc.html")           , emit: fastqc_report_html

    """
    fastqc --extract ${reads_1} ${reads_2} 

    """
}

process post_fastqc {
    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}*_fastqc.html", mode:'copy'

    input:
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    tuple val(sample_id), path("${sample_id}*_fastqc.html")           , emit: fastqc_report_html

    """
    fastqc --extract ${reads_1} ${reads_2} 

    """
}

process cutadapter {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}*.out.fastq.gz", mode:'copy'
    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}.cutadapt.log", mode:'copy'

    input:
    tuple val(sample_id), path(reads_1), path(reads_2),path(adapters)

    output:
    tuple val(sample_id), path("${sample_id}_out_R1.fastq.gz"), path("${sample_id}_out_R2.fastq.gz"), emit: out_reads
    path("${sample_id}.cutadapt.log"), emit: log
    tuple val(sample_id), path("${sample_id}_cutadapt_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: cutadapt\\n" > ${sample_id}_cutadapt_provenance.yml
    printf -- "  tool_name: cutadapt\\n  tool_version: \$(cutadapt --version 2>&1 | cut -d ' ' -f 2)\\n" >> ${sample_id}_cutadapt_provenance.yml

    cutadapt \
      -j ${task.cpus} \
      --nextseq-trim=20 -m 1 \
      -b file:${adapters} \
      -B file:${adapters} \
      -o ${sample_id}_out_R1.fastq.gz \
      -p ${sample_id}_out_R2.fastq.gz \
      ${reads_1}\
      ${reads_2}\
      > ${sample_id}.cutadapt.log
      
    """
}

process bbdukclean {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}*.cleaned.fastq.gz", mode:'copy'

    input:
    tuple val(sample_id), path(reads_1), path(reads_2), path(artifacts)

    output:
    tuple val(sample_id), path("${sample_id}_R1.cleaned.fastq.gz"), path("${sample_id}_R2.cleaned.fastq.gz"), emit: cleaned_reads

    script:
    """
    bbduk.sh in=${reads_1} in2=${reads_2} out=${sample_id}_R1.cleaned.fastq.gz out2=${sample_id}_R2.cleaned.fastq.gz ref=adapter,artifacts tbo tpe hdist=1 ktrim=r mink=11 qtrim=rl trimq=10 trimpolyg=10 entropy=0.7
   
    """
}

process maprawreads {
    
    errorStrategy 'ignore'
    //publishDir "${params.outdir}/${sample_id}/debug", pattern: "${sample_id}_mapped_to_db.bam*", mode:'copy'
    //publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_core_ns5b_mapped_reads.csv", mode:'copy'
    //publishDir "${params.outdir}/${sample_id}/debug", pattern: "${sample_id}_mapped_to_db.depth", mode:'copy'

    input:
    tuple val(sample_id), path(reads_1), path(reads_2), path(ref)

    output:
    //tuple val(sample_id), path("${sample_id}_mapped_to_db.bam*"), emit: readsbam, optional: true
    tuple val(sample_id), path("${sample_id}_core_ns5b_mapped_reads.csv"), emit: mappedreads, optional: true
    tuple val(sample_id), path("${sample_id}_mapped_to_db.depth"), emit: dbdepth, optional: true

    """
    
    bwa index ${ref}
    bwa mem ${ref} ${reads_1} ${reads_2} > ${sample_id}_align.sam
    samtools view -f 1 -F 2316 -h ${sample_id}_align.sam | samtools sort -o ${sample_id}_mapped_to_db.bam
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
    
    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_*_db_depth_plots.png", mode:'copy'

    input:
    tuple val(sample_id), path(mapped_db_depth)

    output:
    tuple val(sample_id), path("${sample_id}_*_db_depth_plots.png"), emit: dbdepthplot, optional: true


    """
    
    Rscript ${projectDir}/bin/plotting_db.R ${mapped_db_depth}

    """

}
