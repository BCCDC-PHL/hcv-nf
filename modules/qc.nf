process fastp {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}*.trim.fastq.gz", mode:'copy'

    input:
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    tuple val(sample_id), path("${sample_id}_R1.trim.fastq.gz"), path("${sample_id}_R2.trim.fastq.gz"), emit: trimmed_reads
    tuple val(sample_id),path("${sample_id}_fastp.json"), emit: json

    script:
    """      
    fastp \
      -t ${task.cpus} \
      -i ${reads_1} \
      -I ${reads_2} \
      -o ${sample_id}_R1.trim.fastq.gz \
      -O ${sample_id}_R2.trim.fastq.gz \
      -j ${sample_id}_fastp.json \
      --detect_adapter_for_pe

    """
}

process fastp_json_to_csv {

  tag { sample_id }

  executor 'local'

  //publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}/${sample_id}", pattern: "${sample_id}_fastp.csv", mode: 'copy'
  publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_fastp.csv", mode:'copy'
  
  input:
  tuple val(sample_id), path(fastp_json)

  output:
  tuple val(sample_id), path("${sample_id}_fastp.csv")

  script:
  """
  fastp_json_to_csv.py -s ${sample_id} ${fastp_json} > ${sample_id}_fastp.csv
  """
}


process cutadapter {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}*.out.fastq.gz", mode:'copy'
    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}.cutadapt.log", mode:'copy'

    input:
    tuple val(sample_id), path(reads_1), path(reads_2),path(adapters)

    output:
    tuple val(sample_id), path("${sample_id}_R1.out.fastq.gz"), path("${sample_id}_R2.out.fastq.gz"), emit: out_reads
    path("${sample_id}.cutadapt.log"), emit: log

    script:
    """
    cutadapt \
      -j ${task.cpus} \
      -b file:${adapters} \
      -B file:${adapters} \
      -o ${sample_id}_R1.out.fastq.gz \
      -p ${sample_id}_R2.out.fastq.gz \
      ${reads_1}\
      ${reads_2}\
      > ${sample_id}.cutadapt.log
    """
}

process bbdukadapter {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}*.cleaned.fastq.gz", mode:'copy'

    input:
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    tuple val(sample_id), path("${sample_id}_R1.cleaned.fastq.gz"), path("${sample_id}_R2.cleaned.fastq.gz"), emit: cleaned_reads

    script:
    """
    bbduk.sh in=${reads_1} in2=${reads_2} out=${sample_id}_R1.cleaned.fastq.gz out2=${sample_id}_R2.cleaned.fastq.gz literal=CTGTCTCTTATACACATCT rcomp=t ktrim=rl k=19 mink=11 hdist=1 tpe tbo 
   
    """
}

process errorcorrect {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}*.corrected.fastq.gz", mode:'copy'

    input:
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    tuple val(sample_id), path("${sample_id}_R1.corrected.fastq.gz"), path("${sample_id}_R2.corrected.fastq.gz"), emit: corrected_reads

    script:
    """
    tadpole.sh in=${reads_1} in2=${reads_2} out=${sample_id}_R1.corrected.fastq.gz out2=${sample_id}_R2.corrected.fastq.gz mode=correct
    
    """
}

process normalize {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}*.norm.fastq.gz", mode:'copy'

    input:
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    tuple val(sample_id), path("${sample_id}_R1.norm.fastq.gz"), path("${sample_id}_R2.norm.fastq.gz"), emit: normed_reads

    script:
    """
    bbnorm.sh in=${reads_1} in2=${reads_2} out=${sample_id}_R1.norm.fastq.gz out2=${sample_id}_R2.norm.fastq.gz target=200 min=5
    """
}