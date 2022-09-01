process fastp {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}*.trim.fastq.gz", mode:'copy'

    input:
    tuple val(sample_id), path(reads_1), path(reads_2), path(adapters)

    output:
    tuple val(sample_id), path("${sample_id}_R1.trim.fastq.gz"), path("${sample_id}_R2.trim.fastq.gz"), emit: trimmed_reads
    path("${sample_id}.fastp.json"), emit: json

    script:
    """
    fastp \
      -t ${task.cpus} \
      -i ${reads_1} \
      -I ${reads_2} \
      -o ${sample_id}_R1.trim.fastq.gz \
      -O ${sample_id}_R2.trim.fastq.gz \
      -j ${sample_id}.fastp.json \
      --detect_adapter_for_pe \
      --adapter_fasta ${adapters}
    """
}

//removed from script - put back in for provenance
//printf -- "- process_name: fastp\\n" > ${sample_id}_fastp_provenance.yml
//printf -- "  tool_name: fastp\\n  tool_version: \$(fastp --version 2>&1 | cut -d ' ' -f 2)\\n" >> ${sample_id}_fastp_provenance.yml

process fastp_json_to_csv {

  tag { sample_id }

  executor 'local'

  publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}/${sample_id}", pattern: "${sample_id}_fastp.csv", mode: 'copy'

  input:
  tuple val(sample_id), path(fastp_json)

  output:
  tuple val(sample_id), path("${sample_id}_fastp.csv")

  script:
  """
  fastp_json_to_csv.py -s ${sample_id} ${fastp_json} > ${sample_id}_fastp.csv
  """
}