process parseQMresults {

  tag { sample_id }

  publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_parsed_genome_results.csv", mode:'copy'

  input:
  tuple val(sample_id), path(qualimap_result_txt)

  output:
  tuple val(sample_id), path("${sample_id}_parsed_genome_results.csv")

  //val(sample_id), path("${sample_id}"_FluViewer_provenance.yml), emit: provenance

  script:
  """

  parse_qualimap_results.py ${qualimap_result_txt} -o ${sample_id}_parsed_genome_results.csv
  """
}