process segcov { 
  errorStrategy 'ignore'
  tag { sample_id }

  publishDir "${params.outdir}/${sample_id}", mode:'copy'


  input:
  tuple val(sample_id), path(bamfile)

  output:
  tuple val(sample_id), path("*_alignment.coverage"),emit: coverage_file
  tuple val(sample_id), path("*_depth_plots.pdf"), emit: coverage_plot

  """
  samtools depth -a ${bamfile} > ${sample_id}_alignment.coverage
  Rscript ${projectDir}/bin/plotting.R ${sample_id}_alignment.coverage
  """


}