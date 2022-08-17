process QualiMap {
  errorStrategy 'ignore'

  tag { sample_id }

  publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_alignment_filtered_sorted_stats", mode:'copy', saveAs: { filename -> filename.split("/").last() }

  input:
  tuple val(sample_id), path(bamfile)

  output:
  tuple val(sample_id), path("${sample_id}_alignment_filtered_sorted_stats/genome_results.txt"), emit: genome_results, optional:true
  tuple val(sample_id), path("${sample_id}_amplicon.bed"), emit: amplicons, optional: true
  tuple val(sample_id), path("${sample_id}_alignment_filtered_sorted_stats")
  

  //val(sample_id), path("${sample_id}"_FluViewer_provenance.yml), emit: provenance
  //samtools view ${bamfile} | cut -f3 | sort -u > genotype_name.txt
  //associate_bed_to_bam.py --bed ${params.bedfile} --name genotype_name.txt
  //qualimap bamqc -bam ${bamfile} -gff ${sample_id}_amplicon.bed

  script:
  """
  samtools view ${bamfile} | cut -f3 | sort -u > genotype_name.txt
  associate_bed_to_bam.py --bed ${params.bedfile} --name genotype_name.txt
  qualimap bamqc -bam ${bamfile} -gff ${sample_id}_amplicon.bed
  """
}