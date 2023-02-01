process mafftraxmltree {
  errorStrategy 'ignore'

  tag { sample_id }

  publishDir "${params.outdir}/${sample_id}", pattern: "RAxML*", mode:'copy', saveAs: { filename -> filename.split("/").last() }
  publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_mafftouput_*", mode:'copy'


  input:
  tuple val(sample_id), path(consensus)

  output:
  tuple val(sample_id), path("RAxML_bestTree.${sample_id}_ns5b"), emit: ns5b_besttree, optional: true
  tuple val(sample_id), path("RAxML_bestTree.${sample_id}_core"), emit: core_besttree, optional: true
  tuple val(sample_id), path("${sample_id}_mafftoutput_core"), emit: core_alignment, optional: true
  tuple val(sample_id), path("${sample_id}_mafftoutput_ns5b"), emit: ns5b_alignment, optional: true
  tuple val(sample_id), path("${sample_id}_core_consensus.fa"), emit: core_consensus, optional: true
  tuple val(sample_id), path("${sample_id}_ns5b_consensus.fa"), emit: ns5b_consensus, optional: true

  script:
  """
  grep -A1 '|core|' ${consensus} > ${sample_id}_core_consensus.fa

  grep -A1 '|ns5b|' ${consensus} > ${sample_id}_ns5b_consensus.fa

  if [ -s ${sample_id}_core_consensus.fa ]; then
    cat ${params.ref_core} ${sample_id}_core_consensus.fa > mafftinput_core.fa
    mafft --reorder --adjustdirection --anysymbol --thread 4 --auto mafftinput_core.fa > ${sample_id}_mafftoutput_core
    raxmlHPC -f d -p 12345 -s ${sample_id}_mafftoutput_core -m GTRCAT -n ${sample_id}_core
  fi

  if [ -s ${sample_id}_ns5b_consensus.fa ]; then
    cat ${params.ref_ns5b} ${sample_id}_ns5b_consensus.fa > mafftinput_ns5b.fa
    mafft --reorder --adjustdirection --anysymbol --thread 4 --auto mafftinput_ns5b.fa > ${sample_id}_mafftoutput_ns5b
    raxmlHPC -f d -p 12345 -s ${sample_id}_mafftoutput_ns5b -m GTRCAT -n ${sample_id}_ns5b
  fi

  """
}