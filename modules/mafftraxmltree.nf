process mafftraxmltree {
  errorStrategy 'ignore'

  tag { sample_id }

  publishDir "${params.outdir}/${sample_id}", pattern: "RAxML*", mode:'copy', saveAs: { filename -> filename.split("/").last() }
  //publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_mafftouput*", mode:'copy'


  input:
  tuple val(sample_id), path(consensus), path(ref_core),path(ref_ns5b)

  output:
  tuple val(sample_id), path("RAxML_bestTree.${sample_id}_ns5b"), emit: ns5b_besttree, optional: true
  tuple val(sample_id), path("RAxML_bestTree.${sample_id}_core"), emit: core_besttree, optional: true
  //tuple val(sample_id), path("${sample_id}_mafftoutput*"), emit: msa, optional: true
  //tuple val(sample_id), path("${sample_id}_mafftoutput_ns5b"), emit: ns5b_alignment, optional: true

  script:
  """
  grep -A1 '|core|' ${consensus} > ${sample_id}_core_consensus.fa

  grep -A1 '|ns5b|' ${consensus} > ${sample_id}_ns5b_consensus.fa

  if [ -s ${sample_id}_core_consensus.fa ]; then
    cat ${ref_core} ${sample_id}_core_consensus.fa > mafftinput_core.fa
    mafft --reorder --adjustdirection --anysymbol --thread 4 --auto mafftinput_core.fa > ${sample_id}_mafftoutput_core
    raxmlHPC -f d -p 12345 -# 10 -s ${sample_id}_mafftoutput_core -m GTRCAT -n ${sample_id}_core -o 7_KU861171
  fi

  if [ -s ${sample_id}_ns5b_consensus.fa ]; then
    cat ${ref_ns5b} ${sample_id}_ns5b_consensus.fa > mafftinput_ns5b.fa
    mafft --reorder --adjustdirection --anysymbol --thread 4 --auto mafftinput_ns5b.fa > ${sample_id}_mafftoutput_ns5b
    raxmlHPC -f d -p 12345 -# 10 -s ${sample_id}_mafftoutput_ns5b -m GTRCAT -n ${sample_id}_ns5b -o 7_KU861171
  fi

  """
}

process plot_tree {
    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "*.png"

    input:
    tuple val(sample_id), path(tree_file)

    output:
    tuple val(sample_id), path("*.png"), emit: treeplots

    script:
    // Extract tree type from filename (core or ns5b)
    def tree_type = tree_file.name.contains("_core") ? "core" :
                    tree_file.name.contains("_ns5b") ? "ns5b" : "unknown"

    """
    Rscript ${projectDir}/bin/plot_tree.R --sample_id ${sample_id} --tree_file ${tree_file} --tree_type ${tree_type}
    """
}
