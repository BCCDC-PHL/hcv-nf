process mafftraxmltree {
  errorStrategy 'ignore'

  tag { sample_id }

  publishDir "${params.outdir}/${sample_id}", pattern: "RAxML*", mode:'copy', saveAs: { filename -> filename.split("/").last() }
  //publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_mafftouput*", mode:'copy'


  input:
  tuple val(sample_id), path(consensus), path(ref_core),path(ref_ns5b),path(rep_strains)

  output:
  tuple val(sample_id), path("RAxML_bestTree.${sample_id}_ns5b*"), emit: ns5b_besttree, optional: true
  tuple val(sample_id), path("RAxML_bestTree.${sample_id}_core*"), emit: core_besttree, optional: true
  //tuple val(sample_id), path("${sample_id}_mafftoutput*"), emit: msa, optional: true
  //tuple val(sample_id), path("${sample_id}_mafftoutput_ns5b"), emit: ns5b_alignment, optional: true

  script:
  """
  grep -A1 '|core|' ${consensus} > ${sample_id}_core_consensus.fa
  coregeno=\$(grep '^>' ${sample_id}_core_consensus.fa | cut -d'|' -f4 | cut -c 1)
  grep -A1 '|ns5b|' ${consensus} > ${sample_id}_ns5b_consensus.fa
  ns5bgeno=\$(grep '^>' ${sample_id}_ns5b_consensus.fa | cut -d'|' -f4 | cut -c 1)

  if [ -s ${sample_id}_core_consensus.fa ]; then
    grep -A 1 -f ${rep_strains} ${ref_core} | grep -v "^--" > reduced_core_ref.fa
    cat reduced_core_ref.fa ${sample_id}_core_consensus.fa > mafftinput_core.fa
    mafft --reorder --adjustdirection --anysymbol --thread 4 --auto mafftinput_core.fa > ${sample_id}_mafftoutput_core
    raxmlHPC -f d -p 12345 -# 10 -s ${sample_id}_mafftoutput_core -m GTRCAT -n ${sample_id}_core -o 7_KU861171

    awk -v g="\$coregeno" '/^>/ {f=(\$0 ~ "^>"g)} f' "${ref_core}" > subtype_core_ref.fa
    cat subtype_core_ref.fa ${sample_id}_core_consensus.fa > mafftinput_core_subtype.fa
    if [ \$coregeno -ne "7" ];then
      grep -A 1 ">7_KU861171" ${ref_core} >> mafftinput_core_subtype.fa
    fi
    mafft --reorder --adjustdirection --anysymbol --thread 4 --auto mafftinput_core_subtype.fa > ${sample_id}_mafftoutput_core_subtype
    raxmlHPC -f d -p 12345 -# 10 -s ${sample_id}_mafftoutput_core_subtype -m GTRCAT -n ${sample_id}_core_subtype -o 7_KU861171

  fi

  if [ -s ${sample_id}_ns5b_consensus.fa ]; then
    grep -A 1 -f ${rep_strains} ${ref_ns5b} | grep -v "^--" > reduced_ns5b_ref.fa
    cat reduced_ns5b_ref.fa ${sample_id}_ns5b_consensus.fa > mafftinput_ns5b.fa
    mafft --reorder --adjustdirection --anysymbol --thread 4 --auto mafftinput_ns5b.fa > ${sample_id}_mafftoutput_ns5b
    raxmlHPC -f d -p 12345 -# 10 -s ${sample_id}_mafftoutput_ns5b -m GTRCAT -n ${sample_id}_ns5b -o 7_KU861171
  
    awk -v g="\$ns5bgeno" '/^>/ {f=(\$0 ~ "^>"g)} f' "${ref_ns5b}"  > subtype_ns5b_ref.fa
    cat subtype_ns5b_ref.fa ${sample_id}_ns5b_consensus.fa > mafftinput_ns5b_subtype.fa
    if [ \$ns5bgeno -ne "7" ];then
      grep -A 1 ">7_KU861171" ${ref_ns5b} >> mafftinput_ns5b_subtype.fa
    fi
    mafft --reorder --adjustdirection --anysymbol --thread 4 --auto mafftinput_ns5b_subtype.fa > ${sample_id}_mafftoutput_ns5b_subtype
    raxmlHPC -f d -p 12345 -# 10 -s ${sample_id}_mafftoutput_ns5b_subtype -m GTRCAT -n ${sample_id}_ns5b_subtype -o 7_KU861171

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
    def tree_type = tree_file.name.contains('_core_subtype') ? 'core_subtype' :
                    tree_file.name.contains('_ns5b_subtype') ? 'ns5b_subtype' :
                    tree_file.name.endsWith('_core') ? 'core' :
                    tree_file.name.endsWith('_ns5b') ? 'ns5b' : 'unknown'


    """
    Rscript ${projectDir}/bin/plot_tree.R --sample_id ${sample_id} --tree_file ${tree_file} --tree_type ${tree_type}
    
    """
}
