process findamplicon {
  //errorStrategy 'ignore'
  //{contig} here is filtered_contig

  tag { sample_id }

  publishDir "${params.outdir}/${sample_id}", pattern: "*.fa", mode:'copy'

  input:
  tuple val(sample_id), path(contigs)

  output:
  tuple val(sample_id), path("${sample_id}_ref_seqs_for_mapping.fa"), emit: ref_seqs_mapping, optional: true
  tuple val(sample_id), path("${sample_id}_mafft_core.fa"), emit: maffcore, optional: true
  tuple val(sample_id), path("${sample_id}_mafft_ns5b.fa"), emit: maffns5b, optional: true
  //tuple val(sample_id), path("RAxML_bestTree.${sample_id}_core"), emit: core_besttree, optional: true
  //tuple val(sample_id), path("${sample_id}_core_consensus.fa"), emit: core_consensus, optional: true
  //tuple val(sample_id), path("${sample_id}_ns5b_consensus.fa"), emit: ns5b_consensus, optional: true

  script:
  """
 
  grep -A1 core ${contigs} > core_contig.fa

  if [ -s core_contig.fa ]; then
    grep '>' core_contig.fa > core_name
    awk -F'|' 'BEGIN{OFS="_"} {print \$4,\$3}' core_name > core_segment
    grep -A1 -f core_segment ${params.ref_core} > best_core_ref.fa
    cat best_core_ref.fa core_contig.fa > core_contig_ref.fa
    mafft --reorder --adjustdirection --anysymbol --thread 4 --auto core_contig_ref.fa > ${sample_id}_mafft_core.fa
    findamplicon.py -i ${sample_id}_mafft_core.fa -o ref_seqs_for_mapping_core.fa
  fi

  grep -A1 ns5b ${contigs} > ns5b_contig.fa

  if [ -s ns5b_contig.fa ]; then
    grep '>' ns5b_contig.fa | awk -F'|' 'BEGIN{OFS="_"} {print \$4,\$3}' > ns5b_segment
    grep -A1 -f ns5b_segment ${params.ref_ns5b} > best_ns5b_ref.fa
    cat best_ns5b_ref.fa ns5b_contig.fa > ns5b_contig_ref.fa
    mafft --reorder --adjustdirection --anysymbol --thread 4 --auto ns5b_contig_ref.fa > ${sample_id}_mafft_ns5b.fa
    findamplicon.py -i ${sample_id}_mafft_ns5b.fa -o ref_seqs_for_mapping_ns5b.fa
  fi

  if [ -f ref_seqs_for_mapping_core.fa ] && [ -f ref_seqs_for_mapping_ns5b.fa ]
  then
    cat ref_seqs_for_mapping_core.fa ref_seqs_for_mapping_ns5b.fa > ${sample_id}_ref_seqs_for_mapping.fa
  elif [ -f ref_seqs_for_mapping_core.fa ]
  then
    cat ref_seqs_for_mapping_core.fa > ${sample_id}_ref_seqs_for_mapping.fa
  else
    cat ref_seqs_for_mapping_ns5b.fa > ${sample_id}_ref_seqs_for_mapping.fa
  fi

  

  """
}