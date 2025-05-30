process findamplicon {

  tag { sample_id }

  publishDir "${params.outdir}/${sample_id}", pattern: "*.fa", mode:'copy'

  input:
  tuple val(sample_id), path(contigs), path(ref_core), path(ref_ns5b)

  output:
  tuple val(sample_id), path("${sample_id}_ref_seqs_for_mapping.fa"), emit: ref_seqs_mapping, optional: true
  tuple val(sample_id), path("${sample_id}_mafft_core.fa"), emit: maffcore, optional: true
  tuple val(sample_id), path("${sample_id}_mafft_ns5b.fa"), emit: maffns5b, optional: true

  script:
  """
 
  grep -A1 '|core|' ${contigs} > core_contig.fa

  if [ -s core_contig.fa ]; then
    grep '>' core_contig.fa > core_name
    awk -F'|' 'BEGIN{OFS="_"} {print \$4,\$3}' core_name > core_segment
    grep -A1 -f core_segment ${ref_core} > best_core_ref.fa
    findamplicon.py -i core_contig.fa -r best_core_ref.fa -o ref_seqs_for_mapping_core.fa -s ${sample_id}
  fi

  grep -A1 '|ns5b|' ${contigs} > ns5b_contig.fa

  if [ -s ns5b_contig.fa ]; then
    grep '>' ns5b_contig.fa | awk -F'|' 'BEGIN{OFS="_"} {print \$4,\$3}' > ns5b_segment
    grep -A1 -f ns5b_segment ${ref_ns5b} > best_ns5b_ref.fa
    findamplicon.py -i ns5b_contig.fa -r best_ns5b_ref.fa -o ref_seqs_for_mapping_ns5b.fa -s ${sample_id}
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