
process genotype {
    errorStrategy 'ignore'

    memory { 55.GB }
    tag {sample_id}

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}/${sample_id}*", mode:'copy', saveAs: { filename -> filename.split("/").last() }
    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}/logs", mode:'copy', saveAs: { filename -> "f-f logs" }


    input:
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    tuple val(sample_id), path("${sample_id}/${sample_id}_blast_results_prefilter.csv"), emit: allblastresult, optional: true
    tuple val(sample_id), path("${sample_id}/${sample_id}_filtered_blast_results.csv"), emit: blastreport, optional: true
    tuple val(sample_id), path("${sample_id}/${sample_id}_contigs.fa"), emit: contigs, optional: true
    tuple val(sample_id), path("${sample_id}/${sample_id}_filtered_contigs.fa"), emit: filtered_contigs, optional: true
    tuple val(sample_id), path("${sample_id}_genotype_provenance.yml"), emit: provenance
    tuple val(sample_id), path("${sample_id}/logs"), emit: fluviewer_logs

    """
    printf -- "- process_name: genotype\\n" > ${sample_id}_genotype_provenance.yml
    printf -- "  tool_name: Spades\\n  tool_version: \$(spades.py --version | cut -d' ' -f4)\\n" >> ${sample_id}_genotype_provenance.yml
    printf -- "  tool_name: blastn\\n  tool_version: \$(blastn -version | cut -d' ' -f2 | head -n 1)\\n" >> ${sample_id}_genotype_provenance.yml
    printf -- "  database used: ${params.db}\\n" >> ${sample_id}_genotype_provenance.yml
    printf -- "  database_path: \$(readlink -f ${params.db})\\n" >> ${sample_id}_genotype_provenance.yml
    printf -- "  database sha256: \$(shasum -a 256 ${params.db} | awk '{print \$1}')\\n" >> ${sample_id}_genotype_provenance.yml


    genotype_v1.py -f ${reads_1} -r ${reads_2} -o ${sample_id} -d ${params.db} -m ${params.mode}

    """

}

process makeconsensus {

    memory { 55.GB }
    tag {sample_id}

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}/${sample_id}*", mode:'copy', saveAs: { filename -> filename.split("/").last() }
    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}/logs", mode:'copy', saveAs: { filename -> "f-f logs" }


    input:
    tuple val(sample_id), path(reads_1), path(reads_2), path(ref_seq_map)

    output:
    tuple val(sample_id), path("${sample_id}/${sample_id}*.bam*"), emit: alignment, optional: true
    tuple val(sample_id), path("${sample_id}/${sample_id}_consensus_seqs_report.tsv"), emit: consensus_seqs_report, optional: true
    tuple val(sample_id), path("${sample_id}/${sample_id}_consensus_seqs.fa"), emit: consensus_seqs, optional: true
    tuple val(sample_id), path("${sample_id}/logs"), emit: fluviewer_logs

    """
    genotype_v1_2.py -f ${reads_1} -r ${reads_2} -s ${ref_seq_map} -o ${sample_id} -d ${params.db} -m ${params.mode}

    """

}


process blastconcensus {
    errorStrategy 'ignore'
    tag {sample_id}

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}*", mode:'copy'
    
    input:
    tuple val(sample_id), path(query), path(db_dir), path(db_name)

    output:
    tuple val(sample_id), path("${sample_id}_consensus_blast.csv"), emit: consensus_blast, optional:true
    tuple val(sample_id), path("${sample_id}_genotype_calls_nt.csv"), emit: genotyperesult, optional:true
    tuple val(sample_id), path("${sample_id}_seq_description.csv"), emit: seq_description, optional:true
    tuple val(sample_id), path("${sample_id}_blastconsensus_provenance.yml"), emit: provenance


    """
    printf -- "- process_name: blastconsensus\\n" > ${sample_id}_blastconsensus_provenance.yml
    printf -- "  tool_name: blastn\\n  tool_version: \$(blastn -version 2>&1 | cut -d' ' -f2 | head -n 1)\\n" >> ${sample_id}_blastconsensus_provenance.yml
    printf -- "  database used: ${db_dir}\\n" >> ${sample_id}_blastconsensus_provenance.yml
    printf -- "  database_path: \$(readlink -f ${db_dir})\\n" >> ${sample_id}_blastconsensus_provenance.yml

    export BLASTDB="${db_dir}"

    echo "query_seq_id,subject_accession,subject_strand,query_length,query_start,query_end,subject_length,subject_start,subject_end,alignment_length,percent_identity,percent_coverage,num_mismatch,num_gaps,e_value,bitscore,subject_taxids" > ${sample_id}_consensus_blast.csv

    blastn \
      -db ${db_name} \
      -num_threads ${task.cpus} \
      -perc_identity ${params.minid} \
      -qcov_hsp_perc ${params.mincov} \
      -query ${query} \
      -outfmt "6 qseqid saccver sstrand qlen qstart qend slen sstart send length pident qcovhsp mismatch gaps evalue bitscore staxids" \
    | tr \$"\\t" "," >> ${sample_id}_consensus_blast.csv

    tail -qn+2 ${sample_id}_consensus_blast.csv | cut -d',' -f2 | sort -u > seqids
    blastdbcmd -db ${db_name} -entry_batch seqids | grep '>' > ${sample_id}_seq_description

    tail -qn+2 ${sample_id}_consensus_blast.csv | cut -d',' -f17 | sed 's/;/\\n/g' | sort -u > taxids

    taxonkit lineage -r -n  taxids > ${sample_id}_taxon_results.txt

    bind_taxon.py -f ${sample_id}_taxon_results.txt -d ${sample_id}_seq_description -b ${sample_id}_consensus_blast.csv -o ${sample_id}_genotype_calls_nt.csv

    """

}
