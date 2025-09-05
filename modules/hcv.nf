process assemble {

    tag { sample_id }
    
    input:
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    tuple val(sample_id), path("${sample_id}/${sample_id}_contigs.fa"), emit: contigs,optional: true
    tuple val(sample_id), path("${sample_id}_assemble_provenance.yml"), emit: provenance
    tuple val(sample_id), path("${sample_id}/logs"), emit: spades_logs

    script:
    /*
    This script does:
    1. Run initial SPAdes assembly (--rnaviral --isolate -k 15 21 25)
    2. Check if reads are sufficient for subsampling
    3. Subsample with reformat.sh if sufficient reads, else skip
    4. Run SPAdes careful assembly on subsampled or original data
    5. Concatenate contigs/scaffolds from initial and careful assemblies
    */

    def spades_out_1 = "${sample_id}/${sample_id}_spades_results"
    def spades_out_2 = "${sample_id}/${sample_id}_spades_results_careful"
    def subsampled_1 = "${sample_id}_sampled_R1.fq"
    def subsampled_2 = "${sample_id}_sampled_R2.fq"
    def contigs_final = "${sample_id}/${sample_id}_contigs.fa"

    // set a minimum number of reads threshold to run reformat.sh, e.g. 10000 reads per file
    def min_reads = 10000

    """
    printf -- "- process_name: assemble\\n" > ${sample_id}_assemble_provenance.yml
    printf -- "  tool_name: Spades\\n  tool_version: \$(spades.py --version | cut -d' ' -f4)\\n" >> ${sample_id}_assemble_provenance.yml


    mkdir -p ${sample_id} ${sample_id}/logs

    echo "Running initial SPAdes assembly..."
    spades.py --rnaviral --isolate -k 15,21,25 -1 ${reads_1} -2 ${reads_2} -o ${spades_out_1} \
      > ${sample_id}/logs/${sample_id}_spades_stdout.txt 2> ${sample_id}/logs/${sample_id}_spades_stderr.txt

    # Count the number of reads in file 1 (assumes FASTQ, 4 lines per read)
    reads_count=\$(expr \$(zcat ${reads_1} | wc -l) / 4)

    echo "Reads count in ${reads_1}: \$reads_count"

    if [ \$reads_count -ge ${min_reads} ]; then
        echo "Enough reads for subsampling, running reformat.sh to sample 10%..."
        reformat.sh in1=${reads_1} in2=${reads_2} out1=${subsampled_1} out2=${subsampled_2} samplerate=0.1 sampleseed=1234
        subsample_success=\$?
    else
        echo "Not enough reads for subsampling, skipping reformat.sh."
        subsample_success=1
    fi

    if [ \$subsample_success -ne 0 ]; then
        echo "Running SPAdes careful assembly without subsampled reads..."
        spades.py -k 15,21,25 --careful --only-assembler -1 ${reads_1} -2 ${reads_2} -o ${spades_out_2} || true
    else
        echo "Running SPAdes careful assembly with subsampled reads..."
        spades.py -k 15,21,25 --careful --only-assembler -1 ${subsampled_1} -2 ${subsampled_2} -o ${spades_out_2} \
          > ${sample_id}/logs/${sample_id}_spades_careful_stdout.txt 2> ${sample_id}/logs/${sample_id}_spades_careful_stderr.txt || true
    fi

    # Determine contigs/scaffolds files
    
    primary_contigs=${spades_out_1}/scaffolds.fasta
    careful_contigs=${spades_out_2}/scaffolds.fasta

    # Fallback to contigs.fasta if scaffolds.fasta not found
    if [ ! -f \$primary_contigs ]; then
        echo "Primary scaffolds not found, fallback to contigs.fasta"
        primary_contigs=${spades_out_1}/contigs.fasta
    fi

    if [ ! -f \$careful_contigs ]; then
        echo "Careful scaffolds not found, fallback to contigs.fasta"
        careful_contigs=${spades_out_2}/contigs.fasta
    fi

    # Concatenate contigs
    if [ -f \$primary_contigs ] && [ -f \$careful_contigs ]; then
        cat \$primary_contigs \$careful_contigs > ${contigs_final}
    elif [ -f \$primary_contigs ]; then
        cat \$primary_contigs > ${contigs_final}
    elif [ -f \$careful_contigs ]; then
        cat \$careful_contigs > ${contigs_final}
    else
        echo "None of the input files exist."
    fi
    echo "Assembly completed. Output: ${contigs_final}"
    """
}

process assemble_shovill {
    errorStrategy 'ignore'
    tag { sample_id }
    
    input:
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    tuple val(sample_id), path("${sample_id}/contigs.fa"), emit: contigs,optional: true
    tuple val(sample_id), path("${sample_id}_assemble_provenance.yml"), emit: provenance

    script:


    """   
    printf -- "- process_name: assemble\\n" > ${sample_id}_assemble_provenance.yml
    printf -- "  tool_name: shovill\\n  tool_version: \$(shovill --version | cut -d' ' -f2)\\n" >> ${sample_id}_assemble_provenance.yml

    shovill --outdir ${sample_id} --R1 $reads_1 --R2 $reads_2

    """

}


process blastn_and_filter {
    //errorStrategy 'ignore'
    memory { 55.GB }
    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}*", mode:'copy'

    input:
    tuple val(sample_id), path(contigs_fasta),path(ref_seqs_db)             


    output:
    tuple val(sample_id), path("${sample_id}_blast_results_prefilter.csv"), emit: allblastresult, optional: true
    tuple val(sample_id), path("${sample_id}_filtered_blast_results.csv"), emit: blastreport, optional: true
    tuple val(sample_id), path("${sample_id}_filtered_contigs.fa"), emit: filtered_contigs, optional: true
    tuple val(sample_id), path("${sample_id}_blastn_provenance.yml"), emit: provenance
    //tuple val(sample_id), path("${sample_id}/logs"), emit: fluviewer_logs
    script:


    """
    printf -- "- process_name: blastn_and_filter\\n" > ${sample_id}_blastn_provenance.yml
    printf -- "  tool_name: blastn\\n  tool_version: \$(blastn -version | cut -d' ' -f2 | head -n 1)\\n" >> ${sample_id}_blastn_provenance.yml
    printf -- "  database used: ${ref_seqs_db}\\n" >> ${sample_id}_blastn_provenance.yml
    printf -- "  database_path: \$(readlink -f ${ref_seqs_db})\\n" >> ${sample_id}_blastn_provenance.yml
    printf -- "  database sha256: \$(shasum -a 256 ${ref_seqs_db} | awk '{print \$1}')\\n" >> ${sample_id}_blastn_provenance.yml

    makeblastdb -in ${ref_seqs_db} -dbtype nucl


    # Run blastn alignment
    blastn -query ${contigs_fasta} -db ${ref_seqs_db} -outfmt "6 qseqid sseqid pident qlen slen mismatch gapopen qstart qend sstart send bitscore" > ${sample_id}_blast_results.tsv
    if [ -s ${sample_id}_blast_results.tsv ]; then
    filter_best_alignments.py -o ${sample_id} -f ${sample_id}_blast_results.tsv -d ${ref_seqs_db} -m ${params.mode} -t ${contigs_fasta}
    else
    echo "blast result is empty"
    fi
    
    """
}


process makeconsensus {

    memory { 55.GB }
    tag {sample_id}

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}/${sample_id}*", mode:'copy', saveAs: { filename -> filename.split("/").last() }

    input:
    tuple val(sample_id), path(reads_1), path(reads_2), path(ref_seq_map), path(db)

    output:
    tuple val(sample_id), path("${sample_id}/${sample_id}*.bam*"), emit: alignment, optional: true
    tuple val(sample_id), path("${sample_id}/${sample_id}_consensus_seqs_report.tsv"), emit: consensus_seqs_report, optional: true
    tuple val(sample_id), path("${sample_id}/${sample_id}_consensus_seqs.fa"), emit: consensus_seqs, optional: true
    tuple val(sample_id), path("${sample_id}/logs"), emit: fluviewer_logs

    """
    makeconsensus.py -f ${reads_1} -r ${reads_2} -s ${ref_seq_map} -o ${sample_id} -d ${db} -m ${params.mode}

    """

}


process blastconsensus {

    errorStrategy 'ignore'

    memory '256GB'
    cpus 16

    tag {sample_id}

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_genotype_calls_nt.csv", mode:'copy'
    
    input:
    tuple val(sample_id), path(query), path(db_dir), val(db_name)

    output:
    tuple val(sample_id), path("${sample_id}_consensus_blast.csv"), emit: consensus_blast, optional:true
    tuple val(sample_id), path("${sample_id}_genotype_calls_nt.csv"), emit: genotyperesult, optional:true
    tuple val(sample_id), path("${sample_id}_seq_description"), emit: seq_description, optional:true
    tuple val(sample_id), path("${sample_id}_blastconsensus_provenance.yml"), emit: provenance


    """
    printf -- "- process_name: blastconsensus\\n" > ${sample_id}_blastconsensus_provenance.yml
    printf -- "  tool_name: blastn\\n  tool_version: \$(blastn -version 2>&1 | cut -d' ' -f2 | head -n 1)\\n" >> ${sample_id}_blastconsensus_provenance.yml
    printf -- "  database used: ${db_dir}\\n" >> ${sample_id}_blastconsensus_provenance.yml
    printf -- "  database_path: \$(readlink -f ${db_dir})\\n" >> ${sample_id}_blastconsensus_provenance.yml

    export BLASTDB="${db_dir}"

    echo "query_seq_id,subject_accession,subject_strand,query_length,query_start,query_end,subject_length,subject_start,subject_end,alignment_length,percent_identity,percent_coverage,num_mismatch,num_gaps,e_value,bitscore,subject_taxids" > ${sample_id}_consensus_blast.csv

    blastn \
      -db ${db_dir}/${db_name} \
      -num_threads ${task.cpus} \
      -perc_identity ${params.minid} \
      -qcov_hsp_perc ${params.mincov} \
      -query ${query} \
      -outfmt "6 qseqid saccver sstrand qlen qstart qend slen sstart send length pident qcovhsp mismatch gaps evalue bitscore staxids" \
    | tr \$"\\t" "," >> ${sample_id}_consensus_blast.csv

    tail -qn+2 ${sample_id}_consensus_blast.csv | cut -d',' -f2 | sort -u > seqids
    blastdbcmd -db ${db_dir}/${db_name} -entry_batch seqids | grep '>' > ${sample_id}_seq_description

    tail -qn+2 ${sample_id}_consensus_blast.csv | cut -d',' -f17 | sed 's/;/\\n/g' | sort -u > taxids

    taxonkit lineage -r -n  taxids > ${sample_id}_taxon_results.txt

    bind_taxon.py -f ${sample_id}_taxon_results.txt -d ${sample_id}_seq_description -b ${sample_id}_consensus_blast.csv -o ${sample_id}_genotype_calls_nt.csv

    """

}
