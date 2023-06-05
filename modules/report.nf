process report {

    publishDir "${params.outdir}", pattern: "*run_summary_report.csv", mode:'copy'

    input:
    tuple path(fastqlist), path(combined_consensus),path(combined_demix), path(combined_qc), path(mapped_reads_counts), path(genotype_nt)

    output:
    path("*run_summary_report.csv")

    """
    report.py --fastqlist ${fastqlist} --consensus_report ${combined_consensus} \
    --demix_report ${combined_demix} --qc_report ${combined_qc} --reads_count ${mapped_reads_counts} --genotype_nt ${genotype_nt} --prefix ${params.prefix}

    """

}
