process report {

    publishDir "${params.outdir}", pattern: "run_summary_report.csv", mode:'copy'

    input:
    tuple path(fastqlist), path(combined_genotype), path(combined_consensus),path(combined_demix), path(combined_qc),path(basic_qc),path(abundance_top_n)

    output:
    path("run_summary_report.csv")

    """
    report.py --fastqlist ${fastqlist} --genotype_calls ${combined_genotype} --consensus_report ${combined_consensus} \
    --demix_report ${combined_demix} --qc_report ${combined_qc} --basic_qc ${basic_qc} --abundance_top_n ${abundance_top_n} --prefix ${params.prefix}

    """

}
