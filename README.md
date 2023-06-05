# hcv_nf

The process is partially adapted from [FluViewer tool](https://github.com/KevinKuchinski/FluViewer) for genotype HCV amplicon sequencing data. Samples are only amplied in the core (361-764) and ns5b region (8803-9191). This is a purely assembly-based approach. The assembly is done using SPades. Top 10 genotypes are produced when blast (blastn) the contigs to the database. 

The workflow is captured in the diagram below![diagram](pics/workflow.PNG).

## Usage

When using nextflow pipeline, specify the environment by adding ```-profile conda --cache ~/.conda/envs```

```
nextflow run BCCDC-PHL/hcv_nf \
  --fastq_input <path/to/fastq/dirs> \
  --db <path/to/ref/db> \
  --ref_core <path/to/ref_core/db> \
  --ref_ns5b <path/to/ref_ns5b/db> \
  --outdir <path/to/output_dir> \ 
```
### Input

The required inputs are:
- fastq input directory. 
- path to the full length HCV reference database
- path to reference database that have core side extraced
- path to reference database that have ns5b side extraced
- path to basic_qc_stats.csv which is the output of the pipeline BCCDC-PHL basic_qc
- path to abundance top 5 which is the output of the routine sequence qc pipeline maintained by BCCDC-PHL
- outdir directory to store the results


### Output

| outputs  | description |
| ------------- | ------------- |
| run_summary_report.csv  | the combined summary for consensus report, genotype, qc stats, demixming results  and check column |
| consensus_seqs.fa | consensus sequences for core and/or ns5b  |
| genotype_calls.csv | blastn results after blast the consensus sequences to the nt database, some columns are in the run_summary_report.csv|
| demix.csv | proportions of different subtypes present in the sample, are also in the run_summary_report.csv |
| parsed_genome_results.csv | qc stats for mean coverage, total mapped reads, median coverage, depth, percent completeness at different depth. also in the run_summary_report|
| mapped_to_db.bam | mapping raw reads to all references in the database|
| mapped_to_ref.bam | mapping raw reads to the assembly | 
| RAxML_bestTree.1Ao4_core | Tree with sample of interests and the core references|
| RAxML_bestTree.1Ao4_ns5b | Tree with sample of interests and the ns5b references|