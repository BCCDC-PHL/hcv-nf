# hcv_nf

The process for genotyping HCV is partially adapted from [FluViewer tool](https://github.com/KevinKuchinski/FluViewer).

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

The required inputs are fastq input directory, path to the full length HCV reference database, pathes to reference database that have core side extraced, and ns5b extracted.

`--mode` can be either align or assemble, default is assemble. It is recommended to use assemble mode for HCV genotyping.

