### Introduction

**nfc/rnaseq** is a bioinformatics analysis pipeline used for RNA sequencing data.

The workflow processes raw data from
 FastQ inputs ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/),
 [Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)),
  aligns the reads
   ([STAR](https://github.com/alexdobin/STAR) or
    [HiSAT2](https://ccb.jhu.edu/software/hisat2/index.shtml)),
     generates counts relative to genes
      ([featureCounts](http://bioinf.wehi.edu.au/featureCounts/),
       [StringTie](https://ccb.jhu.edu/software/stringtie/)) or transcripts
        ([Salmon](https://combine-lab.github.io/salmon/),
         [tximport](https://bioconductor.org/packages/release/bioc/html/tximport.html)) and performs extensive quality-control on the results
          ([RSeQC](http://rseqc.sourceforge.net/),
           [Qualimap](http://qualimap.bioinfo.cipf.es/),
            [dupRadar](https://bioconductor.org/packages/release/bioc/html/dupRadar.html),
             [Preseq](http://smithlabresearch.org/software/preseq/),
              [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html),
               [MultiQC](http://multiqc.info/)). See the [output documentation](docs/output.md) for more details of the results.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Quick Start

i. Install [`nextflow`](https://nf-co.re/usage/installation)

ii. Install one of [`docker`](https://docs.docker.com/engine/installation/), [`singularity`](https://www.sylabs.io/guides/3.0/user-guide/) or [`conda`](https://conda.io/miniconda.html)

iii. Download the pipeline and test it on a minimal dataset with a single command

```bash
nextflow run nfc/rnaseq -profile test,<docker/singularity/conda>
```

iv. Start running your own analysis!

```bash
nextflow run nfc/rnaseq -profile <docker/singularity/conda> --input samples.tsv
```

See [usage docs](docs/usage.md) for all of the available options when running the pipeline.

### Documentation

The nfc/rnaseq pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](https://nf-co.re/usage/installation)
2. Pipeline configuration
    * [Local installation](https://nf-co.re/usage/local_installation)
    * [Adding your own system config](https://nf-co.re/usage/adding_own_config)
    * [Reference genomes](https://nf-co.re/usage/reference_genomes)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](https://nf-co.re/usage/troubleshooting)

### Credits

These scripts were originally written for use at the [National Genomics Infrastructure](https://ngisweden.scilifelab.se), part of [SciLifeLab](http://www.scilifelab.se/) in Stockholm, Sweden, by Phil Ewels ([@ewels](https://github.com/ewels)) and Rickard Hammarén ([@Hammarn](https://github.com/Hammarn)).

## Citation

