## Nextflow RNA-Seq pipeline

Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html)

    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh

Intall [Nextflow](https://github.com/nextflow-io/nextflow) (in a new conda environment)

    conda create -n nf python=3
    conda activate nf
    conda install nextflow
    conda list

Clone the repo to local:

    git clone https://github.com/orionzhou/nf.git

(or if you cloned the repo a while ago) update the local repo (if there are new changes):

    cd nf
    git pull

[Create a conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file) in order to run the pipeline:

    conda env create -f nf/configs/environments/rnaseq.yml
    conda env list

A test pipeline is under `nf/demo/rnaseq`:

    cd nf/demo/rnaseq

Make necessary changes to `reads.tsv` and `nextflow.config`, as well as the genome index configuration file (`genomes.yml`), then we can run the test pipeline using existing genome database on MSI-mesabi queue:

    nextflow run ../../rnaseq -params-file genomes.yml -profile mesabi
