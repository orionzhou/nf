## Nextflow RNA-Seq pipeline

Install Conda

    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh

Intall Nextflow (in a new environment)

    conda create -n nf python=3
    conda activate nf
    conda install nextflow
    conda list

Clone the repo to local:

    git clone https://github.com/orionzhou/nf.git

[Create a conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file) in order to run the pipeline:

    conda env create -f nf/configs/environments/rnaseq.yml
    conda env list

Run a test pipeline using existing genome database on MSI-mesabi queue

    cd nf/demo/rnaseq
    nextflow run ../../rnaseq -params-file /home/springer/zhoux379/projects/genome/nf/genomes.yml -profile mesabi
