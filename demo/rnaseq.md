## Nextflow RNA-Seq pipeline

Install Conda

    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh

Intall Nextflow (in a new environment)

    conda create -n nf python=3
    conda activate nf
    conda install nextflow

Clone the repo

    git clone https://github.com/orionzhou/nf.git

Run a test pipeline using existing genome database on MSI-mesabi queue

    cd nf/demo/rnaseq
    nextflow run ../../rnaseq -params-file /home/springer/zhoux379/projects/genome/nf/genomes.yml -profile mesabi
