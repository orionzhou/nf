# Nextflow RNA-Seq pipeline

* Install Conda

* Intall Nextflow

    conda install nextflow
    conda activate nextflow

* Clone the repo

    git clone https://github.com/orionzhou/nf.git

* Run a test pipeline using existing genome database on MSI-mesabi queue

    cd nf/demo/rnaseq
    nextflow run ../../rnaseq -params-file /home/springer/zhoux379/projects/genome/nf/genomes.yml -profile mesabi
