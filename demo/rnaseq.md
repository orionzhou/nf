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

    cd /home/springer/zhoux379/git
    git clone https://github.com/orionzhou/nf.git

(or if you cloned the repo a while ago) update the local repo (if there are new changes):

    cd /home/springer/zhoux379/git/nf
    git pull

[Create a conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file) in order to run the pipeline:

    conda env create -n rnaseq -f nf/configs/environments/rnaseq.yml
    conda env list
    # you should now see a new environment named "rnaseq"

Add these environmental variables (with necessary modification) to your `~/.bashrc` or `~/.bash_profile`

    export NXF_HOME=/home/springer/zhoux379/git/nf
    export NXF_CACHE=/scratch.global/zhoux379/nf
    export NXF_EXECUTOR=slurm
    export NXF_CONDA_CACHEDIR=/home/springer/zhoux379/miniconda3/envs
    export NXF_WORK=$NXF_CACHE/work
    export NXF_TEMP=$NXF_CACHE/tmp
    export NXF_SINGULARITY_CACHEDIR=$NXF_CACHE/singularity
    export NXF_OPTS='-Xms1g -Xmx10g'

Log out and log in again (or run `source ~/.bashrc`) to make these variables into effect

Copy the test pipeline (`nf/demo/rnaseq`) to wherever you'd like to take a try:

    cp -rf /home/springer/zhoux379/git/nf/demo/rnaseq /home/springer/zhoux379/rnaseq_test
    cd /home/springer/zhoux379/rnaseq_test
    ls
    # genomes.yml  nextflow.config  reads.tsv  reads.xlsx

Make necessary changes to `reads.tsv` and `nextflow.config`, as well as the genome index configuration file (`genomes.yml`), then we can run the test pipeline using existing genome database on MSI-mesabi queue:

    nextflow run $NXF_HOME/rnaseq -params-file genomes.yml -profile mangi
