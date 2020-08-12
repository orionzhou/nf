from snakemake import shell
input, output, params, threads, w, config = snakemake.input, snakemake.output, snakemake.params, snakemake.threads, snakemake.wildcards, snakemake.config

genome = w.genome
params.hybrid = config['x'][genome]['hybrid']
opt = params.opt

shell("""
    rm -rf {output.fna}* {output.fai}*
    rm -rf {output.chrom_bed} {output.chrom_size} {output.gap}

    mkdir -p {params.wdir}/{params.odir}
    cd {params.wdir}/{params.odir}
    rm -rf raw.fna.* renamed* map* raw.sizes
""")

merge_tag = '--merge_short' if w.genome != 'Mt_R108' else ''
if params.hybrid:
    shell("""
        cat {input} > {params.wdir}/{params.odir}/renamed.fna
        cd {params.wdir}/{params.odir}
        fasta.py size renamed.fna > renamed.sizes
        touch mapf.chain mapb.chain
    """)
else:
    params.gap = int(config['x'][genome]['gap'])
    params.prefix = config['x'][genome]['prefix']
    shell("""
        cd {params.wdir}/{params.odir}
        ln -sf ../download/raw.fna raw.fna
        fasta.py size raw.fna > raw.sizes

        fasta.py rename raw.fna renamed.fna mapf.bed mapb.bed \
            --opt {params.opt} {merge_tag} \
            --gap {params.gap} --prefix_chr {params.prefix}

        fasta.py size renamed.fna > renamed.sizes
        chain.py fromBed mapf.bed raw.sizes renamed.sizes > mapf.chain
        chainSwap mapf.chain mapb.chain
    """)

shell("""
    cd {params.wdir}
    ln -sf {params.odir}/renamed.fna 10_genome.fna
    cd ..

    samtools faidx {output.fna}
    fasta.py size --bed {output.fna} > {output.chrom_bed}
    cut -f1,3 {output.chrom_bed} > {output.chrom_size}
    fasta.py gaps {output.fna} > {output.gap}
""")
