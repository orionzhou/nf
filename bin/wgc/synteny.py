from snakemake import shell
input, output, params, threads, w, config = snakemake.input, snakemake.output, snakemake.params, snakemake.threads, snakemake.wildcards, snakemake.config

shell("""
    mkdir -p {params.odir}
    cd {params.odir}
    ln -sf {input.tgt_faa} t.pep
    ln -sf {input.qry_faa} q.pep
    ln -sf {input.tgt_bed} t.bed
    ln -sf {input.qry_bed} q.bed
    ln -sf {params.in_last} q.t.last
    python -m jcvi.compara.blastfilter q.t.last \
        --cscore={params.cscore} --no_strip_names
""")

quota, Nm = params.quota, params.Nm
if w.qry == 'Sbicolor' and w.tgt == 'B73':
    quota = '1:2'
    Nm = 30
    quota_str = quota.replace(":","x")
    shell("""
        cd {params.odir}
        python -m jcvi.compara.synteny scan q.t.last.filtered q.t.anchors \
            --liftover=q.t.last \
            --dist={params.dist} --no_strip_names
        python -m jcvi.compara.quota q.t.lifted.anchors \
            --quota={quota} --screen --Nm={Nm}
        ln -sf q.t.lifted.{quota_str}.anchors 01.anchors

        reconstruct.py mergechrom 01.anchors --qbed=q.bed --sbed=t.bed
        anchor.py 2tsv 01.mergechrom.anchors > 05.pairs
        python -m jcvi.compara.synteny mcscan q.bed 01.mergechrom.anchors \
            --iter=2 --trackids --Nm={Nm} -o 06.q.blocks
    """)
else:
    if w.qry == 'Athaliana' and w.tgt == 'Zmays_B73':
        quota = '3:4'
    quota_str = quota.replace(":","x")
    shell("""
        cd {params.odir}
        python -m jcvi.compara.synteny scan q.t.last.filtered q.t.anchors \
            --dist={params.dist} --no_strip_names
        python -m jcvi.compara.quota q.t.anchors \
            --quota={quota} --screen --Nm={Nm}
        python -m jcvi.compara.synteny liftover \
            --no_strip_names q.t.last q.t.{quota_str}.anchors
        ln -sf q.t.{quota_str}.lifted.anchors 01.anchors

        anchor.py 2tsv 01.anchors > 05.pairs
        python -m jcvi.compara.synteny mcscan q.bed 01.anchors \
            --iter=1 --Nm={Nm} -o 06.q.blocks
        python -m jcvi.compara.synteny mcscan t.bed 01.anchors \
            --iter=1 --Nm={Nm} -o 06.t.blocks

        python -m jcvi.formats.blast cscore --no_strip_names q.t.last > q.t.rbh
        reconstruct.py fillrbh 06.q.blocks q.t.rbh 07.q.ortholog
        reconstruct.py fillrbh 06.t.blocks q.t.rbh 07.t.ortholog
    """)

shell("""
    cd {params.odir}
    python -m jcvi.compara.synteny simple 01.anchors --rich \
        --qbed=q.bed --sbed=t.bed
    python -m jcvi.graphics.dotplot 01.anchors --nostdpf --notex \
        --qbed=q.bed --sbed=t.bed

    cp 05.pairs {output.pair}
    cp 01.pdf {output.plot}
    cp 07.q.ortholog {output.q}
    cp 07.t.ortholog {output.t}
""")
#rm {params.tfaa} {params.qfaa} {params.last}
