
process blastp {
  label 'mid_memory'
  tag "${params.genomes[qry].alias}-${params.genomes[tgt].alias}"

  input:
  tuple val(qry), val(tgt), path('q.faa'), path('t.faa'), path(blastp)

  output:
  tuple val(qry), val(tgt), path("blastp.tsv")

  script:
  """
  blastp -db $blastp/db -query q.faa -num_threads ${task.cpus} -outfmt 6 -out out.tsv
  mv out.tsv blastp.tsv
  """
}

process syntelog {
  label 'low_memory'
  tag "${params.genomes[qry].alias}-${params.genomes[tgt].alias}"
  publishDir "${params.outdir}/${qry}-${tgt}", mode:'copy', overwrite: true

  input:
  tuple val(qry), val(tgt), path('q.t.last'), path('q.pep'), path('q0.bed'), path('t.pep'), path('t0.bed')

  output:
  tuple val(qry), val(tgt), path("xref.pairs"), path('xref.q.tsv'), path('xref.t.tsv')

  script:
  cscore = 0.7
  quota = '1:1'
  dist = 20
  Nm = 10
  quota_str = quota.replaceAll(":","x")
  """
  sort -k1,1 -k2,2n -k3,3n t0.bed > t.bed
  sort -k1,1 -k2,2n -k3,3n q0.bed > q.bed
  python -m jcvi.compara.blastfilter q.t.last --cscore=${cscore} --no_strip_names

  python -m jcvi.compara.synteny scan q.t.last.filtered q.t.anchors --dist=${dist} --no_strip_names
  $baseDir/bin/wgc/quota.py q.t.anchors --quota=${quota} --screen --Nm=${Nm}
  python -m jcvi.compara.synteny liftover --no_strip_names q.t.last q.t.${quota_str}.anchors
  ln -sf q.t.${quota_str}.lifted.anchors 01.anchors

  anchor.py 2tsv 01.anchors > xref.pairs
  python -m jcvi.compara.synteny mcscan q.bed 01.anchors --iter=1 --Nm=$Nm -o 06.q.blocks
  python -m jcvi.compara.synteny mcscan t.bed 01.anchors --iter=1 --Nm=$Nm -o 06.t.blocks

  python -m jcvi.formats.blast cscore --no_strip_names q.t.last > q.t.rbh
  $baseDir/bin/wgc/reconstruct.py fillrbh 06.q.blocks q.t.rbh xref.q.tsv
  $baseDir/bin/wgc/reconstruct.py fillrbh 06.t.blocks q.t.rbh xref.t.tsv
  """
  //python -m jcvi.compara.quota q.t.anchors --quota=${quota} --screen --Nm=${Nm}
  //python -m jcvi.compara.synteny simple 01.anchors --rich --qbed=q.bed --sbed=t.bed
  //python -m jcvi.graphics.dotplot 01.anchors --nostdpf --notex --qbed=q.bed --sbed=t.bed
}

workflow syn {
  take:
    comps
    qry_pfaa
    qry_pbed
    tgt_pfaa
    tgt_pbed
    tgt_blastp
  main:
    // set up channels
      tgt_pfaa1 = tgt_pfaa.map {r -> [r[1], r[0]]}
      tgt_pbed1 = tgt_pbed.map {r -> [r[1], r[0]]}
      tgt_blastp1 = tgt_blastp.map {r -> [r[1], r[0]]}
    in1 = comps.combine(qry_pfaa,by:0).combine(tgt_pfaa1,by:1)
      .combine(tgt_blastp,by:0)
      .map {r -> [r[1],r[0],r[2],r[3],r[4]]}
    blastp(in1)
    in2 = blastp.out.combine(qry_pfaa,by:0).combine(qry_pbed,by:0)
      .combine(tgt_pfaa1,by:1).combine(tgt_pbed,by:0)
      .map {r -> [r[1],r[0],r[2],r[3],r[4],r[5],r[6]]}
    syntelog(in2)
}
