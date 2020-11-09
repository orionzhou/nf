process prep_qry {
  label 'low_memory'
  tag "${params.genomes[qry].alias}"

  input:
  tuple val(qry), path(fas), path(gbed), path(sizes), path(gap)

  output:
  tuple val(qry), path("${qry}.chain"), emit: chain
  tuple val(qry), path("${qry}_split/p*.fas"), emit: seq

  script:
  """
  bed.py filter -min 5000 ${gap} > 01.qry.gap.bed
  subtractBed -nonamecheck -sorted -a ${gbed} \
      -b ${gap} | bed.py filter -min 100 - | \
      bed.py makewindow -w 1000000 -s 995000 - \
      > 02.qry.clean.bed
  bed.py size 02.qry.clean.bed
  bed.py binpacking 02.qry.clean.bed 03.bed ${qry}_split --N ${params.npieces} --pre p
  cut -f1,3 03.bed > 03.sizes
  chain.py fromBed 03.bed 03.sizes ${sizes} > ${qry}.chain
  ls ${qry}_split/*.bed | \\
    parallel -j ${task.cpus} "fasta.py extract ${fas} {} > {.}.fas"
  """
  //rm ${qry}_split/*.bed
}

process aln {
  label 'high_memory'
  tag "${params.genomes[qry].alias}-${params.genomes[tgt].alias}-${qfas.baseName}"

  input:
  tuple val(qry), val(tgt), path(qfas), path(tfas)

  output:
  tuple val(qry), val(tgt), path("*.psl")

  script:
  """
  minimap2 -x asm20 -a -t ${task.cpus} \\
    ${tfas} ${qfas} | \\
    sam.py 2psl - ${qfas.baseName}.psl
  """
}

process merge {
  label 'mid_memory'
  tag "${params.genomes[qry].alias}-${params.genomes[tgt].alias}"

  input:
  tuple val(qry), val(tgt), path(psls), path('q.sizes'), path('t.sizes')

  output:
  tuple val(qry), val(tgt), path("03.coord.pass.psl")

  script:
  """
  pslCat -nohead ${psls} > 01.psl
  psl.py coordQ 01.psl q.sizes > 02.coord.psl
  pslCheck 02.coord.psl \\
    -pass=02.coord.pass.psl \\
    -fail=02.coord.fail.psl || echo non_success
  psl.py coordT 02.coord.psl t.sizes > 03.coord.psl
  pslCheck 03.coord.psl \\
    -pass=03.coord.pass.psl \\
    -fail=03.coord.fail.psl || echo non_success
  """
}

process chain {
  label 'mid_memory'
  tag "${params.genomes[qry].alias}-${params.genomes[tgt].alias}"

  input:
  tuple val(qry), val(tgt), path(psl), path('q.sizes'), path('q.2bit'), path('t.sizes'), path('t.2bit')

  output:
  tuple val(qry), val(tgt), path("15.chain")

  script:
  """
  axtChain -linearGap=medium -psl $psl t.2bit q.2bit 10.chain
  chainPreNet 10.chain t.sizes q.sizes 11.chain
  chainSwap 11.chain 11.q.chain
  chainNet 11.chain t.sizes q.sizes 13.t.net 13.q.net
  netChainSubset 13.t.net 11.chain stdout |\\
          chainSort stdin 13.t.chain
  netChainSubset 13.q.net 11.q.chain stdout |\\
          chainSort stdin 13.q.chain
  chainNet 13.q.chain q.sizes t.sizes /dev/null 15.net
  netChainSubset 15.net 11.chain 15.chain
  """
}

process post {
  label 'mid_memory'
  tag "${params.genomes[qry].alias}-${params.genomes[tgt].alias}"
  publishDir "${params.outdir}/${qry}-${tgt}", mode:'copy', overwrite: true

  input:
  tuple val(qry), val(tgt), path(fi), path('q.fas'), path('q.sizes'), path('q.ref'), path('t.fas'), path('t.sizes'), path('t.ref')

  output:
  tuple val(qry), val(tgt), path("aln.bed"), path("vnt.bed"), path("q.chain"), path("t.chain"), path("t.vcf.gz"), path("q.vcf.gz")

  script:
  """
  chainStitchId $fi 01.stitched.chain
  chainFilter -minGapless=1000 01.stitched.chain > 02.filtered.chain
  chain.py 2bed 02.filtered.chain > 02.bed
  $baseDir/bin/wgc/chainBedVnt.R 02.bed 05.itv.bed
  $baseDir/bin/wgc/wgc.py callvnt 05.itv.bed t.fas q.fas --vnt 05.vnt.bed > 05.bed
  $baseDir/bin/wgc/chainBedFilter.R 05.bed aln.bed 05.vnt.bed vnt.bed
  chain.py fromBed aln.bed t.sizes q.sizes > t.chain
  chainSwap t.chain q.chain
  $baseDir/bin/wgc/wgc.py bed2vcf --tgt $tgt --qry $qry vnt.bed t.1.vcf q.1.vcf
  sortBed -header -i t.1.vcf > t.1.s.vcf
  sortBed -header -i q.1.vcf > q.1.s.vcf
  gatk UpdateVCFSequenceDictionary -V t.1.s.vcf -R t.ref/db.fasta -O t.2.vcf
  gatk UpdateVCFSequenceDictionary -V q.1.s.vcf -R q.ref/db.fasta -O q.2.vcf
  bcftools norm -f t.fas -c w -d all t.2.vcf -Ou | bcftools sort -Oz -o t.vcf.gz
  bcftools norm -f q.fas -c w -d all q.2.vcf -Ou | bcftools sort -Oz -o q.vcf.gz
  """
}

workflow wgc {
  take:
    comps
    qry_fas
    qry_gbed
    qry_sizes
    qry_gap
    qry_2bit
    qry_gff
    qry_pgff
    qry_gatk
    tgt_fas
    tgt_gbed
    tgt_sizes
    tgt_gap
    tgt_2bit
    tgt_gff
    tgt_pgff
    tgt_gatk
  main:
    // set up channels
      qrys0 = comps.map {r -> [r[0]]}
      tgts0 = comps.map {r -> [r[1]]}
      tgt_fas1 = tgt_fas.map {r -> [r[1], r[0]]}
      tgt_gbed1 = tgt_gbed.map {r -> [r[1], r[0]]}
      tgt_sizes1 = tgt_sizes.map {r -> [r[1], r[0]]}
      tgt_gap1 = tgt_gap.map {r -> [r[1], r[0]]}
      tgt_gatk1 = tgt_gatk.map {r -> [r[1], r[0]]}
      tgt_2bit1 = tgt_2bit.map {r -> [r[1], r[0]]}
      tgt_gff1 = tgt_gff.map {r -> [r[1], r[0]]}
      tgt_pgff1 = tgt_pgff.map {r -> [r[1], r[0]]}
    prep_qry(qry_fas.join(qry_gbed).join(qry_sizes).join(qry_gap))
    qchain = prep_qry.out.chain
    qseqs = prep_qry.out.seq
    in1 = comps.combine(qseqs,by:0).combine(tgt_fas1,by:1).transpose()
      .map {r -> [r[1],r[0],r[2],r[3]]}
    aln(in1)
    in2 = aln.out.groupTuple(by:[0,1])
      .combine(qry_sizes,by:0).combine(tgt_sizes1,by:1)
      .map {r -> [r[1],r[0],r[2],r[3],r[4]]}
    merge(in2)
    in3 = merge.out
      .combine(qry_sizes,by:0).combine(qry_2bit,by:0)
      .combine(tgt_sizes1,by:1).combine(tgt_2bit,by:0)
      .map {r -> [r[1],r[0],r[2],r[3],r[4],r[5],r[6]]}
    chain(in3)
    in4 = chain.out
      .combine(qry_fas,by:0).combine(qry_sizes,by:0).combine(qry_gatk,by:0)
      .combine(tgt_fas1,by:1).combine(tgt_sizes,by:0).combine(tgt_gatk,by:0)
      .map {r -> [r[1],r[0],r[2],r[3],r[4],r[5],r[6],r[7],r[8]]}
    post(in4)

//  emit:
}
