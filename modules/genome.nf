process unzip {
  label 'low_memory'
  tag "${params.genomes[id].alias}"

  input:
  tuple val(id), val(source), path("raw.fasta.gz"), path("raw.gff.gz"), val(url_fas), val(url_gff)

  output:
  tuple val(id), path("raw.fasta"), emit: seq
  tuple val(id), path("raw.gff"), emit: gff

  script:
  //species = species.replaceAll('\\s','_')
  //assembly = assembly.replaceAll('\\s','_')
  if( source == 'local' ) 
    """
    ln -sf ${url_fas} raw.fasta
    ln -sf ${url_gff} raw.gff
    """
  else
    """
    gunzip -c raw.fasta.gz > raw.fasta
    gunzip -c raw.gff.gz > raw.gff
    """
}

process seqfmt {
  label 'low_memory'
  tag "${params.genomes[id].alias}"
  publishDir "${params.outdir}/${id}", mode:"copy", overwrite:'true',
    saveAs: {fn ->
      if (fn.indexOf("10.fasta") >= 0) "$fn"
      else "15_intervals/$fn"
      }

  input:
  tuple val(id), path(fi), val(gap), val(chr_prefix)

  output:
  tuple val(id), path("10.fasta"), path("10.fasta.fai"), emit: seq
  tuple val(id), path("01.chrom.sizes"), emit: chrom_size
  tuple val(id), path("01.chrom.bed"), emit: chrom_bed
  tuple val(id), path("11.gap.bed"), emit: gap
  tuple val(id), path("forward.bed"), path("forward.chain"), emit: fchain
  tuple val(id), path("reverse.bed"), path("reverse.chain"), emit: rchain

  script:
  def merge_tag = id=='Mt_R108' ? '' : '--merge_short'
  """
  mv $fi raw0.fasta
  fasta.py clean raw0.fasta raw.fasta
  fasta.py size $fi raw.sizes

  fasta.py rename $fi renamed.fna forward.bed reverse.bed \\
    --opt ${id} ${merge_tag} \\
    --gap ${gap.replaceAll(/\.[0-9]+$/,'')} --prefix_chr ${chr_prefix}

  fasta.py size renamed.fna renamed.sizes
  chain.py fromBed forward.bed raw.sizes renamed.sizes forward.chain
  chainSwap forward.chain reverse.chain
  mv renamed.fna 10.fasta

  samtools faidx 10.fasta
  fasta.py size --bed 10.fasta 01.chrom.bed
  cut -f1,3 01.chrom.bed > 01.chrom.sizes
  fasta.py gaps 10.fasta 11.gap.bed
  """
}

process gff_cln {
  label 'low_memory'
  tag "${params.genomes[id].alias}"
  publishDir "${params.outdir}/$id/50_annotation", mode:"copy", overwrite:'true'

  input:
  tuple val(id), path(fi), val(fixopt), path(fbed), path(fchain)

  output:
  tuple val(id), path("10.gff")

  script:
  """
  gff.py fix --opt ${fixopt} $fi > 01.fixed.gff
  liftOver -gff 01.fixed.gff ${fchain} 10.gff unmapped
  """
}

process gff_idx {
  label 'mid_memory'
  tag "${params.genomes[id].alias}"
  publishDir "${params.outdir}/$id/50_annotation", mode:"copy", overwrite:'true'

  input:
  tuple val(id), path(ref), path(fai), path(gff), path(sizes)

  output:
  tuple val(id), path("10.gff.db"), emit: db
  tuple val(id), path("10.tsv"), emit: tsv
  tuple val(id), path("10.gtf"), emit: gtf
  tuple val(id), path("10.bed"), emit: bed
  tuple val(id), path("10.desc.tsv"), emit: desc
  tuple val(id), path("10.nt.fasta"), emit: fna
  tuple val(id), path("10.aa.fasta"), emit: faa
  tuple val(id), path("10.sqlite"), emit: txdb
  tuple val(id), path("15.gff"), emit: pgff
  tuple val(id), path("15.gff.db"), emit: pdb
  tuple val(id), path("15.tsv"), emit: ptsv
  tuple val(id), path("15.gtf"), emit: pgtf
  tuple val(id), path("15.bed"), emit: pbed
  tuple val(id), path("15.desc.tsv"), emit: pdesc
  tuple val(id), path("15.nt.fasta"), emit: pfna
  tuple val(id), path("15.aa.fasta"), emit: pfaa

  script:
  """
	gff.py index $gff 10.gff.db
	gff.py 2tsv $gff > 10.tsv
	gff.py 2gtf $gff > 10.gtf
	gff.py 2bed12 $gff > 10.bed
	gff.py note --attribute note1,note2 $gff > 10.desc.tsv
	python -m jcvi.formats.gff load $gff $ref -o 10.nt.fasta
	fasta.py translate 10.nt.fasta > 10.aa.fasta
	$baseDir/bin/genome/gff2txdb.R $gff $sizes 10.sqlite

	gff.py picklong $gff > 15.gff
	gff.py index 15.gff 15.gff.db
	gff.py 2tsv 15.gff > 15.tsv
	gff.py 2gtf 15.gff > 15.gtf
	gff.py 2bed12 15.gff > 15.bed
	gff.py note --attribute note1,note2 15.gff > 15.desc.tsv
	python -m jcvi.formats.gff load 15.gff $ref -o 15.nt.fasta
	fasta.py translate 15.nt.fasta > 15.aa.fasta
	$baseDir/bin/genome/gff2txdb.R 15.gff $sizes 15.sqlite
	"""
}

// sequence db
process i_blat {
  label 'low_memory'
  tag "${params.genomes[id].alias}"
  publishDir "${params.outdir}/$id/21_dbs/blat", mode:"copy", overwrite:'true'

  when: tag == 'T'

  input:
  tuple val(id), path(ref), path(fai), val(tag)

  output:
  tuple val(id), path("db.2bit"), path("db.2bit.tile11.ooc")

  script:
  """
  faToTwoBit $ref db.2bit
  blat db.2bit tmp.fas tmp.out -makeOoc=db.2bit.tile11.ooc
  """
}

process i_gatk {
  label 'low_memory'
  tag "${params.genomes[id].alias}"
  publishDir "${params.outdir}/$id/21_dbs", mode:"copy", overwrite:'true'

  when: tag == 'T'

  input:
  tuple val(id), path(ref), path(fai), val(tag)

  output:
  tuple val(id), val(pre), path("gatk/*")

  script:
  pre = "gatk/db"
  """
  mkdir gatk
  cp -f $ref gatk/db.fasta
  gatk CreateSequenceDictionary -R gatk/db.fasta
  samtools faidx gatk/db.fasta
  """
}

process i_bwa {
  label 'low_memory'
  tag "${params.genomes[id].alias}"
  publishDir "${params.outdir}/$id/21_dbs", mode:"copy", overwrite:'true'

  when: tag == 'T'

  input:
  tuple val(id), path(ref), path(fai), val(tag)

  output:
  tuple val(id), val(pre), path("bwa/*")

  script:
  pre = "bwa/db"
  """
  mkdir bwa
  ln -f $ref bwa/db.fasta
  bwa index -a bwtsw -p bwa/db $ref
  """
}

process i_bismark {
  label 'low_memory'
  tag "${params.genomes[id].alias}"
  publishDir "${params.outdir}/$id/21_dbs", mode:"copy", overwrite:'true'

  when: tag == 'T'

  input:
  tuple val(id), path(ref), path(fai), val(tag)

  output:
  tuple val(id), val(pre), path("bismark/*")

  script:
  pre = "bismark"
  """
  mkdir bismark
  ln -f $ref bismark/db.fasta
  bismark_genome_preparation --bowtie2 bismark/
  """
}

// fna + faa
process i_blastn {
  label 'low_memory'
  tag "${params.genomes[id].alias}"
  publishDir "${params.outdir}/$id/21_dbs", mode:"copy", overwrite:'true'

  when: tag == 'T'

  input:
  tuple val(id), path(fna), val(tag)

  output:
  tuple val(id), val(pre), path("blastn/*")

  script:
  pre = "blastn/db"
  """
  mkdir blastn
  makeblastdb -dbtype nucl -in $fna -title db -out blastn/db
  """
}

process i_blastp {
  label 'low_memory'
  tag "${params.genomes[id].alias}"
  publishDir "${params.outdir}/$id/21_dbs", mode:"copy", overwrite:'true'

  when: tag == 'T'

  input:
  tuple val(id), path(faa), val(tag)

  output:
  tuple val(id), val(pre), path("blastp/*")

  script:
  pre = "blastp/db"
  """
  mkdir blastp
  makeblastdb -dbtype prot -in $faa -title db -out blastp/db
  """
}

process i_lastn {
  label 'low_memory'
  tag "${params.genomes[id].alias}"
  publishDir "${params.outdir}/$id/21_dbs", mode:"copy", overwrite:'true'

  when: tag == 'T'

  input:
  tuple val(id), path(fna), val(tag)

  output:
  tuple val(id), val(pre), path("lastn/*")

  script:
  pre = "lastn/db"
  """
  mkdir lastn
  lastdb lastn/db $fna
  """
}

process i_lastp {
  label 'low_memory'
  tag "${params.genomes[id].alias}"
  publishDir "${params.outdir}/$id/21_dbs", mode:"copy", overwrite:'true'

  when: tag == 'T'

  input:
  tuple val(id), path(faa), val(tag)

  output:
  tuple val(id), val(pre), path("lastp/*")

  script:
  pre = "lastp/db"
  """
  mkdir lastp
  lastdb -p lastp/db $faa
  """
}

process i_salmon {
  label 'low_memory'
  tag "${params.genomes[id].alias}"
  publishDir "${params.outdir}/$id/21_dbs", mode:"copy", overwrite:'true'

  when: tag == 'T'

  input:
  tuple val(id), path(fna), val(tag)

  output:
  tuple val(id), val(pre), path("salmon/*")

  script:
  pre = "salmon/db"
  """
  mkdir salmon
  salmon index -p ${task.cpus} -t $fna --gencode -i salmon/db
  """
}

// fasta + gtf
process i_star {
  label 'low_memory'
  tag "${params.genomes[id].alias}"
  publishDir "${params.outdir}/$id/21_dbs", mode:"copy", overwrite:'true'

  when: tag == 'T'

  input:
  tuple val(id), path(ref), path(fai), path(gtf), val(tag)

  output:
  tuple val(id), val(pre), path("star/*")

  script:
  pre = "star"
  """
  mkdir star
  STAR --runThreadN ${task.cpus} --runMode genomeGenerate \\
    --genomeDir star/ \\
    --genomeFastaFiles $ref --sjdbGTFfile $gtf
  """
}

process i_hisat2 {
  label 'low_memory'
  tag "${params.genomes[id].alias}"
  publishDir "${params.outdir}/$id/21_dbs", mode:"copy", overwrite:'true'

  when: tag == 'T'

  input:
  tuple val(id), path(ref), path(fai), path(gtf), val(tag)

  output:
  tuple val(id), val(pre), path("hisat2/*")

  script:
  pre = "hisat2/db"
  """
  mkdir hisat2
  hisat2_extract_exons.py $gtf > hisat2/db.exon
  hisat2_extract_splice_sites.py $gtf > hisat2/db.ss
  hisat2-build -p ${task.cpus} --ss hisat2/db.ss \\
    --exon hisat2/db.exon $ref hisat2/db
  """
}

process i_snpeff {
  label 'low_memory'
  tag "${params.genomes[id].alias}"
  publishDir "${params.outdir}/$id/21_dbs", mode:"copy", overwrite:'true'

  when: tag == 'T'

  input:
  tuple val(id), path(ref), path(fai), path(gtf), val(tag)

  output:
  tuple val(id), path("snpeff/$cfg"), path("snpeff/*")

  script:
  cfg = 'snpEff.config'
  """
  mkdir -p snpeff/$id
  echo 'data.dir = .' > snpeff/$cfg
  echo '${id}.genome : Zea mays' >> snpeff/$cfg
  ln -f $ref snpeff/$id/sequences.fa
  ln -f $gtf snpeff/$id/genes.gtf
  snpEff -Xmx${task.memory.toGiga()}G build -c snpeff/$cfg -gtf22 -v $id
  """
}

process i_tandup {
  label 'mid_memory'
  tag "${params.genomes[id].alias}"
  publishDir "${params.outdir}/$id/52_tandup", mode:"copy", overwrite:'true'

  when: tag == 'T'

  input:
  tuple val(id), path(fna), path(faa), path(gbed), val(blastn_pre), path("blastn/*"), val(blastp_pre), path("blastp/*"), val(tag)

  output:
  tuple val(id), path("05.tandup.cds.tsv"), path("05.tandup.pro.tsv")

  script:
  """
  blastn -db ${blastn_pre} -query $fna -num_threads ${task.cpus} \
      -outfmt 6 -out blastn.tsv
  blastp -db ${blastp_pre} -query $faa -num_threads ${task.cpus} \
      -outfmt 6 -out blastp.tsv
  python -m jcvi.compara.catalog tandem --strip_gene_name=None blastn.tsv \
      $fna $gbed > 05.tandup.cds.tsv
  python -m jcvi.compara.catalog tandem --strip_gene_name=None blastp.tsv \
      $faa $gbed > 05.tandup.pro.tsv
  """
}

process i_rcfg {
  label 'mid_memory'
  tag "${params.genomes[id].alias}"
  publishDir "${params.outdir}/$id", mode:"copy", overwrite:'true'
  //conda "$NXF_CONDA_CACHEDIR/r"

  when: tag == 'T'

  input:
  tuple val(id), path(chrom_size), path(chrom_bed), path(gap_bed), path(ptsv), path(pdes), val(tag)

  output:
  tuple val(id), path("55.rds")

  script:
  """
  $baseDir/bin/genome/genome.prep.R --dirg ${params.outdir} $id 55.rds
  """
}

workflow genome {
  take:
    design
    pick
  main:
    // stage channels
      gtable = design
        .splitCsv(header:true, sep:'\t')
        .filter { it.genome in pick }
      gtable1 = gtable
        .map {r -> [r.genome, r.source,
          file("${params.outdir}/${r.genome}/raw/raw.fasta.gz", checkIfExists:false),
          file("${params.outdir}/${r.genome}/raw/raw.gff.gz", checkIfExists:false),
          r.url_fas, r.url_gff
          ]}
      gtable2 = gtable.map {r -> [r.genome, r.gap, r.prefix]}
      gtable3 = gtable.map {r -> [r.genome, r.gffopt]}
      gt_blat = gtable.map {r -> [r.genome, r.blat]}
      gt_bwa = gtable.map {r -> [r.genome, r.bwa]}
      gt_star = gtable.map {r -> [r.genome, r.star]}
      gt_gatk = gtable.map {r -> [r.genome, r.gatk]}
      gt_hisat2 = gtable.map {r -> [r.genome, r.hisat2]}
      gt_salmon = gtable.map {r -> [r.genome, r.salmon]}
      gt_snpeff = gtable.map {r -> [r.genome, r.snpeff]}
      gt_blast = gtable.map {r -> [r.genome, r.blast]}
      gt_last = gtable.map {r -> [r.genome, r.last]}
      gt_bismark = gtable.map {r -> [r.genome, r.bismark]}
      gt_tandup = gtable.map {r -> [r.genome, r.tandup]}
      gt_rcfg = gtable.map {r -> [r.genome, r.rcfg]}

    unzip(gtable1)
    seq = unzip.out.seq; gff = unzip.out.gff
    seqfmt(seq.combine(gtable2, by:0))
    seq = seqfmt.out.seq
    chrom_size = seqfmt.out.chrom_size; chrom_bed = seqfmt.out.chrom_bed
    gap = seqfmt.out.gap

    gff_cln(gff.combine(gtable3, by:0).combine(seqfmt.out.fchain, by:0))
    gff_idx(seq.combine(gff_cln.out, by:0).combine(chrom_size, by:0))
    gff = gff_cln.out; pgff = gff_idx.out.pgff; gtf = gff_idx.out.gtf
    pfna = gff_idx.out.pfna; pfaa = gff_idx.out.pfaa; pbed = gff_idx.out.pbed
    ptsv = gff_idx.out.ptsv; pdesc = gff_idx.out.pdesc

    seq.combine(gt_blat, by:0) | i_blat
    seq.combine(gt_gatk, by:0) | i_gatk
    seq.combine(gt_bwa, by:0) | i_bwa
    seq.combine(gt_bismark, by:0) | i_bismark
    pfna.combine(gt_blast, by:0) | i_blastn
    pfna.combine(gt_last, by:0) | i_lastn
    pfna.combine(gt_salmon, by:0) | i_salmon
    pfaa.combine(gt_blast, by:0) | i_blastp
    pfaa.combine(gt_last, by:0) | i_lastp
    seq.combine(gtf, by:0).combine(gt_star, by:0) | i_star
    seq.combine(gtf, by:0).combine(gt_hisat2, by:0) | i_hisat2
    seq.combine(gtf, by:0).combine(gt_snpeff, by:0) | i_snpeff
    pfna.combine(pfaa,by:0).combine(pbed,by:0)
      .combine(i_blastn.out,by:0).combine(i_blastp.out,by:0)
      .combine(gt_tandup,by:0) | i_tandup
    chrom_size.combine(chrom_bed,by:0).combine(gap,by:0)
      .combine(ptsv, by:0).combine(pdesc, by:0)
      .combine(gt_rcfg, by:0) | i_rcfg

  emit:
    seq = seqfmt.out.seq
    chrom_size = seqfmt.out.chrom_size
    chrom_bed = seqfmt.out.chrom_bed
    gap = seqfmt.out.gap
    fchain = seqfmt.out.fchain
    rchain = seqfmt.out.rchain
    gff = gff_cln.out
    pgff = gff_idx.out.pgff
}

