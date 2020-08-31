process download {
  label 'low_memory'
  tag "$id"

  input:
  tuple id, species, source, version, assembly, url_fas, url_gff

  output:
  tuple id, path("raw.fasta"), emit: seq
  tuple id, path("raw.gff"), emit: gff

  script:
  species = species.replaceAll('\\s','_')
  assembly = assembly.replaceAll('\\s','_')
  url_pre = "ftp://ftp.ensemblgenomes.org/pub/plants/release-${version.replaceAll(/\.[0-9]+$/,'')}"
  url_fas = source=='ensembl_plants' ? "${url_pre}/fasta/${species.toLowerCase()}/dna/${species}.${assembly}.dna.toplevel.fa.gz" : url_fas
  url_gff = source=='ensembl_plants' ? "${url_pre}/gff3/${species.toLowerCase()}/${species}.${assembly}.${version.replaceAll(/\.[0-9]+$/,'')}.gff3.gz" : url_gff
  if (source == 'local')
    """
    ln -f ${url_fas} raw.fasta
    ln -f ${url_gff} raw.gff
    """
  else
    """
    wget ${url_fas} -O raw.fasta.gz
    wget ${url_gff} -O raw.gff.gz
    gunzip raw.fasta.gz
    gunzip raw.gff.gz
    """
}

process seqfmt {
  label 'low_memory'
  tag "$id"
  publishDir "${params.outdir}/${id}", mode:"copy", overwrite:'true',
    saveAs: {fn ->
      if (fn.indexOf("10.fasta") >= 0) "$fn"
      else "15_intervals/$fn"
      }

  input:
  tuple id, path(fi), gap, chr_prefix

  output:
  tuple id, path("10.fasta"), path("10.fasta.fai"), emit: seq
  tuple id, path("01.chrom.sizes"), emit: chrom_size
  tuple id, path("01.chrom.bed"), emit: chrom_bed
  tuple id, path("11.gap.bed"), emit: gap
  tuple id, path("forward.bed"), path("forward.chain"), emit: fchain
  tuple id, path("reverse.bed"), path("reverse.chain"), emit: rchain

  script:
  def merge_tag = id=='Mt_R108' ? '' : '--merge_short'
  """
  mv $fi raw0.fasta
  fasta.py clean raw0.fasta > raw.fasta
  fasta.py size $fi > raw.sizes

  fasta.py rename $fi renamed.fna forward.bed reverse.bed \\
    --opt ${id} ${merge_tag} \\
    --gap ${gap} --prefix_chr ${chr_prefix}

  fasta.py size renamed.fna > renamed.sizes
  chain.py fromBed forward.bed raw.sizes renamed.sizes > forward.chain
  chainSwap forward.chain reverse.chain
  mv renamed.fna 10.fasta

  samtools faidx 10.fasta
  fasta.py size --bed 10.fasta > 01.chrom.bed
  cut -f1,3 01.chrom.bed > 01.chrom.sizes
  fasta.py gaps 10.fasta > 11.gap.bed
  """
}

process gff_cln {
  label 'low_memory'
  tag "$id"
  publishDir "${params.outdir}/$id/50_annotation", mode:"copy", overwrite:'true'

  input:
  tuple id, path(fi), fixopt, path(fbed), path(fchain)

  output:
  tuple id, path("10.gff")

  script:
  """
  gff.py fix --opt ${fixopt} $fi > 01.fixed.gff
  liftOver -gff 01.fixed.gff ${fchain} 10.gff unmapped
  """
}

process gff_idx {
  label 'mid_memory'
  tag "$id"
  publishDir "${params.outdir}/$id/50_annotation", mode:"copy", overwrite:'true'

  input:
  tuple id, path(ref), path(fai), path(gff)

  output:
  tuple id, path("10.gff.db"), emit: db
  tuple id, path("10.tsv"), emit: tsv
  tuple id, path("10.gtf"), emit: gtf
  tuple id, path("10.bed"), emit: bed
  tuple id, path("10.desc.tsv"), emit: desc
  tuple id, path("10.nt.fasta"), emit: fna
  tuple id, path("10.aa.fasta"), emit: faa
  tuple id, path("15.gff"), emit: pgff
  tuple id, path("15.gff.db"), emit: pdb
  tuple id, path("15.tsv"), emit: ptsv
  tuple id, path("15.gtf"), emit: pgtf
  tuple id, path("15.bed"), emit: pbed
  tuple id, path("15.desc.tsv"), emit: pdesc
  tuple id, path("15.nt.fasta"), emit: pfna
  tuple id, path("15.aa.fasta"), emit: pfaa

  script:
  """
	gff.py index $gff 10.gff.db
	gff.py 2tsv $gff > 10.tsv
	gff.py 2gtf $gff > 10.gtf
	gff.py 2bed12 $gff > 10.bed
	gff.py note --attribute note1,note2 $gff > 10.desc.tsv
	gff.py 2fas $gff $ref > 10.nt.fasta
	fasta.py translate 10.nt.fasta > 10.aa.fasta

	gff.py picklong $gff > 15.gff
	gff.py index 15.gff 15.gff.db
	gff.py 2tsv 15.gff > 15.tsv
	gff.py 2gtf 15.gff > 15.gtf
	gff.py 2bed12 15.gff > 15.bed
	gff.py note --attribute note1,note2 15.gff > 15.desc.tsv
	gff.py 2fas 15.gff $ref > 15.nt.fasta
	fasta.py translate 15.nt.fasta > 15.aa.fasta
	"""
}

// sequence db
process i_blat {
  label 'low_memory'
  tag "$id"
  publishDir "${params.outdir}/$id/21_dbs/blat", mode:"copy", overwrite:'true'

  when: tag == 'T'

  input:
  tuple id, path(ref), path(fai), tag

  output:
  tuple id, path("db.2bit"), path("db.2bit.tile11.ooc")

  script:
  """
  faToTwoBit $ref db.2bit
  blat db.2bit tmp.fas tmp.out -makeOoc=db.2bit.tile11.ooc
  """
}

process i_gatk {
  label 'low_memory'
  tag "$id"
  publishDir "${params.outdir}/$id/21_dbs", mode:"copy", overwrite:'true'

  when: tag == 'T'

  input:
  tuple id, path(ref), path(fai), tag

  output:
  tuple id, val(pre), path("gatk/*")

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
  tag "$id"
  publishDir "${params.outdir}/$id/21_dbs", mode:"copy", overwrite:'true'

  when: tag == 'T'

  input:
  tuple id, path(ref), path(fai), tag

  output:
  tuple id, pre, path("bwa/*")

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
  tag "$id"
  publishDir "${params.outdir}/$id/21_dbs", mode:"copy", overwrite:'true'

  when: tag == 'T'

  input:
  tuple id, path(ref), path(fai), tag

  output:
  tuple id, pre, path("bismark/*")

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
  tag "$id"
  publishDir "${params.outdir}/$id/21_dbs", mode:"copy", overwrite:'true'

  when: tag == 'T'

  input:
  tuple id, path(fna), tag

  output:
  tuple id, pre, path("blastn/*")

  script:
  pre = "blastn/db"
  """
  mkdir blastn
  makeblastdb -dbtype nucl -in $fna -title db -out blastn/db
  """
}

process i_blastp {
  label 'low_memory'
  tag "$id"
  publishDir "${params.outdir}/$id/21_dbs", mode:"copy", overwrite:'true'

  when: tag == 'T'

  input:
  tuple id, path(faa), tag

  output:
  tuple id, pre, path("blastp/*")

  script:
  pre = "blastp/db"
  """
  mkdir blastp
  makeblastdb -dbtype prot -in $faa -title db -out blastp/db
  """
}

process i_lastn {
  label 'low_memory'
  tag "$id"
  publishDir "${params.outdir}/$id/21_dbs", mode:"copy", overwrite:'true'

  when: tag == 'T'

  input:
  tuple id, path(fna), tag

  output:
  tuple id, val(pre), path("lastn/*")

  script:
  pre = "lastn/db"
  """
  mkdir lastn
  lastdb lastn/db $fna
  """
}

process i_lastp {
  label 'low_memory'
  tag "$id"
  publishDir "${params.outdir}/$id/21_dbs", mode:"copy", overwrite:'true'

  when: tag == 'T'

  input:
  tuple id, path(faa), tag

  output:
  tuple id, val(pre), path("lastp/*")

  script:
  pre = "lastp/db"
  """
  mkdir lastp
  lastdb -p lastp/db $faa
  """
}

process i_salmon {
  label 'low_memory'
  tag "$id"
  publishDir "${params.outdir}/$id/21_dbs", mode:"copy", overwrite:'true'

  when: tag == 'T'

  input:
  tuple id, path(fna), tag

  output:
  tuple id, val(pre), path("salmon/*")

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
  tag "$id"
  publishDir "${params.outdir}/$id/21_dbs", mode:"copy", overwrite:'true'

  when: tag == 'T'

  input:
  tuple id, path(ref), path(fai), path(gtf), tag

  output:
  tuple id, pre, path("star/*")

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
  tag "$id"
  publishDir "${params.outdir}/$id/21_dbs", mode:"copy", overwrite:'true'

  when: tag == 'T'

  input:
  tuple id, path(ref), path(fai), path(gtf), tag

  output:
  tuple id, pre, path("hisat2/*")

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
  tag "$id"
  publishDir "${params.outdir}/$id/21_dbs", mode:"copy", overwrite:'true'

  when: tag == 'T'

  input:
  tuple id, path(ref), path(fai), path(gtf), tag

  output:
  tuple id, path("snpeff/$cfg"), path("snpeff/*")

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
  tag "$id"
  publishDir "${params.outdir}/$id/52_tandup", mode:"copy", overwrite:'true'

  when: tag == 'T'

  input:
  tuple id, path(fna), path(faa), path(gbed), blastn_pre, path("blastn/*"), blastp_pre, path("blastp/*"), tag

  output:
  tuple id, path("05.tandup.cds.tsv"), path("05.tandup.pro.tsv")

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
  tag "$id"
  publishDir "${params.outdir}/$id", mode:"copy", overwrite:'true'
  conda '/home/springer/zhoux379/software/miniconda3/envs/r'

  when: tag == 'T'

  input:
  tuple id, path(chrom_size), path(chrom_bed), path(gap_bed), path(ptsv), path(pdes), tag

  output:
  tuple id, path("55.rds")

  script:
  """
  genome.prep.R --dirg ${params.outdir} $id 55.rds
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
        .map {r -> [r.genome, r.species, r.source, r.version, r.assembly, r.url_fas, r.url_gff]}
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

    download(gtable1)
    seq = download.out.seq; gff = download.out.gff
    seqfmt(seq.join(gtable2, by:0))
    seq = seqfmt.out.seq
    chrom_size = seqfmt.out.chrom_size; chrom_bed = seqfmt.out.chrom_bed
    gap = seqfmt.out.gap

    gff_cln(gff.join(gtable3, by:0).join(seqfmt.out.fchain, by:0))
    gff_idx(seq.join(gff_cln.out, by:0))
    gff = gff_cln.out; pgff = gff_idx.out.pgff; gtf = gff_idx.out.gtf
    pfna = gff_idx.out.pfna; pfaa = gff_idx.out.pfaa; pbed = gff_idx.out.pbed
    ptsv = gff_idx.out.ptsv; pdesc = gff_idx.out.pdesc

    seq.join(gt_blat, by:0) | i_blat
    seq.join(gt_gatk, by:0) | i_gatk
    seq.join(gt_bwa, by:0) | i_bwa
    seq.join(gt_bismark, by:0) | i_bismark
    pfna.join(gt_blast, by:0) | i_blastn
    pfna.join(gt_last, by:0) | i_lastn
    pfna.join(gt_salmon, by:0) | i_salmon
    pfaa.join(gt_blast, by:0) | i_blastp
    pfaa.join(gt_last, by:0) | i_lastp
    seq.join(gtf, by:0).join(gt_star, by:0) | i_star
    seq.join(gtf, by:0).join(gt_hisat2, by:0) | i_hisat2
    seq.join(gtf, by:0).join(gt_snpeff, by:0) | i_snpeff
    pfna.join(pfaa,by:0).join(pbed,by:0)
      .join(i_blastn.out,by:0).join(i_blastp.out,by:0)
      .join(gt_tandup,by:0) | i_tandup
    chrom_size.join(chrom_bed,by:0).join(gap,by:0)
      .join(ptsv, by:0).join(pdesc, by:0)
      .join(gt_rcfg, by:0) | i_rcfg

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

