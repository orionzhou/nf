process cage_gtf {
  label "low_memory"
  tag "$genome"
  publishDir "${params.outdir}/30_gtf", mode:'copy', overwrite:'true',
    saveAs: { fn -> "${genome}.gtf"}

  when:
  params.cage

  input:
  tuple val(genome), path(gtf), path(sizes)

  output:
  tuple val(genome), path("gene.gtf")

  script:
	left = 1000
	right = 0
  """
	bioawk -t '{if(\$3=="gene") {split(\$9, a, ";"); split(a[1],b,"="); print \$1, \$4-1, \$5, b[2], ".", \$7}}' $gtf > gene.bed
	bedtools slop -i gene.bed -g $sizes -s -l $left -r $right > gene.2.bed
	bioawk -t '{print \$1, ".", "gene", \$2+1, \$3, ".", \$6, ".", "gene_id \\""\$4"\\";"}' gene.2.bed > gene.gtf
  """
}

process fcnt {
  label 'low_memory'
  tag "$id"
  publishDir "${params.outdir}/31_featureCounts", mode:'copy', overwrite:'true',
    saveAs: {fn ->
      if (fn.indexOf(".tsv") > 0) "$fn"
      else if (fn.indexOf(".summary.txt") > 0) "summary/$fn"
      else fn
    }

  input:
  tuple val(sid), val(genome), path(bam), path(bai), path(gtf), val(paired), path(reads)

  output:
  path "${id}.tsv", emit: tsv
  path "${id}.summary.txt", emit: log

  script:
  id = "${sid}-${genome}"
  mq = params.mapQuality
  flag_srd = params.stranded == 'no' ? 0 : params.stranded=='reverse' ? 2 : 1
  flag_pe = paired == 'PE' ? "-p" : ""
  flag_long = params.read_type == 'nanopore' ? "-L" : ""
  flag_cage = params.cage ? "--read2pos 5" : ""
  flag_multi = params.count_multi ? "-M -O --fraction" : ""
  """
  featureCounts --primary -Q $mq -T ${task.cpus} \\
    -a $gtf -g gene_id -t gene \\
    $flag_pe $flag_long -s $flag_srd $flag_multi $flag_cage \\
    -o ${id}.tsv \\
    $bam
  mv ${id}.tsv.summary ${id}.summary.txt
  """
}

process stie {
  label "mid_memory"
  tag "$id"
  publishDir "${params.outdir}/32_stringtieFPKM", mode:'copy', overwrite:'true',
    saveAs: { fn ->
      if (fn.indexOf("transcripts.gtf") > 0) "transcripts/$fn"
      else if (fn.indexOf("cov_refs.gtf") > 0) "cov_refs/$fn"
      else if (fn.indexOf("ballgown") > 0) "ballgown/$fn"
      else "$fn"
    }

  when:
  params.run_stringtie

  input:
  tuple val(id), val(genome), path(bam), path(bai), path(gtf)

  output:
  path "${id}_transcripts.gtf", emit: gtf
  path "${id}.gene_abund.txt", emit: abund
  path "${id}.cov_refs.gtf", emit: cov
  path "${id}_ballgown", emit: ballgown

  script:
  id = "${sid}-${genome}"
  def st_direction = params.stranded == 'no' ? '' : params.stranded=='reverse' ? "--rf" : "--fr"
  def ignore_gtf = params.stringtie_ignore_gtf ? "" : "-e"
  """
  stringtie $bam $st_direction \\
      -o ${id}_transcripts.gtf \\
      -v -G $gtf \\
      -A ${id}.gene_abund.txt \\
      -C ${id}.cov_refs.gtf \\
      -b ${id}_ballgown \\
      $ignore_gtf
  """
}

process salm {
  label 'mid_memory'
  tag "$id"
  publishDir "${params.outdir}/34_salmon", mode:'copy', overwrite:'true'

  when:
  params.run_salmon

  input:
  tuple val(sid), val(genome), val(paired), path(reads), path(index), path(gtf), path("tx2gene.csv")

  output:
  id = "${sid}-${genome}"
  path "${id}_salmon_gene_counts.csv", emit: gcnt
  path "${id}_salmon_gene_tpm.csv", emit: gtpm
  path "${id}_salmon_transcript_counts.csv", emit: tcnt
  path "${id}_salmon_transcript_tpm.csv", emit: ttpm

  script:
  def rnastrandness = paired == 'SE' ? 'U' : 'IU'
  if (params.stranded == 'forward') {
      rnastrandness = paired == 'SE' ? 'SF' : 'ISF'
  } else if (params.stranded == 'reverse') {
      rnastrandness = paired == 'SE' ? 'SR' : 'ISR'
  }
  def input = paired=='SE' ? "-r ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
  unmapped = params.save_unmapped ? "--writeUnmappedNames" : ''
  """
  salmon quant --validateMappings \\
               --seqBias --useVBOpt --gcBias \\
               --geneMap ${gtf} \\
               --threads ${task.cpus} \\
               --libType=${rnastrandness} \\
               --index ${index} \\
               $input $unmapped\\
               -o ${id}
  $baseDir/bin/rnaseq/tximport.r NULL ${id} ${id}
  """
}

process bigwig {
  label 'mid_memory'
  tag "$id"
  publishDir "${params.outdir}/23_bigwig", mode:'copy', overwrite:'true'

  input:
  tuple val(sid), val(genome), path(bam), path(bai)

  output:
  tuple val(id), path("${id}.bw")

  script:
  id = "${sid}-${genome}"
  binsize = 50
  //--effectiveGenomeSize ${egs} 
  """
  bamCoverage --binSize ${binsize} -p ${task.cpus} -b $bam -o ${id}.bw
  """
}

workflow rnaseq {
  take:
    bams0
    reads
    design
    genomes
  main:
    // channel set up
      sizes = genomes
        .map {r -> [r[0], file(r[1].genome_sizes, checkIfExists: true)]}
      gff = genomes
        .map {r -> [r[0], file(r[1].gff, checkIfExists: true)]}
      gtf = genomes
        .map {r -> [r[0], file(r[1].gtf, checkIfExists: true)]}
      bed = genomes
        .map {r -> [r[0], file(r[1].bed, checkIfExists: true)]}
      pbed = genomes
        .map {r -> [r[0], file(r[1].pbed, checkIfExists: true)]}
		
		if (params.cage) {
			cage_gtf(gff.join(sizes))
      gtf = cage_gtf.out
		}
    bams = bams0
      .map {r -> [r[0].split("-")[0], r[0].split("-")[1], r[1], r[2]]}
    in1 = bams
      .map {r -> [r[1], r[0], r[2], r[3]]}
      .combine(gtf, by:0)
      .map {r -> [r[1], r[0], r[2], r[3], r[4]]}
      .combine(reads, by:0)
    in1 | fcnt

  emit:
    fcnt_tsv = fcnt.out.tsv
    fcnt_log = fcnt.out.log
}

process mg {
  label "mid_memory"
  tag "${params.name}"
  //conda "$NXF_CONDA_CACHEDIR/r"
  publishDir "${params.outdir}/50_final", mode:'copy', overwrite:'true'
  //publishDir "${params.qcdir}/${params.genome}/${params.name}", mode:'copy', overwrite:'true'

  input:
  path meta
  path bamstats
  path fcnts

  output:
  path "bamstats.tsv"
  path "featurecounts.rds"

  script:
  yid = params.name
  cmds = []
  opts = []
  """
  $baseDir/bin/merge.stats.R --opt bam_stat -o bamstats.tsv $bamstats
  $baseDir/bin/merge.stats.R --opt featurecounts -o featurecounts.rds $fcnts
  """
}


