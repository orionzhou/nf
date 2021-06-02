process version {
  publishDir "${params.outdir}/pipeline_info", mode:'copy', overwrite:'true'

  output:
  path 'software_versions_mqc.yaml', emit: yml
  path "software_versions.csv", emit: csv

  script:
  //salmon --version &> v_salmon.txt
  """
  echo $workflow.manifest.version &> v_ngi_rnaseq.txt
  echo $workflow.nextflow.version &> v_nextflow.txt
  fastqc --version &> v_fastqc.txt
  cutadapt --version &> v_cutadapt.txt
  trim_galore --version &> v_trim_galore.txt
  sortmerna --version &> v_sortmerna.txt
  STAR --version &> v_star.txt
  hisat2 --version &> v_hisat2.txt
  stringtie --version &> v_stringtie.txt
  preseq &> v_preseq.txt
  read_duplication.py --version &> v_rseqc.txt
  bamCoverage --version &> v_deeptools.txt || true
  featureCounts -v &> v_featurecounts.txt
  picard MarkDuplicates --version &> v_markduplicates.txt  || true
  samtools --version &> v_samtools.txt
  multiqc --version &> v_multiqc.txt
  Rscript -e "library(edgeR); write(x=as.character(packageVersion('edgeR')), file='v_edgeR.txt')"
  Rscript -e "library(dupRadar); write(x=as.character(packageVersion('dupRadar')), file='v_dupRadar.txt')"
  unset DISPLAY && qualimap rnaseq &> v_qualimap.txt || true
  $baseDir/bin/rnaseq/versions.py &> software_versions_mqc.yaml
  """
}

process fcnt {
  label 'low_memory'
  tag "$id"
  publishDir "${params.outdir}/31_featureCounts", mode:'copy', overwrite:'true',
    saveAs: {fn ->
      if (fn.indexOf("biotype_counts") > 0) "biotype_counts/$fn"
      else if (fn.indexOf("_gene.featureCounts.txt.summary") > 0) "gene_count_summaries/$fn"
      else if (fn.indexOf("_gene.featureCounts.txt") > 0) "gene_counts/$fn"
      else fn
    }

  input:
  tuple val(sid), val(genome), path(bam), path(bai), path(gtf), path(biotypes_header), val(paired), path(reads)

  output:
  path "${id}_gene.featureCounts.txt", emit: txt
  path "${id}_gene.featureCounts.txt.summary", emit: log
  path "${id}_biotype_counts*mqc.{txt,tsv}", emit: biotype optional true

  script:
  id = "${sid}-${genome}"
  def mq = params.mapQuality
  def flag_attr = params.fc_extra_attributes ? "--extraAttributes ${params.fc_extra_attributes}" : ''
  def flag_srd = params.stranded == 'no' ? 0 : params.stranded=='reverse' ? 2 : 1
  def flag_pe = paired == 'PE' ? "-p" : ""
  def flag_long = params.read_type == 'nanopore' ? "-L" : ""
  biotype_qc = params.skip_biotype_qc ? '' : """
    featureCounts -a $gtf -g ${params.biotype} \\
    -o ${id}_biotype.featureCounts.txt -p \\
    -s $flag_srd $bam
    """.stripIndent()
  mod_biotype = params.skip_biotype_qc ? '' : """
    cut -f 1,7 ${id}_biotype.featureCounts.txt | \\
    tail -n +3 | cat $biotypes_header - >> ${id}_biotype_counts_mqc.txt &&\\
    mqc_features_stat.py ${id}_biotype_counts_mqc.txt -s $id \\
    -f rRNA -o ${id}_biotype_counts_gs_mqc.tsv
    """.stripIndent()
  //$biotype_qc
  //$mod_biotype
  """
  featureCounts -a $gtf --primary -Q $mq -T ${task.cpus} \\
    -g ${params.fc_group_features} -t ${params.fc_count_type} \\
    -o ${id}_gene.featureCounts.txt \\
    $flag_attr $flag_pe $flag_long -s $flag_srd \\
    $bam
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

include {cage1; cage2} from './cage.nf'

workflow rnaseq {
  take:
    bams0
    reads
    design
    genomes
  main:
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

    bams = bams0
      .map {r -> [r[0].split("-")[0], r[0].split("-")[1], r[1], r[2]]}
    in1 = bams
      .map {r -> [r[1], r[0], r[2], r[3]]}
      .join(gtf)
      .map {r -> [r[1], r[0], r[2], r[3], r[4]]}
      .join(reads)
    fcnt(in1)

  emit:
    fcnt_txt = fcnt.out.txt
    fcnt_log = fcnt.out.log
}

process mg {
  label "mid_memory"
  tag "${params.name}"
  //conda "$NXF_CONDA_CACHEDIR/r"
  publishDir "${params.outdir}/50_final", mode:'copy', overwrite:'true'
  publishDir "${params.qcdir}/${params.genome}/${params.name}", mode:'copy', overwrite:'true'

  input:
  path meta
  path bamstats
  path fcnts

  output:
  path "meta.tsv"
  path "bamstats.tsv"
  path "featurecounts.rds"

  script:
  yid = params.name
  cmds = []
  opts = []
  """
  cp $meta meta.tsv
  $baseDir/bin/merge.stats.R --opt bam_stat -o bamstats.tsv $bamstats
  $baseDir/bin/merge.stats.R --opt featurecounts -o featurecounts.rds $fcnts
  """
}




