process bamsort {
  label 'mid_memory'
  tag "$id"
  publishDir "${params.outdir}/20_bam", mode:'copy', overwrite:'true',
    saveAs: { fn -> params.saveBAM && params.skip_markdup ? fn : null }

  input:
  tuple val(id), path(ibams)

  output:
  tuple val(id), path("${id}.sorted.bam"), path("${id}.sorted.bam.bai")

  script:
  def suff_mem = ("${(task.memory.toBytes() - 6000000000) / task.cpus}" > 2000000000) ? 'true' : 'false'
  def avail_mem = (task.memory && suff_mem) ? "-m" + "${(task.memory.toBytes() - 6000000000) / task.cpus}" : ''
  """
  bam.py sort ${id}.sorted.bam ${ibams} \\
    --tmpdir ${TMPDIR} --threads ${task.cpus}
  """
}

process markdup {
  label 'mid_memory'
  tag "$id"
  publishDir "${params.outdir}/20_bam", mode:'copy', overwrite:'true',
    saveAs: { fn ->
      if (fn.indexOf(".txt") > 0) "markdup_metrics/$fn"
      else params.saveBAM ? fn : null
    }

  input:
  tuple val(id), path(bam), path(bai)
  
  output:
  tuple val(id), path("${id}.bam"), path("${id}.bam.bai"), emit: bam
  path "${id}.txt", emit: metric

  script:
  markdup_java_options = (task.memory.toGiga() > 8) ? params.markdup_java_options : "\"-Xms" +  (task.memory.toGiga() - 2 )+"g "+ "-Xmx" + (task.memory.toGiga() - 2)+ "g\""
  """
  picard ${markdup_java_options} MarkDuplicates \\
      INPUT=$bam \\
      OUTPUT=${id}.bam \\
      METRICS_FILE=${id}.txt \\
      REMOVE_DUPLICATES=false \\
      ASSUME_SORTED=true \\
      PROGRAM_RECORD_ID='null' \\
      VALIDATION_STRINGENCY=LENIENT
  samtools index ${id}.bam
  """
}

process bamstat {
  label 'process_low'
  tag "$id"
  publishDir "${params.outdir}/20_bam/stats", mode:'copy', overwrite:'true'

  input:
  tuple val(id), path(bam), path(bai)

  output:
  path "${id}.stat.tsv"

  script:
  //samtools flagstat ${id}.bam > ${id}.flagstat.txt
  //samtools idxstats ${id}.bam > ${id}.idxstats.txt
  //samtools stats ${id}.bam > ${id}.samstats.txt
  """
  $baseDir/bin/bam.py stat $bam > ${id}.stat.tsv
  """
}

process pseq {
  label 'mid_memory'
  tag "$id"
  publishDir "${params.outdir}/23_preseq", mode:'copy', overwrite:'true'

  input:
  tuple val(id), path(bam), path(bai), val(spots)

  when:
  !params.skip_qc && !params.skip_preseq && spots >= params.preseq_min_reads

  output:
  path "${id}.ccurve.txt"

  script:
  """
  preseq lc_extrap -v -B $bam -o ${id}.ccurve.txt
  """
}

workflow bam {
  take:
    raw_bams
    ch_read_num
  main:
    raw_bams | bamsort | bamstat
    bams = bamsort.out
    markup = Channel.empty()
    if ( !params.skip_markdup ) {
      bams | markdup
      bams = markdup.out.bam
      markdup = markdup.out.metric
    }
    pseq(bams.join(ch_read_num))
  emit:
    bams = bams
    markdup = markdup
    stats = bamstat.out
    pseq = pseq.out
}



