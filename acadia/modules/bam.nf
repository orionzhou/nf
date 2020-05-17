process bamsort {
  label 'mid_memory'
  tag "$id"
  //publishDir "${params.outdir}/20_bam_sorted", mode:'link', overwrite:'true',
  //  saveAs: { fn -> params.saveBAM ? fn : null }

  input:
  tuple id, path(ibams)

  output:
  tuple id, path("${id}.sorted.bam"), path("${id}.sorted.bam.bai")

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
  publishDir "${params.outdir}/20_bam", mode:'link', overwrite:'true',
    saveAs: { fn ->
      if (fn.indexOf(".txt") > 0) "markdup_metrics/$fn"
      else params.saveBAM ? fn : null
    }

  input:
  tuple id, path(bam), path(bai)

  output:
  tuple id, path("${id}.bam"), path("${id}.bam.bai"), emit: bam
  path "${id}.txt", emit: metric

  script:
  markdup_java_options = (task.memory.toGiga() > 8) ? params.markdup_java_options : "\"-Xms" +  (task.memory.toGiga() / 2 )+"g "+ "-Xmx" + (task.memory.toGiga() - 1)+ "g\""
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
  publishDir "${params.outdir}/20_bam/stats", mode:'link', overwrite:'true'

  input:
  tuple id, path(bam), path(bai)

  output:
  path "${id}.stat.tsv"

  script:
  //samtools flagstat ${id}.bam > ${id}.flagstat.txt
  //samtools idxstats ${id}.bam > ${id}.idxstats.txt
  //samtools stats ${id}.bam > ${id}.samstats.txt
  """
  bam.py stat $bam > ${id}.stat.tsv
  """
}

process pseq {
  label 'mid_memory'
  tag "$id"
  publishDir "${params.outdir}/23_preseq", mode:'link', overwrite:'true'

  input:
  tuple id, path(bam), path(bai), spots

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
    raw_bams | bamsort | (bamstat & markdup)
    bams = markdup.out.bam
    pseq(bams.join(ch_read_num))
  emit:
    bams = bams
    markdup = markdup.out.metric
    stats = bamstat.out
    pseq = pseq.out
}



