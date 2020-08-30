
process cage1 {
  label "process_medium"
  tag "$id"
  publishDir "${params.outdir}/34_ctss", mode:'copy', overwrite:'true'
  conda '/home/springer/zhoux379/software/miniconda3/envs/cage'

  input:
  tuple id, path(bam), path(bai)

  output:
  tuple id, path("${id}.plus.bw"), path("${id}.minus.bw")

  when:
  params.cage

  script:
  mq = params.mapQuality
  """
  bam2bw.R $bam -o1 ${id}.plus.bw -o2 ${id}.minus.bw
  """
}



