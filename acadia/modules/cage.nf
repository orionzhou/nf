
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

process cage2 {
  label "process_medium"
  publishDir "${params.outdir}/50_final", mode:'copy', overwrite:'true'
  //conda '/home/springer/zhoux379/software/miniconda3/envs/cage'

  input:
  path(bws)
  path(bed)

  output:
  tuple path("metaplot.bed"), path("metaplot.tab"), path("metaplot.tab.gz")

  when:
  params.cage

  script:
  """
  computeMatrix scale-regions -p ${task.cpus} \
    -R $bed \
    -S $bws \
    -b 2500 -a 2500 \
    --binSize 50 \
    --regionBodyLength 5000 \
    --skipZeros -o metaplot.tab.gz \
    --outFileNameMatrix metaplot.tab \
    --outFileSortedRegions metaplot.bed
  """
}


