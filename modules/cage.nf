
process cage1 {
  label "process_medium"
  tag "$id"
  publishDir "${params.outdir}/34_ctss", mode:'copy', overwrite:'true'
  conda "$NXF_CONDA_CACHEDIR/cage"

  input:
  tuple val(id), path(bam), path(bai)

  output:
  tuple val(id), path("${id}.plus.bw"), path("${id}.minus.bw")

  when:
  params.cage

  script:
  mq = params.mapQuality
  genome = params.genome.replaceAll("_",".")
  """
  $baseDir/bin/rnaseq/bam2bw.R $bam -o1 ${id}.plus.bw -o2 ${id}.minus.bw -genome BSgenome.${genome}
  """
}

process cage2 {
  label "process_medium"
  publishDir "${params.outdir}/50_final", mode:'copy', overwrite:'true'
  //conda "$NXF_CONDA_CACHEDIR/cage"

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


