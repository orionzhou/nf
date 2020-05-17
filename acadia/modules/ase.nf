
process ase1 {
  label "process_medium"
  tag "$name"
  conda '/home/springer/zhoux379/software/miniconda3/envs/alfred'

  input:
  tuple name, path(bam), path(bai), path(bcf), path(csi), path(ref)

  output:
  tuple name, path("${name}.tsv.gz"), path("${name}.1.bam"), path("${name}.2.bam")

  when:
  params.ase && bcf.exists()

  script:
  mq = params.mapQuality
  """
  alfred ase -r $ref -s sample -v $bcf -m $mq -p -a ${name}.tsv.gz $bam
  alfred split -r $ref -s sample -v $bcf -m $mq \\
      -p ${name}.1.bam -q ${name}.2.bam $bam
  """
}

process ase2 {
  label "process_medium"
  tag "$name"
  //publishDir "${params.outdir}/37_ase", mode:'link', overwrite:'true'

  when:
  params.ase

  input:
  tuple name, path(tsv), path(bam1), path(bam2), path(gtf)

  output:
  path "$tsv", emit: snp
  path "${name}.?.tsv", emit: gene

  script:
  mq = params.mapQuality
  def extraAttributes = params.fc_extra_attributes ? "--extraAttributes ${params.fc_extra_attributes}" : ''
  def featureCounts_direction = params.stranded ? params.stranded=='reverse' ? 2 : 1 : 0
  extra = "--primary -Q $mq -t exon -g gene_id --byReadGroup -s $featureCounts_direction"
  """
  featureCounts -a $gtf -Q ${params.mapQuality} -T ${task.cpus} \\
      -g ${params.fc_group_features} -t ${params.fc_count_type} \\
      $extraAttributes -p -s $featureCounts_direction \\
      -o ${name}.1.tsv $bam1
  featureCounts -a $gtf -Q ${params.mapQuality} -T ${task.cpus} \\
      -g ${params.fc_group_features} -t ${params.fc_count_type} \\
      $extraAttributes -p -s $featureCounts_direction \\
      -o ${name}.2.tsv $bam2
  """
}


