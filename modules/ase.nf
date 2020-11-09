
process ase1 {
  label "process_medium"
  tag "$name"
  //conda "$NXF_CONDA_CACHEDIR/alfred"

  input:
  tuple val(name), path(bam), path(bai), path(bcf), path(csi), path(ref)

  output:
  tuple val(name), path("${name}.tsv.gz"), path("${name}.1.bam"), path("${name}.2.bam")

  when:
  params.ase && bcf.exists()

  script:
  mq = params.mapQuality
  """
  echo hkk
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
  tuple val(name), path(tsv), path(bam1), path(bam2), path(gtf)

  output:
  path "$tsv", emit: snp
  path "${name}.?.tsv", emit: gene

  script:
  def mq = params.mapQuality
  def flag_attr = params.fc_extra_attributes ? "--extraAttributes ${params.fc_extra_attributes}" : ''
  def flag_srd = params.stranded == 'no' ? 0 : params.stranded=='reverse' ? 2 : 1
  def flag_pe = ''//paired == 'PE' ? "-p" : ""
  def flag_long = params.read_type == 'nanopore' ? "-L" : ""
  extra = "--primary -Q $mq -t exon -g gene_id --byReadGroup -s $flag_srd"
  """
  featureCounts -a $gtf --primary -Q $mq -T ${task.cpus} \\
    -g ${params.fc_group_features} -t ${params.fc_count_type} \\
    $flag_attr $flag_pe $flag_long -s $flag_srd \\
    -o ${name}.1.tsv $bam1
  featureCounts -a $gtf --primary -Q $mq -T ${task.cpus} \\
    -g ${params.fc_group_features} -t ${params.fc_count_type} \\
    $flag_attr $flag_pe $flag_long -s $flag_srd \\
    -o ${name}.2.tsv $bam2
  """
}


