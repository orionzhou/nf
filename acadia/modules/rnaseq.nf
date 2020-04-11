process version {
  output:
  path 'software_versions_mqc.yaml', emit: yml
  path "software_versions.csv", emit: csv

  script:
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
  salmon --version &> v_salmon.txt
  picard MarkDuplicates --version &> v_markduplicates.txt  || true
  samtools --version &> v_samtools.txt
  multiqc --version &> v_multiqc.txt
  Rscript -e "library(edgeR); write(x=as.character(packageVersion('edgeR')), file='v_edgeR.txt')"
  Rscript -e "library(dupRadar); write(x=as.character(packageVersion('dupRadar')), file='v_dupRadar.txt')"
  unset DISPLAY && qualimap rnaseq &> v_qualimap.txt || true
  scrape_software_versions.py &> software_versions_mqc.yaml
  """
}

process fqc {
  label 'mid_memory'
  tag "${params.name}.$name"

  when:
  !params.skipQC && !params.skipFastQC

  input:
  tuple name, paired, path(reads)

  output:
  path "*_fastqc.zip", emit: zip
  path "*_fastqc.html", emit: html

  script:
  """
  fastqc --quiet --threads $task.cpus $reads
  """
}

process trim {
  label 'low_memory'
  tag "${params.name}.$name"

  input:
  tuple name, paired, path(reads)

  output:
  tuple name, paired, path("*fq.gz"), emit: reads
  path "*trimming_report.txt", emit: log
  path "*_fastqc.{html,zip}", emit: fastqc

  script:
  c_r1 = params.clip_r1 > 0 ? "--clip_r1 ${clip_r1}" : ''
  c_r2 = params.clip_r2 > 0 ? "--clip_r2 ${clip_r2}" : ''
  tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${three_prime_clip_r1}" : ''
  tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${three_prime_clip_r2}" : ''
  nextseq = params.trim_nextseq > 0 ? "--nextseq ${params.trim_nextseq}" : ''
  pair_str = paired ? "--paired" : ""
  clip_str = paired ? "$c_r1 $c_r2" : c_r1
  tpc_str = paired ? "$tpc_r1 $tpc_r2" : tpc_r1
  """
  trim_galore $pair_str --fastqc --gzip $clip_str $tpc_str $nextseq $reads
  """
}

process hs2 {
  label 'high_memory'
  tag "${params.name}.$name"

  input:
  tuple name, paired, path(reads), pre, path(hs2_indices), path(splicesites)

  output:
  tuple name, path("${name}.bam"), emit: bam
  path "${name}.hisat2_summary.txt", emit: log
  path "unmapped.hisat2*" optional true

  script:
  seq_center = params.seq_center ? 'CN:$params.seq_center' : ''
  rg = "--rg-id ${name} --rg ${seq_center.replaceAll('\\s','_')} SM:${name}"
  def strandness = ''
  if (params.stranded == 'forward') {
    strandness = !paired ? '--rna-strandness F' : '--rna-strandness FR'
  } else if (params.stranded == 'reverse') {
    strandness = !paired ? '--rna-strandness R' : '--rna-strandness RF'
  }
  input = paired ? "-1 ${reads[0]} -2 ${reads[1]}" : "-U ${reads}"
  unaligned = params.saveUnaligned ? !paired ? "--un-gz unmapped.hisat2.gz" : "--un-conc-gz unmapped.hisat2.gz" : ''
  // --dta --no-mixed --no-discordant
  // --known-splicesite-infile $splicesites \\
  """
  hisat2 -x ${params.hisat2_index} \\
     ${input} $strandness \\
     -p ${task.cpus} ${unaligned} \\
     --met-stderr --new-summary \\
     --summary-file ${name}.hisat2_summary.txt $rg \\
     | samtools view -bSh -o ${name}.bam -
  """
}

process star {
  label 'high_memory'
  tag name

  input:
  tuple name, paired, path(reads), path(index), path(gtf)

  output:
  tuple path("*Log.final.out"), path('*.bam'), emit: bam
  path "*.out", emit: log
  path "*SJ.out.tab"
  path "*Log.out"
  path "*Unmapped*" optional true

  script:
  def avail_mem = "--limitBAMsortRAM ${task.memory.toBytes() - 100000000}"
  seq_center = params.seq_center ? 'CN:$params.seq_center' : ''
  unaligned = params.saveUnaligned ? "--outReadsUnmapped Fastx" : ''
  extra = """
      --outSAMmapqUnique 60 \\
      --outFilterType BySJout \\
      --outFilterMultimapNmax 10 \\
      --outFilterMismatchNmax 999 \\
      --outFilterMismatchNoverLmax 1 \\
      --outFilterMismatchNoverReadLmax 1.0 \\
      --outFilterMatchNminOverLread 0 \\
      --outFilterMatchNmin 0 \\
      --outFilterScoreMinOverLread 0 \\
      --alignSJoverhangMin 8 \\
      --alignSJDBoverhangMin 1 \\
      --alignIntronMin 20 \\
      --alignIntronMax 1000000 \\
      --alignMatesGapMax 1000000 \\
      --outSAMunmapped Within KeepPairs \\
      --outSAMattrRGline ID:$name $seq_center 'SM:$name'
      """.stripIndent()
  """
  STAR --genomeDir $index \\
    --sjdbGTFfile $gtf \\
    --runThreadN ${task.cpus} \\
    --readFilesIn $reads  \\
    --readFilesCommand zcat \\
    --twopassMode Basic \\
    --outWigType bedGraph \\
    --outSAMtype SAM Unsorted $avail_mem \\
    --runDirPerm All_RWX $unaligned \\
    --outFileNamePrefix $name \\
    --outStd Log $extra
  """
}

process bamsort {
  label 'mid_memory'
  tag "${params.name}.$name"

  input:
  tuple name, path(ibams)

  output:
  tuple name, path("${name}.sorted.bam"), path("${name}.sorted.bam.bai")

  script:
  def suff_mem = ("${(task.memory.toBytes() - 6000000000) / task.cpus}" > 2000000000) ? 'true' : 'false'
  def avail_mem = (task.memory && suff_mem) ? "-m" + "${(task.memory.toBytes() - 6000000000) / task.cpus}" : ''
  """
  bam.py sort ${name}.sorted.bam ${ibams} \\
    --tmpdir ${TMPDIR} --threads ${task.cpus}
  """
}

process bamstat {
  label 'process_low'
  tag "${params.name}.$name"

  input:
  tuple name, path(bam), path(bai)

  output:
  path "${name}.tsv"

  script:
  """
  bam.py stat $bam > ${name}.tsv
  """
}

process rseqc {
  label 'mid_memory'
  tag "${params.name}.$name"

  when:
  !params.skipQC && !params.skipRseQC

  input:
  tuple name, path(bam), path(bai), path(bed12)

  output:
  path "*bam_stat.txt", emit: bamstat
  path "*infer_experiment.txt", emit: infer
  path "*read_distribution.txt", emit: distr
  path "*read_duplication.DupRate_plot.pdf", emit: dup
  path "*read_duplication.DupRate_plot.r", emit: dup_r
  path "*read_duplication.pos.DupRate.xls", emit: dup_pos
  path "*read_duplication.seq.DupRate.xls", emit: dup_seq
  path "*RPKM_saturation.saturation.pdf", emit: sat optional true
  path "*RPKM_saturation.saturation.r", emit: sat_r optional true
  path "*RPKM_saturation.eRPKM.xls", emit: sat_rpkm optional true
  path "*RPKM_saturation.rawCount.xls", emit: sat_cnt optional true
  path "*inner_distance.txt", emit: inndst optional true
  path "*inner_distance_freq.txt", emit: inndst_data optional true
  path "*inner_distance_plot.r", emit: inndst_r optional true
  path "*inner_distance_plot.pdf", emit: inndst_plot optional true
  path "*junction_annotation_log.txt", emit: jct
  path "*junction_plot.r", emit: jct_r
  path "*junction.xls", emit: jct_data
  path "*splice_events.pdf", emit: jct_event
  path "*splice_junction.pdf", emit: jct_jct
  path "*junctionSaturation_plot.pdf", emit: jctsat
  path "*junctionSaturation_plot.r", emit: jctsat_r

  script:
  """
  infer_experiment.py -i $bam -r $bed12 > ${name}.infer_experiment.txt
  junction_annotation.py -i $bam -o ${name}.rseqc -r $bed12 2> ${name}.junction_annotation_log.txt
  bam_stat.py -i $bam 2> ${name}.bam_stat.txt
  junction_saturation.py -i $bam -o ${name}.rseqc -r $bed12
  inner_distance.py -i $bam -o ${name}.rseqc -r $bed12
  read_distribution.py -i $bam -r $bed12 > ${name}.read_distribution.txt
  read_duplication.py -i $bam -o ${name}.read_duplication
  """
}

process pseq {
  label 'mid_memory'
  tag "${params.name}.$name"

  input:
  tuple name, path(bam), path(bai), spots

  when:
  !params.skipQC && !params.skipPreseq && spots >= params.preseq_min_reads

  output:
  path "${name}.ccurve.txt"

  script:
  """
  preseq lc_extrap -v -B $bam -o ${name}.ccurve.txt
  """
}

process markdup {
  label 'mid_memory'
  tag "${params.name}.$name"

  when:
  !params.skipQC && !params.skipDupRadar

  input:
  tuple name, path(bam), path(bai)

  output:
  tuple name, path("${name}.bam"), path("${name}.bam.bai"), emit: bam
  path "${name}.txt", emit: metric

  script:
  markdup_java_options = (task.memory.toGiga() > 8) ? params.markdup_java_options : "\"-Xms" +  (task.memory.toGiga() / 2 )+"g "+ "-Xmx" + (task.memory.toGiga() - 1)+ "g\""
  """
  picard ${markdup_java_options} MarkDuplicates \\
      INPUT=$bam \\
      OUTPUT=${name}.bam \\
      METRICS_FILE=${name}.txt \\
      REMOVE_DUPLICATES=false \\
      ASSUME_SORTED=true \\
      PROGRAM_RECORD_ID='null' \\
      VALIDATION_STRINGENCY=LENIENT
  samtools index ${name}.bam
  """
}

process qmap {
  label 'low_memory'
  tag "${params.name}.$name"

  when:
  !params.skipQC && !params.skipQualimap

  input:
  tuple name, path(bam), path(bai), path(gtf), paired, path(reads)

  output:
  path "${name}"

  script:
  def drc = params.stranded ? params.stranded=='reverse' ? 'strand-specific-reverse' : 'strand-specific-forward' : 'non-strand-specific'
  def pair_str = paired ? '-pe' : ''
  memory = task.memory.toGiga() + "G"
  """
  unset DISPLAY
  qualimap --java-mem-size=${memory} rnaseq -p $drc $pair_str \\
    -bam $bam -gtf $gtf -outdir ${name}
  """
}

process duprad {
  label 'low_memory'
  tag "${params.name}.$name"

  when:
  !params.skipQC && !params.skipDupRadar

  input:
  tuple name, path(bam), path(bai), path(gtf), paired, path(reads)

  output:
  path "*_duprateExpDens.pdf", emit: expDen
  path "*_duprateExpBoxplot.pdf", emit: expBox
  path "*_expressionHist.pdf", emit: hist
  path "*_dupMatrix.txt", emit: data
  path "*_duprateExpDensCurve_mqc.txt", emit: curve
  path "*_intercept_slope.txt", emit: slope

  script: // This script is bundled with the pipeline, in bin/
  def drc = params.stranded ? params.stranded=='reverse' ? 2:1:0
  def pair_str = paired ? 'paired':'single'
  """
  dupRadar.r $bam $gtf $drc $pair_str ${task.cpus}
  """
}

process fcnt {
  label 'low_memory'
  tag "${params.name}.$name"

  input:
  tuple name, path(bam), path(bai), path(gtf), path(biotypes_header)

  output:
  path "${name}_gene.featureCounts.txt", emit: txt
  path "${name}_gene.featureCounts.txt.summary", emit: log
  path "${name}_biotype_counts*mqc.{txt,tsv}", emit: biotype optional true

  script:
  def extra = params.fc_extra_attributes ? "--extraAttributes ${params.fc_extra_attributes}" : ''
  def drc = params.stranded ? params.stranded=='reverse' ? 2 : 1 : 0
  biotype_qc = params.skipBiotypeQC ? '' : """
    featureCounts -a $gtf -g ${params.biotype} \\
    -o ${name}_biotype.featureCounts.txt -p \\
    -s $drc $bam
    """.stripIndent()
  mod_biotype = params.skipBiotypeQC ? '' : """
    cut -f 1,7 ${name}_biotype.featureCounts.txt | \\
    tail -n +3 | cat $biotypes_header - >> ${name}_biotype_counts_mqc.txt &&\\
    mqc_features_stat.py ${name}_biotype_counts_mqc.txt -s $name \\
    -f rRNA -o ${name}_biotype_counts_gs_mqc.tsv
    """.stripIndent()
  """
  featureCounts -a $gtf -Q ${params.mapQuality} -T ${task.cpus} \\
    -g ${params.fc_group_features} -t ${params.fc_count_type} \\
    -o ${name}_gene.featureCounts.txt \\
    $extra -p -s $drc $bam
  $biotype_qc
  $mod_biotype
  """
}

process stie {
  label "mid_memory"
  tag "${params.name}.$name"

  input:
  tuple name, path(bam), path(bai), path(gtf)

  output:
  path "${name}_transcripts.gtf", emit: gtf
  path "${name}.gene_abund.txt", emit: abund
  path "${name}.cov_refs.gtf", emit: cov
  path "${name}_ballgown", emit: ballgown

  script:
  def st_direction = params.stranded ? params.stranded=='reverse' ? "--rf" : "--fr" : ''
  def ignore_gtf = params.stringTieIgnoreGTF ? "" : "-e"
  """
  stringtie $bam $st_direction \\
      -o ${name}_transcripts.gtf \\
      -v -G $gtf \\
      -A ${name}.gene_abund.txt \\
      -C ${name}.cov_refs.gtf \\
      -b ${name}_ballgown \\
      $ignore_gtf
  """
}

process corr {
  label 'low_memory'
  tag "${params.name}"

  when:
  !params.skipQC && !params.skipEdgeR

  input:
  path fcnts
  val num_bams
  path mdsplot_header
  path heatmap_header

  output:
  path "*.{txt,pdf,csv}"

  when:
  num_bams >= 3 && (!params.sampleLevel)

  script: // This script is bundled with the pipeline, in bin/
  """
  edgeR_heatmap_MDS.r $fcnts
  cat $mdsplot_header edgeR_MDS_Aplot_coordinates_mqc.csv >> tmp_file
  mv tmp_file edgeR_MDS_Aplot_coordinates_mqc.csv
  cat $heatmap_header log2CPM_sample_correlation_mqc.csv >> tmp_file
  mv tmp_file log2CPM_sample_correlation_mqc.csv
  """
}

process salm {
  label 'mid_memory'
  tag "${params.name}.$name"

  input:
  tuple name, paired, path(reads), path(index), path(gtf), path("tx2gene.csv")

  output:
  path "${name}_salmon_gene_counts.csv", emit: gcnt
  path "${name}_salmon_gene_tpm.csv", emit: gtpm
  path "${name}_salmon_transcript_counts.csv", emit: tcnt
  path "${name}_salmon_transcript_tpm.csv", emit: ttpm

  script:
  def rnastrandness = !paired ? 'U' : 'IU'
  if (params.stranded == 'forward') {
      rnastrandness = !paired ? 'SF' : 'ISF'
  } else if (params.stranded == 'reverse') {
      rnastrandness = !paired ? 'SR' : 'ISR'
  }
  def input = paired ? "-1 ${reads[0]} -2 ${reads[1]}" : "-r ${reads[0]}"
  unmapped = params.saveUnaligned ? "--writeUnmappedNames" : ''
  """
  salmon quant --validateMappings \\
               --seqBias --useVBOpt --gcBias \\
               --geneMap ${gtf} \\
               --threads ${task.cpus} \\
               --libType=${rnastrandness} \\
               --index ${index} \\
               $input $unmapped\\
               -o ${name}
  tximport.r NULL ${name} ${name}
  """
}

process ase1 {
  label "process_medium"
  tag "${params.name}.$name"
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
  tag "${params.name}.$name"

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

process outdoc {
  tag "${params.name}"

  input:
  path output_docs

  output:
  path "results_description.html"

  script:
  """
  markdown_to_html.r $output_docs results_description.html
  """
}

process mg {
  label "mid_memory"
  tag "$yid"
  conda '/home/springer/zhoux379/software/miniconda3/envs/r'

  input:
  path meta
  path rcfg
  path bamstats
  path fcnts
  path salm_gene_cnt
  path salm_gene_tpm
  path salm_tx_cnt
  path salm_tx_tpm
  path ase_gene
  path ase_snp
  path ril
  path ril_txt

  output:
  path "00.raw.rds"

  script:
  yid = params.name
  //rc2cpm.R $meta featurecounts.rds cpm.rds --opt featurecounts --config $cfg
  cmds = []
  opts = []
  if (params.ase) {
    cmds += "merge.stats.R --opt ase_gene -o ase_gene.rds $ase_gene"
    cmds += "merge.stats.R --opt ase_snp -o ase_snp.rds $ase_snp"
    opts += '--ase_gene ase_gene.rds'
    opts += '--ase_snp ase_snp.rds'
  }
  if (params.ril) {
    cmds += "merge.stats.R --opt snpbinner -o ril.rds $ril"
    opts += '--ril ril.rds'
  }
  cmd = cmds.join("\n")
  opt = opts.join(" ")
  """
  merge.stats.R --opt bam_stat -o bamstats.tsv $bamstats
  merge.stats.R --opt featurecounts -o featurecounts.rds $fcnts
  merge.stats.R --opt salmon -o salmon_gene_cnt.rds $salm_gene_cnt
  merge.stats.R --opt salmon -o salmon_gene_tpm.rds $salm_gene_tpm
  merge.stats.R --opt salmon -o salmon_tx_cnt.rds $salm_tx_cnt
  merge.stats.R --opt salmon -o salmon_tx_tpm.rds $salm_tx_tpm
  $cmd

  nf_rnaseq.R $meta merged.rds --bamstat bamstats.tsv \\
    --fcnt featurecounts.rds \\
    --salmon_gcnt salmon_gene_cnt.rds --salmon_gtpm salmon_gene_tpm.rds \\
    --salmon_tcnt salmon_tx_cnt.rds   --salmon_ttpm salmon_tx_tpm.rds \\
    $opt
  nf_rnaseq_norm.R merged.rds 00.raw.rds --rcfg $rcfg
  """
}

process renorm {
  label 'process_medium'
  tag "${params.name}"
  conda '/home/springer/zhoux379/software/miniconda3/envs/r'

  input:
  path meta
  path raw_rds
  path rcfg

  output:
  path '01.rds'

  when:
  meta.exists()

  script:
  """
  nf_rnaseq_norm.R --meta $meta $raw_rds 01.rds --rcfg $rcfg
  """
}

process mqc {
  label 'low_memory'
  tag "${params.name}"

  when:
  !params.skipMultiQC

  input:
  path multiqc_config
  path ('fastqc/*')
  path ('trimgalore/*')
  path ('alignment/*')
  path ('rseqc/*')
  path ('preseq/*')
  path ('qualimap/*')
  path ('dupradar/*')
  path ('featureCounts/*')
  path ('featureCounts_biotype/*')
  //path ('salmon/*')
  path ('sample_correlation_results/*')
  //path ('sortmerna/*')
  path ('software_versions/*')
  path workflow_summary

  output:
  path "*.html", emit: html
  path "*_data", emit: data
  path "multiqc_plots", emit: plot

  script:
  def run = params.name
  rtitle = "--title \"$run\""
  rfilename = "--filename " + run.replaceAll('\\W','_').replaceAll('_+','_')
  """
  multiqc . -f $rtitle $rfilename --config $multiqc_config \\
      -m custom_content -m picard -m preseq -m rseqc -m featureCounts \\
      -m hisat2 -m star -m cutadapt -m sortmerna -m fastqc -m qualimap -m salmon
  """
}





