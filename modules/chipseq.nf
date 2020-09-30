process version {
  label 'low_memory'
  tag "${params.name}"

  output:
  path 'software_versions_mqc.yaml', emit: yml
  path "software_versions.csv", emit: csv

  script:
  """
  echo $workflow.manifest.version > v_pipeline.txt
  echo $workflow.nextflow.version > v_nextflow.txt
  fastqc --version > v_fastqc.txt
  trim_galore --version > v_trim_galore.txt
  echo \$(bwa 2>&1) > v_bwa.txt
  samtools --version > v_samtools.txt
  bedtools --version > v_bedtools.txt
  echo \$(bamtools --version 2>&1) > v_bamtools.txt
  echo \$(plotFingerprint --version 2>&1) > v_deeptools.txt || true
  picard MarkDuplicates --version &> v_picard.txt  || true
  echo \$(R --version 2>&1) > v_R.txt
  python -c "import pysam; print(pysam.__version__)" > v_pysam.txt
  echo \$(macs2 --version 2>&1) > v_macs2.txt
  touch v_homer.txt
  echo \$(featureCounts -v 2>&1) > v_featurecounts.txt
  preseq &> v_preseq.txt
  multiqc --version > v_multiqc.txt
  $baseDir/bin/chipseq/versions.py &> software_versions_mqc.yaml
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
  $baseDir/bin/markdown_to_html.r $output_docs results_description.html
  """
}

process fmt {
  label 'low_memory'
  tag "${params.name}"

  input:
  path("01.tsv")

  output:
  path "design.tsv", emit: design
  path "control_mapping.tsv", emit: mapping

  script:
  """
  $baseDir/bin/nf.fmt.design.py chipseq 01.tsv design.tsv control_mapping.tsv
  """
}

process fbam {
  label 'process_medium'
  tag "$id"
  publishDir "${params.outdir}/21_bam_clean", mode:'link', overwrite:'true',
    saveAs: { fn ->
      if (fn.indexOf("stat") > 0) "stats/$fn"
      else params.saveBAM ? fn : null
    }

  input:
  tuple val(id), path(bam), path(bai)
  path genome_bed
  path bamtools_filter_config

  output:
  tuple val(id), path("${pre}.bam"), path("${pre}.bam.bai"), emit: bam
  tuple val(id), path("${pre1}.bam"), emit: nbam
  tuple val(id), path("${pre}.flagstat"), emit: flagstat
  tuple val(id), path("${pre}.{idxstats,stats}"), emit: stat

  script:
  pre1 = "${id}.mLb.clN"
  pre = "${id}.mLb.flT"
  filter_params = params.single_end ? "-F 0x004" : "-F 0x004 -F 0x0008 -f 0x001"
  dup_params = params.keep_dups ? "" : "-F 0x0400"
  multimap_params = params.keep_multi_map ? "" : "-q 1"
  blacklist_params = params.blacklist ? "-L $bed" : ""
  cmd = "cp ${pre1}.bam ${pre}.bam"
  if (!params.single_end) {
    cmd = """
      bampe_rm_orphan.py ${pre1}.bam ${pre}.bam --only_fr_pairs
      samtools sort -@ $task.cpus -o ${pre}.bam -T $pre ${pre1}.bam
      """.stripIndent()
  }
  name_sort_bam = params.single_end ? "" : "samtools sort -n -@ $task.cpus -o ${pre}.namesorted.bam -T $pre ${pre}.bam"
  //$name_sort_bam
  //| bamtools filter -out ${pre1}.bam -script $bamtools_filter_config
  """
  samtools view \\
      $filter_params \\
      $dup_params \\
      $multimap_params \\
      $blacklist_params \\
      -b ${bam} -O bam -o ${pre1}.bam
  samtools sort -@ $task.cpus -o ${pre}.bam -T $pre ${pre1}.bam
  samtools index ${pre}.bam
  samtools flagstat ${pre}.bam > ${pre}.flagstat
  samtools idxstats ${pre}.bam > ${pre}.idxstats
  samtools stats ${pre}.bam > ${pre}.stats
  """
}

process picard {
  label 'process_medium'
  tag "$id"
  publishDir "${params.outdir}/22_bam_metric", mode:'link', overwrite:'true',
    saveAs: { fn ->
      if (fn.indexOf(".tsv") > 0) "bamstats/$fn"
      else if (fn.indexOf("_metrics") > 0) "picard_metrics/$fn"
      else if (fn.indexOf(".pdf") > 0) "pdf/$fn"
      else null
    }

  when:
  !params.skip_picard_metrics

  input:
  tuple val(id), path(bam), path(bai)
  path genome_fasta

  output:
  path "*_metrics", emit: metric
  tuple val(id), path("${id}.tsv"), emit: stat
  path "*.pdf", emit: pdf

  script:
  prefix = "${id}.mLb.clN"
  def avail_mem = 3
  if (!task.memory) {
    log.info "[Picard CollectMultipleMetrics] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this."
  } else {
    avail_mem = task.memory.toGiga()
  }
  """
  bam.py stat $bam > ${id}.tsv
  picard -Xmx${avail_mem}g CollectMultipleMetrics \\
      INPUT=${bam} \\
      OUTPUT=${id} \\
      REFERENCE_SEQUENCE=${genome_fasta} \\
      VALIDATION_STRINGENCY=LENIENT \\
      TMP_DIR=${TMPDIR}
  """
}

process mgstat {
  label "low_memory"
  tag "${params.name}"
  conda "$NXF_CONDA_CACHEDIR/r"

  input:
  path bamstats

  output:
  path "bamstats.tsv"

  script:
  """
  $baseDir/bin/merge.stats.R --opt bam_stat -o bamstats.tsv $bamstats
  """
}

process bigwig {
  label 'process_medium'
  tag "$id"
  publishDir "${params.outdir}/23_bigwig", mode:'link', overwrite:'true',
    saveAs: { fn ->
      if (fn.indexOf("scale_factor.txt") > 0) "scale/$fn"
      else if (fn.indexOf(".bigWig") > 0) "$fn"
      else null
    }

  input:
  tuple val(id), path(bam), path(bai), path(flagstat)
  path genome_sizes

  output:
  tuple val(id), path("${prefix}.bigWig"), emit: bw
  path "*scale_factor.txt", emit: scale
  path "*igv.txt", emit: igv

  script:
  //prefix = "${id}.mLb.clN"
  prefix = "${id}"
  pe_fragment = params.single_end ? "" : "-pc"
  extend = (params.single_end && params.fragment_size > 0) ? "-fs ${params.fragment_size}" : ''
  """
  SCALE_FACTOR=\$(grep 'mapped (' $flagstat | awk '{print 1000000/\$1}')
  echo \$SCALE_FACTOR > ${prefix}.scale_factor.txt
  genomeCoverageBed -ibam ${bam} -bg -scale \$SCALE_FACTOR $pe_fragment $extend | sort -k1,1 -k2,2n >  ${prefix}.bedGraph

  bedGraphToBigWig ${prefix}.bedGraph ${genome_sizes} ${prefix}.bigWig

  find * -type f -name "*.bigWig" -exec echo -e "bigwig/"{}"\\t0,0,178" \\; > ${prefix}.bigWig.igv.txt
  """
}

process metaplot {
  label 'process_high'
  tag "$id"
  publishDir "${params.outdir}/24_meta_plot", mode:'link', overwrite:'true'

  when:
  !params.skip_meta_plot

  input:
  tuple val(id), path(bw)
  path bed

  output:
  path '*.{gz,pdf}', emit: result
  path '*.plotProfile.tab', emit: mqc

  script:
  """
  computeMatrix scale-regions \\
      --regionsFileName $bed \\
      --scoreFileName $bw \\
      --outFileName ${id}.computeMatrix.mat.gz \\
      --outFileNameMatrix ${id}.computeMatrix.vals.mat.gz \\
      --regionBodyLength 1000 \\
      --beforeRegionStartLength 3000 \\
      --afterRegionStartLength 3000 \\
      --skipZeros \\
      --smartLabels \\
      --numberOfProcessors $task.cpus

  plotProfile --matrixFile ${id}.computeMatrix.mat.gz \\
      --outFileName ${id}.plotProfile.pdf \\
      --outFileNameData ${id}.plotProfile.tab
  """
}

process spp {
  label 'process_medium'
  tag "$id"
  publishDir "${params.outdir}/25_spp", mode:'link', overwrite:'true'

  when:
  !params.skip_spp

  input:
  tuple val(id), path(bam), path(bai)
  path spp_correlation_header
  path spp_nsc_header
  path spp_rsc_header

  output:
  path '*.pdf', emit: plot
  path '*.spp.out', emit: data
  path '*_mqc.tsv', emit: mqc

  script:
  """
  RUN_SPP=`which run_spp.R`
  Rscript -e "library(caTools); source(\\"\$RUN_SPP\\")" -c="${bam}" -savp="${id}.spp.pdf" -savd="${id}.spp.Rdata" -out="${id}.spp.out" -p=$task.cpus
  cp $spp_correlation_header ${id}_spp_correlation_mqc.tsv
  Rscript -e "load('${id}.spp.Rdata'); write.table(crosscorr\\\$cross.correlation, file=\\"${id}_spp_correlation_mqc.tsv\\", sep=",", quote=FALSE, row.names=FALSE, col.names=FALSE,append=TRUE)"

  awk -v OFS='\t' '{print "${id}", \$9}' ${id}.spp.out | cat $spp_nsc_header - > ${id}_spp_nsc_mqc.tsv
  awk -v OFS='\t' '{print "${id}", \$10}' ${id}.spp.out | cat $spp_rsc_header - > ${id}_spp_rsc_mqc.tsv
  """
}

process finger {
  label 'process_high'
  publishDir "${params.outdir}/26_fingerprint", mode:'link', overwrite:'true'
  tag "${ip}"

  when:
  !params.skip_plot_fingerprint

  input:
  tuple val(antibody), val(replicatesExist), val(multipleGroups), val(ip), path(ipbam), path(ipbai), control, path(controlbam), path(controlbai), path(ipflagstat)

  output:
  path '*.{txt,pdf}', emit: result
  path '*.raw.txt', emit: mqc

  script:
  extend = (params.single_end && params.fragment_size > 0) ? "--extendReads ${params.fragment_size}" : ''
  """
  plotFingerprint \\
      --bamfiles ${ipbam} ${controlbam} \\
      --plotFile ${ip}.plotFingerprint.pdf \\
      $extend \\
      --labels $ip $control \\
      --outRawCounts ${ip}.plotFingerprint.raw.txt \\
      --outQualityMetrics ${ip}.plotFingerprint.qcmetrics.txt \\
      --skipZeros \\
      --JSDsample ${controlbam} \\
      --numberOfProcessors $task.cpus \\
      --numberOfSamples $params.fingerprint_bins
  """
}

process macs2 {
  label 'process_medium'
  tag "${ip}"
  publishDir "${params.outdir}/27_macs2", mode:'link', overwrite:'true',
    saveAs: { fn ->
      if (fn.indexOf(".tsv") > 0) "qc/$fn"
      else if (fn.indexOf(".igv.txt") > 0) null
      else fn
    }

  when:
  params.macs_gsize

  input:
  tuple val(antibody), val(replicatesExist), val(multipleGroups), val(ip), path(ipbam), path(ipbai), control, path(controlbam), path(controlbai), path(ipflagstat)
  path peak_count_header
  path frip_score_header

  output:
  tuple antibody, replicatesExist, multipleGroups, ip, control, path("${ip}.bed"), emit: peak
  //tuple ip, path("*.{bed,xls,gappedPeak,bdg}"), emit: peak
  tuple ip, path("${ip}.gapped.bed"), emit: peak_gapped optional true
  tuple ip, path("${ip}.xls"), emit: peak_xls
  path "*igv.txt", emit: igv
  path "*_mqc.tsv", emit: mqc

  script:
  broad = params.narrow_peak ? '' : "--broad --broad-cutoff ${params.broad_cutoff}"
  format = params.single_end ? "BAM" : "BAMPE"
  pileup = params.save_macs_pileup ? "-B --SPMR" : ""
  """
  MAPPED=\$(grep 'mapped (' $ipflagstat | awk -v OFS='\t' '{print \$1}')
  if [ \$MAPPED -lt ${params.min_mapped_reads} ]; then extra="--nomodel"; else extra=""; fi
  macs2 callpeak \\
      -t ${ipbam[0]} \\
      -c ${controlbam[0]} \\
      $broad \$extra \\
      -f $format \\
      -g $params.macs_gsize \\
      -n $ip \\
      $pileup \\
      --keep-dup all

  mv ${ip}_peaks.${params.peak_type} ${ip}.bed
  [ -f ${ip}_peaks.gappedPeak ] && mv ${ip}_peaks.gappedPeak ${ip}.gapped.bed
  mv ${ip}_peaks.xls ${ip}.xls

  cat ${ip}.bed | wc -l | awk -v OFS='\t' '{ print "${ip}", \$1 }' | cat $peak_count_header - > ${ip}_peaks.count_mqc.tsv

  READS_IN_PEAKS=\$(intersectBed -a ${ipbam[0]} -b ${ip}.bed -bed -c -f 0.20 | awk -F '\t' '{sum += \$NF} END {print sum}')
  grep 'mapped (' $ipflagstat | awk -v a="\$READS_IN_PEAKS" -v OFS='\t' '{print "${ip}", a/\$1}' | cat $frip_score_header - > ${ip}_peaks.FRiP_mqc.tsv

  find * -type f -name "${ip}.bed" -exec echo -e "peaks/"{}"\\t0,0,178" \\; > ${ip}_peaks.igv.txt
  """
}

process homer {
  label 'process_medium'
  tag "${ip}"
  publishDir "${params.outdir}/27_macs2", mode:'link', overwrite:'true'

  when:
  params.macs_gsize

  input:
  tuple val(antibody), val(replicatesExist), val(multipleGroups), val(ip), val(control), path(peak)
  path fasta
  path gtf

  output:
  path "*.txt"

  script:
  """
  annotatePeaks.pl \\
      $peak \\
      $fasta \\
      -gid \\
      -gtf $gtf \\
      -cpu $task.cpus \\
      > ${ip}.txt
  """
}

process pqc {
  label "process_medium"
  tag "${params.name}"
  publishDir "${params.outdir}/27_macs2/qc", mode:'link', overwrite:'true'

  when:
  params.macs_gsize

  input:
  path peaks
  path annos
  path peak_annotation_header

  output:
  path "*.{txt,pdf}", emit: result
  path "*.tsv", emit: mqc

  script:  // This script is bundled with the pipeline, in nf-core/chipseq/bin/
  """
  $baseDir/bin/chipseq/plot_macs_qc.r \\
      -i ${peaks.join(',')} \\
      -s ${peaks.join(',').replaceAll("_peaks.${params.peak_type}","")} \\
      -o ./ \\
      -p macs_peak

  $baseDir/bin/chipseq/plot_homer_annotatepeaks.r \\
      -i ${annos.join(',')} \\
      -s ${annos.join(',').replaceAll("_peaks.annotatePeaks.txt","")} \\
      -o ./ \\
      -p macs_annotatePeaks

  cat $peak_annotation_header macs_annotatePeaks.summary.txt > macs_annotatePeaks.summary_mqc.tsv
  """
}

process consensus {
  label 'process_long'
  tag "${antibody}"
  publishDir "${params.outdir}/28_consensus/${antibody}", mode:'link', overwrite: 'true',
    saveAs: { fn -> if (fn.endsWith(".igv.txt")) null else fn }

  when:
  params.macs_gsize && (replicatesExist || multipleGroups) && peaks.collect().size > 1

  input:
  tuple val(antibody), val(replicatesExist), val(multipleGroups), path(peaks)

  output:
  tuple val(antibody), val(replicatesExist), val(multipleGroups), path("*.bed"), emit: bed
  tuple val(antibody), file("*.saf"), emit: saf
  path "*.boolean.txt", emit: bool
  path "*.intersect.{txt,plot.pdf}", emit: ovlp
  path "*igv.txt", emit: igv

  script: // scripts are bundled with the pipeline, in nf-core/chipseq/bin/
  //prefix = "${antibody}.consensus_peaks"
  prefix = "${antibody}"
  mergecols = params.narrow_peak ? (2..10).join(',') : (2..9).join(',')
  collapsecols = params.narrow_peak ? (["collapse"]*9).join(',') : (["collapse"]*8).join(',')
  expandparam = params.narrow_peak ? "--is_narrow_peak" : ""
  """
  sort -k1,1 -k2,2n ${peaks.collect{it.toString()}.sort().join(' ')} \\
      | mergeBed -c $mergecols -o $collapsecols > ${prefix}.txt

  $baseDir/bin/chipseq/macs2_merged_expand.py ${prefix}.txt \\
      ${peaks.collect{it.toString()}.sort().join(',').replaceAll("_peaks.${params.peak_type}","")} \\
      ${prefix}.boolean.txt \\
      --min_replicates $params.min_reps_consensus \\
      $expandparam

  awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print \$1, \$2, \$3, \$4, "0", "+" }' ${prefix}.boolean.txt > ${prefix}.bed

  echo -e "GeneID\tChr\tStart\tEnd\tStrand" > ${prefix}.saf
  awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print \$4, \$1, \$2, \$3,  "+" }' ${prefix}.boolean.txt >> ${prefix}.saf

  $baseDir/bin/chipseq/plot_peak_intersect.r -i ${prefix}.boolean.intersect.txt -o ${prefix}.boolean.intersect.plot.pdf

  find * -type f -name "${prefix}.bed" -exec echo -e "peaks/"{}"\\t0,0,0" \\; > ${prefix}.bed.igv.txt
  """
}

process homer2 {
  label 'process_medium'
  tag "${antibody}"
  publishDir "${params.outdir}/28_consensus/${antibody}", mode:'link', overwrite: 'true'

  when:
  params.macs_gsize && (replicatesExist || multipleGroups)

  input:
  tuple val(antibody), val(replicatesExist), val(multipleGroups), path(bed)
  path bool
  path fasta
  path gtf

  output:
  path "${antibody}*.txt"

  script:
  prefix = "${antibody}"
  """
  annotatePeaks.pl \\
      $bed \\
      $fasta \\
      -gid \\
      -gtf $gtf \\
      -cpu $task.cpus \\
      > ${prefix}.txt

  cut -f2- ${prefix}.txt | awk 'NR==1; NR > 1 {print \$0 | "sort -k1,1 -k2,2n"}' | cut -f6- > tmp.txt
  paste $bool tmp.txt > ${prefix}.boolean.txt
  """
}

process deseq2 {
  label 'process_medium'
  tag "${antibody}"
  publishDir "${params.outdir}/28_consensus/${antibody}/deseq2", mode:'link', overwrite: 'true',
    saveAs: { fn -> if (fn.endsWith(".igv.txt")) null else fn }

  when:
  params.macs_gsize && replicatesExist && multipleGroups && !params.skip_diff_analysis

  input:
  tuple val(antibody), val(replicatesExist), val(multipleGroups), path(bams), path(saf)
  path deseq2_pca_header
  path deseq2_clustering_header

  output:
  path "*featureCounts.txt", emit: cnt
  path "*featureCounts.txt.summary", emit: cnt_mqc
  path "*.{RData,results.txt,pdf,log}", emit: result
  path "sizeFactors", emit: factor
  path "*vs*/*.{pdf,txt}", emit: comp_result
  path "*vs*/*.bed", emit: comp_bed
  path "*igv.txt", emit: comp_igv
  path "*.tsv", emit: mqc

  script:
  prefix = "${antibody}.consensus_peaks"
  bam_files = bams.findAll { it.toString().endsWith('.bam') }.sort()
  bam_ext = params.single_end ? ".mLb.clN.sorted.bam" : ".mLb.clN.bam"
  pe_params = params.single_end ? '' : "-p --donotsort"
  """
  featureCounts \\
      -F SAF \\
      -O \\
      --fracOverlap 0.2 \\
      -T $task.cpus \\
      $pe_params \\
      -a $saf \\
      -o ${prefix}.featureCounts.txt \\
      ${bam_files.join(' ')}

  $baseDir/bin/chipseq/featurecounts_deseq2.r -i ${prefix}.featureCounts.txt -b '$bam_ext' -o ./ -p $prefix -s .mLb

  sed 's/deseq2_pca/deseq2_pca_${task.index}/g' <$deseq2_pca_header >tmp.txt
  sed -i -e 's/DESeq2:/${antibody} DESeq2:/g' tmp.txt
  cat tmp.txt ${prefix}.pca.vals.txt > ${prefix}.pca.vals_mqc.tsv

  sed 's/deseq2_clustering/deseq2_clustering_${task.index}/g' <$deseq2_clustering_header >tmp.txt
  sed -i -e 's/DESeq2:/${antibody} DESeq2:/g' tmp.txt
  cat tmp.txt ${prefix}.sample.dists.txt > ${prefix}.sample.dists_mqc.tsv

  find * -type f -name "*.FDR0.05.results.bed" -exec echo -e "bwa/mergedLibrary/macs/${params.peak_type}/consensus/${antibody}/deseq2/"{}"\\t255,0,0" \\; > ${prefix}.igv.txt
  """
}

process igv {
  tag "${params.name}"
  publishDir "${params.outdir}/29_igv", mode:'link', overwrite:'true'

  when:
  !params.skip_igv

  input:
  path fasta
  path bigwigs
  path peaks
  path consensus_peaks
  path differential_peaks

  output:
  path "*.{txt,xml}"

  script: // scripts are bundled with the pipeline, in nf-core/chipseq/bin/
  """
  cat *.txt > igv_files.txt
  $baseDir/bin/igv_files_to_session.py igv_session.xml igv_files.txt ../../reference_genome/${fasta.getName()} --path_prefix './'
  """
}

workflow chipseq {
  take:
    ibams
    reads
    mapping
  main:
    // stage channels
      def PEAK_TYPE = params.narrow_peak ? "narrowPeak" : "broadPeak"
      genome_bed = file(params.genome_bed, checkIfExists: true)
      genome_sizes = file(params.genome_sizes, checkIfExists: true)
      genome_fasta = file(params.fasta, checkIfExists: true)
      gene_bed = file(params.bed12, checkIfExists: true)
      gene_gtf = file(params.gtf, checkIfExists: true)
      cfg_fbam = params.single_end ?
        file(params.bamtools_filter_se_config, checkIfExists: true) :
        file(params.bamtools_filter_pe_config, checkIfExists: true)
      spp_corr = file("$baseDir/assets/multiqc/spp_correlation_header.txt", checkIfExists: true)
      spp_nsc = file("$baseDir/assets/multiqc/spp_nsc_header.txt", checkIfExists: true)
      spp_rsc = file("$baseDir/assets/multiqc/spp_rsc_header.txt", checkIfExists: true)
      peak_count = file("$baseDir/assets/multiqc/peak_count_header.txt", checkIfExists: true)
      frip_score = file("$baseDir/assets/multiqc/frip_score_header.txt", checkIfExists: true)
      peak_anno = file("$baseDir/assets/multiqc/peak_annotation_header.txt", checkIfExists: true)
      deseq2_pca = file("$baseDir/assets/multiqc/deseq2_pca_header.txt", checkIfExists: true)
      deseq2_clustering = file("$baseDir/assets/multiqc/deseq2_clustering_header.txt", checkIfExists: true)

    fbam(ibams, genome_bed, cfg_fbam)
    picard(fbam.out.bam, genome_fasta)
    mgstat(picard.out.stat.collect{ it[1] }.ifEmpty([]))
    good_samples = mgstat.out
      .splitCsv(header:true, sep:"\t")
      .filter { it.pair_map.toInteger() + it.unpair_map.toInteger() >= params.min_mapped_reads }
      .map { it -> [it.sid] }
    bams = fbam.out.bam.join(good_samples, by:0)
    bigwig(bams.join(fbam.out.flagstat, by: [0]), genome_sizes)
    metaplot(bigwig.out.bw, gene_bed)
    spp(bams, spp_corr, spp_nsc, spp_rsc)

    bam_by_bam = bams.combine(bams)
    mapping2 = mapping.combine(bam_by_bam)
      .filter { it[0] == it[5] && it[1] == it[8] }
      .join(fbam.out.flagstat, by:0)
      .map { it ->  it[2..-1] }
    mapping2 | finger
    macs2(mapping2, peak_count, frip_score)
    peaks = macs2.out.peak
    homer(peaks, genome_fasta, gene_gtf)
    pqc(peaks.collect{ it[-1] }, homer.out.collect(), peak_anno)

    macs_consensus = peaks
        .map { it ->  [ it[0], it[1], it[2], it[-1] ] }
        .groupTuple()
        .map { it ->  [ it[0], it[1][0], it[2][0], it[3].sort() ] }
    consensus(macs_consensus)
    homer2(consensus.out.bed, consensus.out.bool, genome_fasta, gene_gtf)
    bam_deseq = mapping2
      .map { it -> [ it[3], [ it[0], it[1], it[2] ] ] }
      .join(fbam.out.nbam)
      .map { it -> [ it[1][0], it[1][1], it[1][2], it[2] ] }
      .groupTuple()
      .map { it -> [ it[0], it[1][0], it[2][0], it[3].flatten().sort() ] }
      .join(consensus.out.saf)
    deseq2(bam_deseq, deseq2_pca, deseq2_clustering)
    igv(genome_fasta, macs2.out.igv.collect().ifEmpty([]),
        consensus.out.igv.collect().ifEmpty([]),
        deseq2.out.comp_igv.collect().ifEmpty([]),
        bigwig.out.igv.collect().ifEmpty([]))
  emit:
    fbam_bam = fbam.out.bam
    fbam_nbam = fbam.out.nbam
    fbam_flagstat = fbam.out.flagstat
    fbam_stat = fbam.out.stat
    picard_metric = picard.out.metric
    picard_stat = picard.out.stat
    picard_pdf = picard.out.pdf
    bigwig_bw = bigwig.out.bw
    bigwig_igv = bigwig.out.igv
    metaplot_result = metaplot.out.result
    metaplot_mqc = metaplot.out.mqc
    spp_plot = spp.out.plot
    spp_data = spp.out.data
    spp_mqc = spp.out.mqc
    finger_result = finger.out.result
    finger_mqc = finger.out.mqc
    macs2_peak = macs2.out.peak
    macs2_peak_gapped = macs2.out.peak_gapped
    macs2_peak_xls = macs2.out.peak_xls
    macs2_igv = macs2.out.igv
    macs2_mqc = macs2.out.mqc
    homer = homer.out
    pqc_result = pqc.out.result
    pqc_mqc = pqc.out.mqc
    consensus_bed = consensus.out.bed
    consensus_saf = consensus.out.saf
    consensus_igv = consensus.out.igv
    consensus_ovlp = consensus.out.ovlp
    homer2 = homer2.out
    deseq2_cnt = deseq2.out.cnt
    deseq2_cnt_mqc = deseq2.out.cnt_mqc
    deseq2_comp_igv = deseq2.out.comp_igv
    deseq2_mqc = deseq2.out.mqc
    igv = igv.out
}

process mg {
  label "mid_memory"
  tag "${params.name}"
  conda "$NXF_CONDA_CACHEDIR/r"
  publishDir "${params.outdir}/50_final", mode:'link', overwrite:'true'

  input:
  path meta
  path bamstats
  path macs2
  path consensus

  output:
  path("{bamstats,macs2_count,macs2_frip}.tsv")

  script:
  yid = params.name
  cmds = []
  opts = []
  macs2_cnt = macs2.findAll {it.toString().endsWith(".count_mqc.tsv")}.join(' ')
  macs2_frip = macs2.findAll {it.toString().endsWith(".FRiP_mqc.tsv")}.join(' ')
  """
  $baseDir/bin/merge.stats.R --opt bam_stat -o bamstats.tsv $bamstats
  grep -v '^#' -h ${macs2_cnt} > macs2_count.tsv
  grep -v '^#' -h ${macs2_frip} > macs2_frip.tsv
  """
}

process mqc {
  label 'low_memory'
  tag "${params.name}"
  publishDir "${params.outdir}/40_multiqc", mode:'link', overwrite:'true'

  when:
  !params.skip_multiqc

  input:
  path multiqc_config

  path ('software_versions/*')
  path ('workflow_summary/*')

  path ('fastqc/*')
  path ('trimgalore/*')
  //path ('trimgalore/fastqc/*')

  path ('alignment/mergedLibrary/*')
  path ('alignment/mergedLibrary/*')
  path ('alignment/mergedLibrary/picard_metrics/*')

  path ('macs/*')
  path ('macs/*')
  path ('macs/consensus/*')
  path ('macs/consensus/*')

  path ('preseq/*')
  path ('deeptools/*')
  path ('deeptools/*')
  path ('phantompeakqualtools/*')
  path ('phantompeakqualtools/*')

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
      -m custom_content -m fastqc -m cutadapt -m samtools -m picard -m preseq -m featureCounts -m deeptools -m phantompeakqualtools
  """
}


