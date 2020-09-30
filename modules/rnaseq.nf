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

process outdoc {
  publishDir "${params.outdir}/pipeline_info", mode:'copy', overwrite:'true'

  input:
  path output_docs

  output:
  path "results_description.html"

  script:
  """
  $baseDir/bin/markdown_to_html.r $output_docs results_description.html
  """
}

process rseqc {
  label 'mid_memory'
  tag "$id"
  publishDir "${params.outdir}/26_rseqc", mode:'copy', overwrite:'true',
    saveAs: { fn ->
           if (fn.indexOf("bam_stat.txt") > 0)                      "bam_stat/$fn"
      else if (fn.indexOf("infer_experiment.txt") > 0)              "infer_experiment/$fn"
      else if (fn.indexOf("read_distribution.txt") > 0)             "read_distribution/$fn"
      else if (fn.indexOf("read_duplication.DupRate_plot.pdf") > 0) "read_duplication/$fn"
      else if (fn.indexOf("read_duplication.DupRate_plot.r") > 0)   "read_duplication/rscripts/$fn"
      else if (fn.indexOf("read_duplication.pos.DupRate.xls") > 0)  "read_duplication/dup_pos/$fn"
      else if (fn.indexOf("read_duplication.seq.DupRate.xls") > 0)  "read_duplication/dup_seq/$fn"
      else if (fn.indexOf("RPKM_saturation.eRPKM.xls") > 0)         "RPKM_saturation/rpkm/$fn"
      else if (fn.indexOf("RPKM_saturation.rawCount.xls") > 0)      "RPKM_saturation/counts/$fn"
      else if (fn.indexOf("RPKM_saturation.saturation.pdf") > 0)    "RPKM_saturation/$fn"
      else if (fn.indexOf("RPKM_saturation.saturation.r") > 0)      "RPKM_saturation/rscripts/$fn"
      else if (fn.indexOf("inner_distance.txt") > 0)                "inner_distance/$fn"
      else if (fn.indexOf("inner_distance_freq.txt") > 0)           "inner_distance/data/$fn"
      else if (fn.indexOf("inner_distance_plot.r") > 0)             "inner_distance/rscripts/$fn"
      else if (fn.indexOf("inner_distance_plot.pdf") > 0)           "inner_distance/plots/$fn"
      else if (fn.indexOf("junction_plot.r") > 0)                   "junction_annotation/rscripts/$fn"
      else if (fn.indexOf("junction.xls") > 0)                      "junction_annotation/data/$fn"
      else if (fn.indexOf("splice_events.pdf") > 0)                 "junction_annotation/events/$fn"
      else if (fn.indexOf("splice_junction.pdf") > 0)               "junction_annotation/junctions/$fn"
      else if (fn.indexOf("junctionSaturation_plot.pdf") > 0)       "junction_saturation/$fn"
      else if (fn.indexOf("junctionSaturation_plot.r") > 0)         "junction_saturation/rscripts/$fn"
      else fn
    }

  when:
  !params.skip_qc && !params.skip_rseqc

  input:
  tuple val(id), path(bam), path(bai), path(bed)

  output:
  path("*.{txt,pdf,r,xls}")

  script:
  """
  infer_experiment.py -i $bam -r $bed > ${id}.infer_experiment.txt
  junction_annotation.py -i $bam -o ${id}.rseqc -r $bed 2> ${id}.junction_annotation_log.txt
  bam_stat.py -i $bam 2> ${id}.bam_stat.txt
  junction_saturation.py -i $bam -o ${id}.rseqc -r $bed
  inner_distance.py -i $bam -o ${id}.rseqc -r $bed
  read_distribution.py -i $bam -r $bed > ${id}.read_distribution.txt
  read_duplication.py -i $bam -o ${id}.read_duplication
  """
}

process qmap {
  label 'low_memory'
  tag "$id"
  publishDir "${params.outdir}/26_qualimap", mode:'copy', overwrite:'true'

  when:
  !params.skip_qc && !params.skip_qualimap

  input:
  tuple val(id), path(bam), path(bai), path(gtf), val(paired), path(reads)

  output:
  path "${id}"

  script:
  def drc = params.stranded ? params.stranded=='reverse' ? 'strand-specific-reverse' : 'strand-specific-forward' : 'non-strand-specific'
  def pair_str = paired ? '-pe' : ''
  memory = task.memory.toGiga() + "G"
  """
  unset DISPLAY
  qualimap --java-mem-size=${memory} rnaseq -p $drc $pair_str \\
    -bam $bam -gtf $gtf -outdir ${id}
  """
}

process duprad {
  label 'low_memory'
  tag "$id"
  publishDir "${params.outdir}/26_dupradar", mode:'copy', overwrite:'true',
    saveAs: { fn ->
      if (fn.indexOf("_duprateExpDens.pdf") > 0) "scatter_plots/$fn"
      else if (fn.indexOf("_duprateExpBoxplot.pdf") > 0) "box_plots/$fn"
      else if (fn.indexOf("_expressionHist.pdf") > 0) "histograms/$fn"
      else if (fn.indexOf("_dupMatrix.txt") > 0) "gene_data/$fn"
      else if (fn.indexOf("_duprateExpDensCurve.txt") > 0) "scatter_curve_data/$fn"
      else if (fn.indexOf("_intercept_slope.txt") > 0) "intercepts_slopes/$fn"
      else fn
    }

  when:
  !params.skip_qc && !params.skip_dupradar

  input:
  tuple val(id), path(bam), path(bai), path(gtf), val(paired), path(reads)

  output:
  path "*.{pdf,txt}"

  script:
  def drc = params.stranded ? params.stranded=='reverse' ? 2:1:0
  def pair_str = paired ? 'paired':'single'
  """
  $baseDir/bin/rnaseq/dupRadar.r $bam $gtf $drc $pair_str ${task.cpus}
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
  tuple val(id), path(bam), path(bai), path(gtf), path(biotypes_header)

  output:
  path "${id}_gene.featureCounts.txt", emit: txt
  path "${id}_gene.featureCounts.txt.summary", emit: log
  path "${id}_biotype_counts*mqc.{txt,tsv}", emit: biotype optional true

  script:
  def extra = params.fc_extra_attributes ? "--extraAttributes ${params.fc_extra_attributes}" : ''
  def drc = params.stranded ? params.stranded=='reverse' ? 2 : 1 : 0
  biotype_qc = params.skip_biotype_qc ? '' : """
    featureCounts -a $gtf -g ${params.biotype} \\
    -o ${id}_biotype.featureCounts.txt -p \\
    -s $drc $bam
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
  featureCounts -a $gtf -Q ${params.mapQuality} -T ${task.cpus} \\
    -g ${params.fc_group_features} -t ${params.fc_count_type} \\
    -o ${id}_gene.featureCounts.txt \\
    $extra -p -s $drc $bam
  """
}

process corr {
  label 'low_memory'
  tag "${params.name}"
  publishDir "${params.outdir}/33_sample_correlation", mode:'copy', overwrite:'true'

  when:
  !params.skip_qc && !params.skip_edger

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
  $baseDir/bin/rnaseq/edgeR_heatmap_MDS.r $fcnts
  cat $mdsplot_header edgeR_MDS_Aplot_coordinates_mqc.csv >> tmp_file
  mv tmp_file edgeR_MDS_Aplot_coordinates_mqc.csv
  cat $heatmap_header log2CPM_sample_correlation_mqc.csv >> tmp_file
  mv tmp_file log2CPM_sample_correlation_mqc.csv
  """
}

process stie {
  label "mid_memory"
  tag "$name"
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
  tuple val(name), path(bam), path(bai), path(gtf)

  output:
  path "${name}_transcripts.gtf", emit: gtf
  path "${name}.gene_abund.txt", emit: abund
  path "${name}.cov_refs.gtf", emit: cov
  path "${name}_ballgown", emit: ballgown

  script:
  def st_direction = params.stranded ? params.stranded=='reverse' ? "--rf" : "--fr" : ''
  def ignore_gtf = params.stringtie_ignore_gtf ? "" : "-e"
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

process salm {
  label 'mid_memory'
  tag "$name"
  publishDir "${params.outdir}/34_salmon", mode:'copy', overwrite:'true'

  when:
  params.run_salmon

  input:
  tuple val(name), val(paired), path(reads), path(index), path(gtf), path("tx2gene.csv")

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
  unmapped = params.save_unmapped ? "--writeUnmappedNames" : ''
  """
  salmon quant --validateMappings \\
               --seqBias --useVBOpt --gcBias \\
               --geneMap ${gtf} \\
               --threads ${task.cpus} \\
               --libType=${rnastrandness} \\
               --index ${index} \\
               $input $unmapped\\
               -o ${name}
  $baseDir/bin/rnaseq/tximport.r NULL ${name} ${name}
  """
}

process vnt1 {
  label "process_medium"
  tag "${params.name}.$rid"

  input:
  tuple path(bams), path(bais), path(ref), rid, region

  output:
  tuple rid, path("${rid}.vcf.gz"), path("${rid}.vcf.gz.tbi")

  when:
  params.ril

  script:
  """
  bcftools mpileup -f $ref -r $region -Ou $bams |\
    bcftools call -c -Oz -o ${rid}.vcf.gz
  bcftools index -t ${rid}.vcf.gz
  """
}

process vnt2 {
  label "process_medium"
  tag "${params.name}"

  input:
  val rids
  path(vcfs)
  path(tbis)

  output:
  tuple path("vnt.vcf.gz"), path("vnt.vcf.gz.tbi")

  when:
  params.ril

  script:
  vcf_str = rids.sort().collect {"${it}.vcf.gz"}.join(' ')
  """
  bcftools concat -n $vcf_str -Oz -o vnt.vcf.gz
  bcftools index -t vnt.vcf.gz
  """
}

include {get_ch_bcf} from "./utils.nf"
include {ase1; ase2} from './ase.nf'
include {ril1; ril2; ril3} from './ril.nf'
include {cage1; cage2} from './cage.nf'

workflow rnaseq {
  take:
    bams
    reads
    design
  main:
    gtf = Channel.fromPath(params.gtf, checkIfExists: true)
       .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
    bed = Channel
      .fromPath(params.bed, checkIfExists: true)
      .ifEmpty { exit 1, "BED annotation file not found: ${params.bed}" }
    pbed = Channel
      .fromPath(params.pbed, checkIfExists: true)
      .ifEmpty { exit 1, "BED annotation file not found: ${params.pbed}" }
    rseqc(bams.combine(bed))
    qmap(bams.combine(gtf).join(reads))
    duprad(bams.combine(gtf).join(reads))

    ch_biotypes_header = Channel.fromPath("$baseDir/assets/biotypes_header.txt", checkIfExists: true)
    ch_mdsplot_header = Channel.fromPath("$baseDir/assets/mdsplot_header.txt", checkIfExists: true)
    ch_heatmap_header = Channel.fromPath("$baseDir/assets/heatmap_header.txt", checkIfExists: true)
    fcnt(bams.combine(gtf).combine(ch_biotypes_header))
    corr(fcnt.out.txt.collect(), fcnt.out.txt.count(), ch_mdsplot_header, ch_heatmap_header)

    stie_out = Channel.empty()
    if( params.run_stringtie ) {
      stie(bams.combine(gtf))
      stie_out = stie.out
    }
    salm_gcnt = Channel.empty(); salm_gtpm = Channel.empty()
    salm_tcnt = Channel.empty(); salm_ttpm = Channel.empty()
    if( params.run_salmon ) {
      salmon_index = Channel
        .fromPath(params.salmon_index, checkIfExists: true)
        .ifEmpty { exit 1, "Salmon index not found: ${params.salmon_index}" }
      fna = Channel
        .fromPath(params.fna, checkIfExists: true)
        .ifEmpty { exit 1, "Transcript fasta file not found: ${params.fna}" }
      tx2gene = Channel
        .fromPath(params.tx2gene, checkIfExists: true)
        .ifEmpty { exit 1, "tx2gene not found: ${params.tx2gene}" }
      salm(reads.combine(salmon_index).combine(gtf).combine(tx2gene))
      salm_gcnt = salm.out.gcnt; salm_gtpm = salm.out.gtpm
      salm_tcnt = salm.out.tcnt; salm_ttpm = salm.out.ttpm
    }

    ase_gene = Channel.empty(); ase_snp = Channel.empty()
    if (params.ase) {
      genome_fasta = Channel.fromPath(params.fasta, checkIfExists: true)
        .ifEmpty { exit 1, "Genome fasta file not found: ${params.fasta}" }
      ch_bcf_ase = get_ch_bcf(design)
      ase1(bams.join(ch_bcf_ase).combine(genome_fasta))
      ase2(ase1.out.combine(gtf))
      ase_gene = ase2.out.gene
      ase_snp = ase2.out.snp
    }

    ril_csv = Channel.empty(); ril_txt = Channel.empty()
    if (params.ril) {
      genome_fasta = Channel.fromPath(params.fasta, checkIfExists: true)
        .ifEmpty { exit 1, "Genome fasta file not found: ${params.fasta}" }
      ril_sites = Channel
        .fromPath(params.ril_sites, checkIfExists: true)
        .ifEmpty { exit 1, "RIL sites file not found: ${params.ril_sites}" }
      ril_sites_idx = Channel
        .fromPath(params.ril_sites_idx, checkIfExists: true)
        .ifEmpty { exit 1, "RIL sites index not found: ${params.ril_sites_idx}" }
      win10 = Channel
        .fromPath(params.win10, checkIfExists: true)
        .ifEmpty { exit 1, "win10 tsv file not found: ${params.win10}" }
        .splitCsv(header:true, sep:"\t")
        .map { row -> [ row.rid, row.region ] }
      win56 = Channel
        .fromPath(params.win56, checkIfExists: true)
        .ifEmpty { exit 1, "win56 tsv file not found: ${params.win56}" }
        .splitCsv(header:true, sep:"\t")
        .map { row -> [ row.rid, row.region ] }

      ibams = bams.collect({it[1]}).toSortedList()
      ibais = bams.collect({it[2]}).toSortedList()
      ril1(ibams.combine(ibais).combine(genome_fasta).combine(ril_sites).combine(ril_sites_idx).combine(win56))
      ril2_in = win56.join(ril1.out, by:0)
        .toSortedList {entry -> entry[0]}
      vcfs = ril2_in.map { it.collect {it[2]} }.collect()
      tbis = ril2_in.map { it.collect {it[3]} }.collect()
      ril2(win56.collect{it[0]}, vcfs, tbis)
      ril3(ril2.out.combine(win10))
      ril_csv = ril3.out.csv; ril_txt = ril3.out.txt
    }

    bigwigs = Channel.empty();
    if (params.cage) {
      cage1(bams)
      bigwigs = cage1.out
      cage2(bigwigs.collect({[it[1],it[2]]}).flatten().toSortedList(),  pbed)
    }

    vnt_vcf = Channel.empty()
    if (params.vntcall) {
      genome_fasta = Channel.fromPath(params.fasta, checkIfExists: true)
        .ifEmpty { exit 1, "Genome fasta file not found: ${params.fasta}" }
      win11 = Channel
        .fromPath(params.win11, checkIfExists: true)
        .ifEmpty { exit 1, "win11 tsv file not found: ${params.win11}" }
        .splitCsv(header:true, sep:"\t")
        .map { row -> [ row.rid, row.region ] }
      win56 = Channel
        .fromPath(params.win56, checkIfExists: true)
        .ifEmpty { exit 1, "win56 tsv file not found: ${params.win56}" }
        .splitCsv(header:true, sep:"\t")
        .map { row -> [ row.rid, row.region ] }

      ibams = bams.collect({it[1]}).toSortedList()
      ibais = bams.collect({it[2]}).toSortedList()
      vnt1(ibams.combine(ibais).combine(genome_fasta).combine(win56))
      vnt2_in = win56.join(vnt1.out, by:0)
        .toSortedList {entry -> entry[0]}
      vcfs = vnt2_in.map { it.collect {it[2]} }.collect()
      tbis = vnt2_in.map { it.collect {it[3]} }.collect()
      vnt2(win56.collect{it[0]}, vcfs, tbis)
      vnt_vcf = vnt2.out
    }

  emit:
    fcnt_txt = fcnt.out.txt
    fcnt_log = fcnt.out.log
    fcnt_biotype = fcnt.out.biotype
    corr = corr.out
    rseqc = rseqc.out
    qmap = qmap.out
    duprad = duprad.out
    stie = stie_out
    salm_gcnt = salm_gcnt
    salm_gtpm = salm_gtpm
    salm_tcnt = salm_tcnt
    salm_ttpm = salm_ttpm
    ase_gene = ase_gene
    ase_snp = ase_snp
    ril_csv = ril_csv
    ril_txt = ril_txt
    vnt_vcf = vnt_vcf
    cage = bigwigs
}

process mg {
  label "mid_memory"
  tag "${params.name}"
  conda "$NXF_CONDA_CACHEDIR/r"
  publishDir "${params.outdir}/50_final", mode:'copy', overwrite:'true'
  publishDir "${params.qcdir}/${params.genome}/${params.name}", mode:'copy', overwrite:'true'

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
  if (params.run_salmon) {
    cmds += "merge.stats.R --opt salmon -o salmon_gene_cnt.rds $salm_gene_cnt"
    cmds += "merge.stats.R --opt salmon -o salmon_gene_tpm.rds $salm_gene_tpm"
    cmds += "merge.stats.R --opt salmon -o salmon_tx_cnt.rds $salm_tx_cnt"
    cmds += "merge.stats.R --opt salmon -o salmon_tx_tpm.rds $salm_tx_tpm"
    opts += "--salmon_gcnt salmon_gene_cnt.rds --salmon_gtpm salmon_gene_tpm.rds"
    opts += "--salmon_tcnt salmon_tx_cnt.rds   --salmon_ttpm salmon_tx_tpm.rds"
  }
  cmd = cmds.join("\n")
  opt = opts.join(" ")
  """
  $baseDir/bin/merge.stats.R --opt bam_stat -o bamstats.tsv $bamstats
  $baseDir/bin/merge.stats.R --opt featurecounts -o featurecounts.rds $fcnts
  $cmd

  $baseDir/bin/rnaseq/nf_rnaseq.R $meta merged.rds --bamstat bamstats.tsv \\
    --fcnt featurecounts.rds \\
    $opt
  $baseDir/bin/rnaseq/nf_rnaseq_norm.R merged.rds 00.raw.rds --rcfg $rcfg
  """
}

process renorm {
  label 'process_medium'
  tag "${params.name}"
  conda "$NXF_CONDA_CACHEDIR/r"
  publishDir "${params.outdir}/50_final", mode:'copy', overwrite:'true'
  publishDir "${params.qcdir}/${params.genome}/${params.name}", mode:'copy', overwrite:'true'
  publishDir "${params.s3dir}/${params.genome.toLowerCase().replaceAll(/_/, "-")}", mode:'copy', overwrite:'true',
    saveAs: { fn ->
      if (fn.indexOf(".rds") > 0) "${params.name}.rds"
      else null
    }

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
  $baseDir/bin/rnaseq/nf_rnaseq_norm.R --meta $meta $raw_rds 01.rds --rcfg $rcfg
  """
}

process mqc {
  label 'low_memory'
  tag "${params.name}"
  publishDir "${params.outdir}/40_multiqc", mode:'copy', overwrite:'true'
  publishDir "${params.s3dir}/${params.genome.toLowerCase().replaceAll(/_/, "-")}", mode:'copy', overwrite:'true',
    saveAs: { fn ->
      if (fn.indexOf(".html") > 0) "$fn"
      else null
    }

  when:
  !params.skip_multiqc

  input:
  path multiqc_config
  path ('fastqc/*')
  path ('trimgalore/*')
  path ('alignment/*')
  path ('preseq/*')
  path ('rseqc/*')
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
  path "*_plots", emit: plot

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





