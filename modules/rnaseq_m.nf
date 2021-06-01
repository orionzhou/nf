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
  tuple val(id), path(bam), path(bai), path(gtf), path(biotypes_header), val(paired), path(reads)

  output:
  path "${id}_gene.featureCounts.txt", emit: txt
  path "${id}_gene.featureCounts.txt.summary", emit: log
  path "${id}_biotype_counts*mqc.{txt,tsv}", emit: biotype optional true

  script:
  def mq = params.mapQuality
  def flag_attr = params.fc_extra_attributes ? "--extraAttributes ${params.fc_extra_attributes}" : ''
  def flag_srd = params.stranded == 'no' ? 0 : params.stranded=='reverse' ? 2 : 1
  def flag_pe = paired == 'PE' ? "-p" : ""
  def flag_long = params.read_type == 'nanopore' ? "-L" : ""
  biotype_qc = params.skip_biotype_qc ? '' : """
    featureCounts -a $gtf -g ${params.biotype} \\
    -o ${id}_biotype.featureCounts.txt -p \\
    -s $flag_srd $bam
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
  featureCounts -a $gtf --primary -Q $mq -T ${task.cpus} \\
    -g ${params.fc_group_features} -t ${params.fc_count_type} \\
    -o ${id}_gene.featureCounts.txt \\
    $flag_attr $flag_pe $flag_long -s $flag_srd \\
    $bam
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
  def st_direction = params.stranded == 'no' ? '' : params.stranded=='reverse' ? "--rf" : "--fr"
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
  def rnastrandness = paired == 'SE' ? 'U' : 'IU'
  if (params.stranded == 'forward') {
      rnastrandness = paired == 'SE' ? 'SF' : 'ISF'
  } else if (params.stranded == 'reverse') {
      rnastrandness = paired == 'SE' ? 'SR' : 'ISR'
  }
  def input = paired=='SE' ? "-r ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
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

process bigwig {
  label 'mid_memory'
  tag "$id"
  publishDir "${params.outdir}/23_bigwig", mode:'copy', overwrite:'true'

  input:
  tuple val(id), path(bam), path(bai)

  output:
  tuple val(id), path("${id}.bw")

  script:
  binsize = 50
  //--effectiveGenomeSize ${egs} 
  """
  bamCoverage --binSize ${binsize} -p ${task.cpus} -b $bam -o ${id}.bw
  """
}

include {cage1; cage2} from './cage.nf'

workflow rnaseq {
  take:
    bams
    reads
    design
  main:
    genome_sizes = Channel.fromPath(params.genome_sizes, checkIfExists: true)
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
    fcnt(bams.combine(gtf).combine(ch_biotypes_header).join(reads))
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
    } else {
      bigwig(bams)
      bigwigs = bigwig.out
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
  //conda "$NXF_CONDA_CACHEDIR/r"
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

  output:
  path "00.raw.rds"

  script:
  yid = params.name
  //rc2cpm.R $meta featurecounts.rds cpm.rds --opt featurecounts --config $cfg
  cmds = []
  opts = []
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
  """
}




