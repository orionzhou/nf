#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include '../modules/utils.nf'
include './helper.nf'

if (params.help) { help(); exit 0 }

prep_params(params, workflow)

// stage channels, variables, files
  nofile = Channel.fromPath("NO_FILE")
  nofile2 = Channel.fromPath("NO_FILE2")
  nofile3 = Channel.fromPath("NO_FILE3")
  nofile4 = Channel.fromPath("NO_FILE4")
  wherefile = Channel.fromPath("$baseDir/assets/where_are_my_files.txt", checkIfExists: true)
  ch_mdsplot_header = Channel.fromPath("$baseDir/assets/mdsplot_header.txt", checkIfExists: true)
  ch_heatmap_header = Channel.fromPath("$baseDir/assets/heatmap_header.txt", checkIfExists: true)
  ch_biotypes_header = Channel.fromPath("$baseDir/assets/biotypes_header.txt", checkIfExists: true)
// Stage config files
  ch_mqc_cfg = file(params.multiqc_config, checkIfExists: true)
  ch_out_doc = file("$baseDir/docs/output.md", checkIfExists: true)
ch_multiqc_config = file("$baseDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_output_docs = file("$baseDir/docs/output.md", checkIfExists: true)

// validate inputs
  nofile = Channel.fromPath("NO_FILE")
  design = Channel
    .fromPath(params.design, checkIfExists: true)
    .ifEmpty { exit 1, "sample design table missing: ${params.design}" }
  design2 = Channel
    .fromPath("${params.qcdir}/${params.name}/01.meta.tsv")
if( params.aligner =~ /bismark/ ){
    assert params.bismark_index || params.fasta : "No reference genome index or fasta file specified"
    ch_wherearemyfiles_for_alignment.set { ch_wherearemyfiles_for_bismark_align }

    if( params.bismark_index ){
        Channel
            .fromPath(params.bismark_index, checkIfExists: true)
            .ifEmpty { exit 1, "Bismark index file not found: ${params.bismark_index}" }
            .into { ch_bismark_index_for_bismark_align; ch_bismark_index_for_bismark_methXtract }
    }
    else if( params.fasta ){
        Channel
            .fromPath(params.fasta, checkIfExists: true)
            .ifEmpty { exit 1, "fasta file not found : ${params.fasta}" }
            .set { ch_fasta_for_makeBismarkIndex }
    }
}
else if( params.aligner == 'bwameth' ){
    assert params.fasta : "No Fasta reference specified! This is required by MethylDackel."
    ch_wherearemyfiles_for_alignment.into { ch_wherearemyfiles_for_bwamem_align; ch_wherearemyfiles_for_samtools_sort_index_flagstat }

    Channel
        .fromPath(params.fasta, checkIfExists: true)
        .ifEmpty { exit 1, "fasta file not found : ${params.fasta}" }
        .into { ch_fasta_for_makeBwaMemIndex; ch_fasta_for_makeFastaIndex; ch_fasta_for_methyldackel }

    if( params.bwa_meth_index ){
        Channel
            .fromPath("${params.bwa_meth_index}*", checkIfExists: true)
            .ifEmpty { exit 1, "bwa-meth index file(s) not found: ${params.bwa_meth_index}" }
            .set { ch_bwa_meth_indices_for_bwamem_align }
        ch_fasta_for_makeBwaMemIndex.close()
    }

    if( params.fasta_index ){
        Channel
            .fromPath(params.fasta_index, checkIfExists: true)
            .ifEmpty { exit 1, "fasta index file not found: ${params.fasta_index}" }
            .set { ch_fasta_index_for_methyldackel }
        ch_fasta_for_makeFastaIndex.close()
    }
}
  if (params.fasta) {
    if (hasExtension(params.fasta, 'gz')) {
      Channel.fromPath(params.fasta, checkIfExists: true)
        .ifEmpty { exit 1, "Genome fasta file not found: ${params.fasta}" }
        .set { genome_fasta_gz }
    } else {
      genome_fasta = Channel.fromPath(params.fasta, checkIfExists: true)
        .ifEmpty { exit 1, "Genome fasta file not found: ${params.fasta}" }
    }
  } else {
    exit 1, "No reference genome files specified!"
  }

  if (params.win11) {
    win11 = Channel
      .fromPath(params.win11, checkIfExists: true)
      .ifEmpty { exit 1, "win11 tsv file not found: ${params.win11}" }
      .splitCsv(header:true, sep:"\t")
      .map { row -> [ row.rid, row.region ] }
  }
  if (params.win56) {
    win56 = Channel
      .fromPath(params.win56, checkIfExists: true)
      .ifEmpty { exit 1, "win56 tsv file not found: ${params.win56}" }
      .splitCsv(header:true, sep:"\t")
      .map { row -> [ row.rid, row.region ] }
  }
  if (params.win127) {
    win127 = Channel
      .fromPath(params.win127, checkIfExists: true)
      .ifEmpty { exit 1, "win127 tsv file not found: ${params.win127}" }
      .splitCsv(header:true, sep:"\t")
      .map { row -> [ row.rid, row.region ] }
  }

include {fqc; trim} from '../modules/fastq.nf'
include {hs2; star} from '../modules/mapping.nf'
include {bamsort; markdup; bamstat; rseqc; pseq; qmap; duprad} from '../modules/bam.nf'
include '../modules/rnaseq.nf'
include {ril1; ril2; ril3} from '../modules/ril.nf'
ch_reads = get_reads(design)
ch_read_num = get_read_num(design)
ch_bcf_ase = get_ch_bcf(design)

def summary = getSummary()
log.info showHeader(summary)
checkHostname()

workflow {
  main:
    version()
    ch_reads | (fqc & trim)
    aln = nofile
    aln_log = nofile
    if (params.aligner == 'hisat2') {
      hs2(trim.out.reads.combine(hs2_indices).combine(splicesites))
      aln = hs2.out.bam
      aln_log = hs2.out.log
    } else if (params.aligner == 'star') {
      star(trim.out.reads.combine(star_index).combine(gtf))
      aln = star.out.bam
      aln_log = star.out.log
    }
    aln | bamsort | (bamstat & markdup)
    pseq(bamsort.out.join(ch_read_num))
    rseqc(bamsort.out.combine(bed12))
    qmap(bamsort.out.combine(gtf).join(ch_reads))
    duprad(bamsort.out.combine(gtf).join(ch_reads))
    fcnt(bamsort.out.combine(gtf).combine(ch_biotypes_header))
    stie(bamsort.out.combine(gtf))
    salm(trim.out.reads.combine(salmon_index).combine(gtf).combine(tx2gene))
    ase_gene = nofile
    ase_snp = nofile2
    if (params.ase) {
      ase1(bamsort.out.join(ch_bcf_ase).combine(genome_fasta))
      ase2(ase1.out.combine(gtf))
      ase_gene = ase2.out.gene
      ase_snp = ase2.out.snp
    }
    ril_csv = nofile3
    ril_txt = nofile4
    if (params.ril) {
      bams = bamsort.out.collect({it[1]}).toSortedList()
      bais = bamsort.out.collect({it[2]}).toSortedList()
      ril1(bams.combine(bais).combine(genome_fasta).combine(ril_sites).combine(ril_sites_idx).combine(win56))
      ril2_in = win56.join(ril1.out, by:0)
        .toSortedList {entry -> entry[0]}
      vcfs = ril2_in.map { it.collect {it[2]} }.collect()
      tbis = ril2_in.map { it.collect {it[3]} }.collect()
      ril2(win56.collect{it[0]}, vcfs, tbis)
      ril3(ril2.out.combine(win11))
      ril_csv = ril3.out.csv
      ril_txt = ril3.out.txt
    }
    corr(fcnt.out.txt.collect(), fcnt.out.txt.count(), ch_mdsplot_header, ch_heatmap_header)
    outdoc(ch_out_doc)
    mqc(ch_mqc_cfg,
      fqc.out.zip.collect().ifEmpty([]),
      trim.out.log.collect().ifEmpty([]),
      aln_log.collect().ifEmpty([]),
      rseqc.out.distr.concat(rseqc.out.infer).concat(rseqc.out.jctsat_r).collect().ifEmpty([]),
      pseq.out.collect().ifEmpty([]),
      qmap.out.collect().ifEmpty([]),
      duprad.out.data.concat(duprad.out.curve).concat(duprad.out.slope).collect().ifEmpty([]),
      fcnt.out.log.collect().ifEmpty([]),
      fcnt.out.biotype.collect().ifEmpty([]),
      corr.out.collect().ifEmpty([]),
      version.out.yml.collect().ifEmpty([]),
      create_workflow_summary(summary)
      )
    mg(design.collect(), rcfg.collect(),
        bamstat.out.collect(), fcnt.out.txt.collect(),
        salm.out.gcnt.collect(), salm.out.gtpm.collect(),
        salm.out.tcnt.collect(), salm.out.ttpm.collect(),
        ase_gene.collect(), ase_snp.collect(),
        ril_csv.collect(), ril_txt.collect())
    renorm(design2.collect(), mg.out, rcfg.collect())
  publish:
    // version
      version.out.csv to: "${params.outdir}/pipeline_info", mode:'link', overwrite:'true'
      version.out.yml to: "${params.outdir}/pipeline_info", mode:'link', overwrite:'true'
      outdoc.out to: "${params.outdir}/pipeline_info", mode:'link', overwrite:'true'
    // fastqc
      fqc.out.html to: "${params.outdir}/01_fastqc", mode:'link', overwrite:'true'
      fqc.out.zip to: "${params.outdir}/01_fastqc/zips", mode:'link', overwrite:'true'
    // trim_galore
      //trim.out.reads to: "${params.outdir}/02_trim_galore", mode:'link', overwrite:'true'
      trim.out.fastqc to: "${params.outdir}/02_trim_galore/fastqc", mode:'link', overwrite:'true'
      trim.out.log to: "${params.outdir}/02_trim_galore/logs", mode:'link', overwrite:'true'
    // hisat2
      //hs2.out.bam to: "${params.outdir}/11_hisat2", mode:'link', overwrite:'true'
      hs2.out.log to: "${params.outdir}/11_hisat2/logs", mode:'link', overwrite:'true'
    // bamsort
    bamsort.out to: "${params.outdir}/20_bam_sorted", mode:'link', overwrite:'true', enabled: params.saveBAM
    bamstat.out to: "${params.outdir}/20_bam_sorted", mode:'link', overwrite:'true'
    // rseqc
      rseqc.out.bamstat to: "${params.outdir}/21_rseqc/bam_stat", mode:'link', overwrite:'true'
      rseqc.out.infer to: "${params.outdir}/21_rseqc/infer_experiment", mode:'link', overwrite:'true'
      rseqc.out.distr to: "${params.outdir}/21_rseqc/read_distribution", mode:'link', overwrite:'true'
      rseqc.out.dup to: "${params.outdir}/21_rseqc/read_duplication", mode:'link', overwrite:'true'
      rseqc.out.dup_r to: "${params.outdir}/21_rseqc/read_duplication/rscripts", mode:'link', overwrite:'true'
      rseqc.out.dup_pos to: "${params.outdir}/21_rseqc/read_duplication/dup_pos", mode:'link', overwrite:'true'
      rseqc.out.dup_seq to: "${params.outdir}/21_rseqc/read_duplication/dup_seq", mode:'link', overwrite:'true'
      rseqc.out.sat to: "${params.outdir}/21_rseqc/RPKM_saturation", mode:'link', overwrite:'true'
      rseqc.out.sat_r to: "${params.outdir}/21_rseqc/RPKM_saturation/rscripts", mode:'link', overwrite:'true'
      rseqc.out.sat_rpkm to: "${params.outdir}/21_rseqc/RPKM_saturation/rpkm", mode:'link', overwrite:'true'
      rseqc.out.sat_cnt to: "${params.outdir}/21_rseqc/RPKM_saturation/counts", mode:'link', overwrite:'true'
      rseqc.out.inndst to: "${params.outdir}/21_rseqc/inner_distance", mode:'link', overwrite:'true'
      rseqc.out.inndst_data to: "${params.outdir}/21_rseqc/inner_distance/data", mode:'link', overwrite:'true'
      rseqc.out.inndst_r to: "${params.outdir}/21_rseqc/inner_distance/rscripts", mode:'link', overwrite:'true'
      rseqc.out.inndst_plot to: "${params.outdir}/21_rseqc/inner_distance/plots", mode:'link', overwrite:'true'
      rseqc.out.jct to: "${params.outdir}/21_rseqc/junction_annotation", mode:'link', overwrite:'true'
      rseqc.out.jct_r to: "${params.outdir}/21_rseqc/junction_annotation/rscripts", mode:'link', overwrite:'true'
      rseqc.out.jct_data to: "${params.outdir}/21_rseqc/junction_annotation/data", mode:'link', overwrite:'true'
      rseqc.out.jct_event to: "${params.outdir}/21_rseqc/junction_annotation/events", mode:'link', overwrite:'true'
      rseqc.out.jct_jct to: "${params.outdir}/21_rseqc/junction_annotation/junctions", mode:'link', overwrite:'true'
      rseqc.out.jctsat to: "${params.outdir}/21_rseqc/junction_saturation", mode:'link', overwrite:'true'
      rseqc.out.jctsat_r to: "${params.outdir}/21_rseqc/junction_saturation/rscripts", mode:'link', overwrite:'true'
    // preseq
    pseq.out to: "${params.outdir}/22_preseq", mode:'link', overwrite:'true'
    // markdup
      //markdup.out.bam to: "${params.outdir}/24_markdup", mode:'link', overwrite:'true'
      markdup.out.metric to: "${params.outdir}/24_markdup/metrics", mode:'link', overwrite:'true'
    // qualimap
    qmap.out to: "${params.outdir}/25_qualimap", mode:'link', overwrite:'true'
    // dupradar
      duprad.out.expDen to: "${params.outdir}/26_dupradar/scatter_plots", mode:'link', overwrite:'true'
      duprad.out.expBox to: "${params.outdir}/26_dupradar/box_plots", mode:'link', overwrite:'true'
      duprad.out.hist to: "${params.outdir}/26_dupradar/histograms", mode:'link', overwrite:'true'
      duprad.out.data to: "${params.outdir}/26_dupradar/gene_data", mode:'link', overwrite:'true'
      duprad.out.curve to: "${params.outdir}/26_dupradar/scatter_curve_data", mode:'link', overwrite:'true'
      duprad.out.slope to: "${params.outdir}/26_dupradar/intercepts_slopes", mode:'link', overwrite:'true'
    // featureCounts
      fcnt.out.txt to: "${params.outdir}/31_featureCounts/gene_counts", mode:'link', overwrite:'true'
      fcnt.out.log to: "${params.outdir}/31_featureCounts/gene_count_summaries", mode:'link', overwrite:'true'
      fcnt.out.biotype to: "${params.outdir}/31_featureCounts/biotype_counts", mode:'link', overwrite:'true'
    // stringtie
      //stie.out.abund to: "${params.outdir}/32_stringtieFPKM", mode:'link', overwrite:'true'
      //stie.out.gtf to: "${params.outdir}/32_stringtieFPKM/transcripts", mode:'link', overwrite:'true'
      //stie.out.cov to: "${params.outdir}/32_stringtieFPKM/cov_refs", mode:'link', overwrite:'true'
      //stie.out.ballgown to: "${params.outdir}/32_stringtieFPKM/ballgown", mode:'link', overwrite:'true'
    corr.out to: "${params.outdir}/33_sample_correlation", mode:'link', overwrite:'true'
    salm.out to: "${params.outdir}/34_salmon", mode:'link', overwrite:'true'
    // ase_gene to: "${params.outdir}/37_ase", mode:'link', overwrite:'true'
    mqc.out to: "${params.outdir}/40_multiqc", mode:'link', overwrite:'true'
    mqc.out.html to: "${params.webdir}", mode:'copy', overwrite:'true'
    mg.out to: "${params.qcdir}/${params.name}", mode:'copy', overwrite:'true'
    renorm.out to: "${params.qcdir}/${params.name}", mode:'copy', overwrite:'true'
}

ch_reads_l = ch_reads.toList()
workflow.onComplete {
  mf = workflow.manifest
  name = params.name ?: mf.name; author = mf.author; url = mf.homePage;
  // Set up the e-mail variables
  def subject = "[$name] Suceeded"
  if (!workflow.success) subject = "[$name] FAILED"
  def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = name
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    if (workflow.container) email_fields['summary']['Docker image'] = workflow.container
    email_fields['percent_aln_skip'] = params.percent_aln_skip
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

  // On success try attach the multiqc report
  def mqc_report = null
  try {
    if (workflow.success && !params.skipMultiQC) {
      mqc_report = mqc.out.html.getVal()
      if (mqc_report.getClass() == ArrayList) {
        log.warn "[$name] Found multiple reports from process 'multiqc', will use only one"
        mqc_report = mqc_report[0]
      }
    }
  } catch (all) {
    log.warn "[$name] Could not attach multiqc report to summary email"
  }

  // Check if we are only sending emails on failure
  email_address = params.email
  if (!params.email && params.email_on_fail && !workflow.success)
    email_address = params.email_on_fail

  // Render the TXT template
  def engine = new groovy.text.GStringTemplateEngine()
  def tf = new File("$baseDir/assets/email_template.txt")
  def txt_template = engine.createTemplate(tf).make(email_fields)
  def email_txt = txt_template.toString()

  // Render the HTML template
  def hf = new File("$baseDir/assets/email_template.html")
  def html_template = engine.createTemplate(hf).make(email_fields)
  def email_html = html_template.toString()

  // Render the sendmail template
  def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.max_multiqc_email_size.toBytes() ]
  def sf = new File("$baseDir/assets/sendmail_template.txt")
  def sendmail_template = engine.createTemplate(sf).make(smail_fields)
  def sendmail_html = sendmail_template.toString()

  // Send the HTML e-mail
  if (email_address) {
    try {
      if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
      // Try to send HTML e-mail using sendmail
      [ 'sendmail', '-t' ].execute() << sendmail_html
      log.info "[$name] Sent summary e-mail to $email_address (sendmail)"
    } catch (all) {
      // Catch failures and try with plaintext
      [ 'mail', '-s', subject, email_address ].execute() << email_txt
      log.info "[$name] Sent summary e-mail to $email_address (mail)"
    }
  }

  // Write summaries to a file
  // Make directory first if needed
  def output_d = new File("${params.outdir}/pipeline_info/")
  if (!output_d.exists()) output_d.mkdirs()
  // Replace the email logo cid with a base64 encoded image
  def logo_b64 = new File("$baseDir/assets/nf-core-rnaseq_logo.png").bytes.encodeBase64().toString()
  email_html = email_html.replace('<img src="cid:nfcorepipelinelogo">', "<img src=\"data:image/png;base64, ${logo_b64}\">")
  // Print to file
  def output_hf = file("${output_d}/pipeline_report.html")
  output_hf.withWriter { w -> w << email_html }
  def output_tf = file("${output_d}/pipeline_report.txt")
  output_tf.withWriter { w -> w << email_txt }

  c_reset = params.monochrome_logs ? '' : "\033[0m";
  c_purple = params.monochrome_logs ? '' : "\033[0;35m";
  c_green = params.monochrome_logs ? '' : "\033[0;32m";
  c_red = params.monochrome_logs ? '' : "\033[0;31m";

  def nsam = ch_reads_l.getVal().size()
  log.info "[${c_purple}${name}${c_reset}] ${c_green}$nsam samples processed${c_reset}"

  if (workflow.stats.ignoredCount > 0 && workflow.success) {
    log.info "- ${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}"
    log.info "- ${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}"
    log.info "- ${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}"
  }

  if (workflow.success) {
    log.info "[${c_purple}${name}${c_reset}] ${c_green}Pipeline completed successfully${c_reset}"
  } else {
    checkHostname()
    log.info "[${c_purple}${name}${c_reset}] ${c_red}Pipeline completed with errors${c_reset}"
  }
}


