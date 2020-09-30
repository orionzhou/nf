def show_header(summary) {
  // Log colors ANSI codes
  c_reset = params.monochrome_logs ? '' : "\033[0m";
  c_dim = params.monochrome_logs ? '' : "\033[2m";
  c_black = params.monochrome_logs ? '' : "\033[0;30m";
  c_green = params.monochrome_logs ? '' : "\033[0;32m";
  c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
  c_blue = params.monochrome_logs ? '' : "\033[0;34m";
  c_purple = params.monochrome_logs ? '' : "\033[0;35m";
  c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
  c_white = params.monochrome_logs ? '' : "\033[0;37m";

  mf = workflow.manifest
  name = mf.name; author = mf.author; url = mf.homePage;
  desc = mf.description; version = mf.version

  header = """
  -${c_dim}--------------------------------------------------${c_reset}-
  ${c_blue}  ${author} ${c_green} ${name} v${version}${c_reset}
  ${c_yellow}  ${desc}${c_reset}
  ${c_purple}  ${url}${c_reset}
  -${c_dim}--------------------------------------------------${c_reset}-
  """.stripIndent()
  sum = summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
  foot = "-${c_dim}--------------------------------------------------${c_reset}-"
  return [header, sum, foot].join("\n")
}

def has_ext(it, extension) {
  it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

def prep_params_genome(params) {
  if (params.genomes && params.genome && !params.genomes.containsKey(params.genome))
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
  // genome parameters
  params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
  params.fasta_idx = params.genome ? params.genomes[ params.genome ].fasta_idx ?: false : false
  params.genome_bed = params.genome ? params.genomes[ params.genome ].genome_bed ?: false : false
  params.genome_sizes = params.genome ? params.genomes[ params.genome ].genome_sizes ?: false : false
  params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
  params.gff = params.genome ? params.genomes[ params.genome ].gff ?: false : false
  params.bed = params.genome ? params.genomes[ params.genome ].bed ?: false : false
  params.fna = params.genome ? params.genomes[ params.genome ].fna ?: false : false
  params.faa = params.genome ? params.genomes[ params.genome ].faa ?: false : false
  params.pgtf = params.genome ? params.genomes[ params.genome ].pgtf ?: false : false
  params.pgff = params.genome ? params.genomes[ params.genome ].pgff ?: false : false
  params.pbed = params.genome ? params.genomes[ params.genome ].pbed ?: false : false
  params.pfna = params.genome ? params.genomes[ params.genome ].pfna ?: false : false
  params.pfaa = params.genome ? params.genomes[ params.genome ].pfaa ?: false : false
  params.win11 = params.genome ? params.genomes[ params.genome ].win11 ?: false : false
  params.win56 = params.genome ? params.genomes[ params.genome ].win56 ?: false : false
  params.win127 = params.genome ? params.genomes[ params.genome ].win127 ?: false : false
  params.rcfg = params.genome ? params.genomes[ params.genome ].rcfg ?: false : false
  params.tx2gene = params.genome ? params.genomes[ params.genome ].tx2gene ?: false : false
  params.bwa_index = params.genome ? params.genomes[ params.genome ].bwa ?: false : false
  params.hisat2_index = params.genome ? params.genomes[ params.genome ].hisat2 ?: false : false
  params.hisat2_ss = params.genome ? params.genomes[ params.genome ].hisat2_ss ?: false : false
  params.star_index = params.genome ? params.genomes[ params.genome ].star ?: false : false
  params.salmon_index = params.genome ? params.genomes[ params.genome ].salmon ?: false : false
  params.bismark_index = params.genome ? params.genomes[ params.genome ].bismark ?: false : false
  params.bwa_meth_index = params.genome ? params.genomes[ params.genome ].bwa_meth ?: false : false
  params.tss_bed = params.genome ? params.genomes[ params.genome ].tss_bed ?: false : false
  params.macs_gsize = params.genome ? params.genomes[ params.genome ].macs_gsize ?: false : false
  params.blacklist = params.genome ? params.genomes[ params.genome ].blacklist ?: false : false
}

def yml_summary(summary) {
  def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
  yaml_file.text  = """
  id: 'rnaseq-summary'
  description: " - this information is collected when the pipeline is started."
  section_name: 'nfc/rnaseq Workflow Summary'
  section_href: 'https://github.com/orionzhou/nfc/rnaseq'
  plot_type: 'html'
  data: |
      <dl class=\"dl-horizontal\">
      ${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
      </dl>
    """.stripIndent()
  return yaml_file
}

def check_host() {
  def c_reset = params.monochrome_logs ? '' : "\033[0m"
  def c_white = params.monochrome_logs ? '' : "\033[0;37m"
  def c_red = params.monochrome_logs ? '' : "\033[1;91m"
  def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
  if (params.hostnames) {
    def hostname = "hostname".execute().text.trim()
    params.hostnames.each { prof, hnames ->
      hnames.each { hname ->
        if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
          log.error "====================================================\n" +
            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
            "============================================================"
        }
      }
    }
  }
}

def get_read_num(design) {
  ch_read_num = design
    .splitCsv(header:true, sep:"\t")
    .map { row -> [ row.SampleID, row.spots.toInteger() ] }
  return ch_read_num
}

def gt_upper(String str) {
  str.replace("x","_").toUpperCase().replace("_","x")
}

def get_ch_bcf(design) {
  if (params.ase) {
    ch_bcf_ase = design
      .splitCsv(header:true, sep:"\t")
      .map { row -> [ row.SampleID,
    file("${params.ase_dir}/${row.Genotype.replace('x','_').toUpperCase().replace('_','x')}.bcf"),
    file("${params.ase_dir}/${row.Genotype.replace('x','_').toUpperCase().replace('_','x')}.bcf.csi")
                    ] }
  } else {
    ch_bcf_ase = Channel.from(false)
  }
  return ch_bcf_ase
}

def on_complete1(workflow, params) {
  mf = workflow.manifest
  name = params.name ?: mf.name; author = mf.author; url = mf.homePage;
  // Set up the e-mail variables def subject = "[$name] Suceeded"
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
    email_fields['summary'] = sum
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    if (workflow.container) email_fields['summary']['Docker image'] = workflow.container
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

  // On success try attach the multiqc report
  def mqc_report = null
  try {
    if (workflow.success && !params.skip_multiqc) {
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
  //def logo_b64 = new File("$baseDir/assets/nf-core-rnaseq_logo.png").bytes.encodeBase64().toString()
  //email_html = email_html.replace('<img src="cid:nfcorepipelinelogo">', "<img src=\"data:image/png;base64, ${logo_b64}\">")
  email_html = email_html.replace('<img src="cid:nfcorepipelinelogo">', "")
  // Print to file
  def output_hf = file("${output_d}/pipeline_report.html")
  output_hf.withWriter { w -> w << email_html }
  def output_tf = file("${output_d}/pipeline_report.txt")
  output_tf.withWriter { w -> w << email_txt }

  c_reset = params.monochrome_logs ? '' : "\033[0m";
  c_purple = params.monochrome_logs ? '' : "\033[0;35m";
  c_green = params.monochrome_logs ? '' : "\033[0;32m";
  c_red = params.monochrome_logs ? '' : "\033[0;31m";

  // def nsam = ch_reads_l.getVal().size()
  log.info "[${c_purple}${name}${c_reset}] ${c_green}All samples processed${c_reset}"

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

def on_complete(workflow, params) {
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
    email_fields['summary'] = sum
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
    if (workflow.success && !params.skip_multiqc) {
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
  //def logo_b64 = new File("$baseDir/assets/nf-core-rnaseq_logo.png").bytes.encodeBase64().toString()
  //email_html = email_html.replace('<img src="cid:nfcorepipelinelogo">', "<img src=\"data:image/png;base64, ${logo_b64}\">")
  email_html = email_html.replace('<img src="cid:nfcorepipelinelogo">', "")
  // Print to file
  def output_hf = file("${output_d}/pipeline_report.html")
  output_hf.withWriter { w -> w << email_html }
  def output_tf = file("${output_d}/pipeline_report.txt")
  output_tf.withWriter { w -> w << email_txt }

  c_reset = params.monochrome_logs ? '' : "\033[0m";
  c_purple = params.monochrome_logs ? '' : "\033[0;35m";
  c_green = params.monochrome_logs ? '' : "\033[0;32m";
  c_red = params.monochrome_logs ? '' : "\033[0;31m";

  //def nsam = ch_reads_l.getVal().size()
  log.info "[${c_purple}${name}${c_reset}] ${c_green}All samples processed${c_reset}"

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

