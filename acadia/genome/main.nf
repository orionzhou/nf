#!/usr/bin/env nextflow
nextflow.preview.dsl = 2
include {show_header; check_host; yml_summary; has_ext} from '../modules/utils.nf'
include {help; prep_params; summary} from './helper.nf'

if (params.help) { help(); exit 0 }

prep_params(params, workflow)

// Stage config files
  ch_out_doc = file("$baseDir/docs/output.md", checkIfExists: true)

// validate inputs
  design = Channel
    .fromPath(params.design, checkIfExists: true)
    .ifEmpty { exit 1, "genome design table missing: ${params.design}" }

include {genome} from '../modules/genome.nf'
//include {version; outdoc} from '../modules/genome.nf'

def sum = summary()
log.info show_header(sum)
check_host()

workflow {
  main:
    pick = file(params.pick).readLines()
    genome(design, pick)
    //version()
    //outdoc(ch_out_doc)
  //publish:
}

workflow.onComplete {
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
  def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
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
  // log.info "[${c_purple}${name}${c_reset}] ${c_green}All samples processed${c_reset}"

  if (workflow.stats.ignoredCount > 0 && workflow.success) {
    log.info "- ${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}"
    log.info "- ${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}"
    log.info "- ${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}"
  }

  if (workflow.success) {
    log.info "[${c_purple}${name}${c_reset}] ${c_green}Pipeline completed successfully${c_reset}"
  } else {
    check_host()
    log.info "[${c_purple}${name}${c_reset}] ${c_red}Pipeline completed with errors${c_reset}"
  }
}


