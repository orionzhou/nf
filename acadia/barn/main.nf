#!/usr/bin/env nextflow
nextflow.preview.dsl = 2
include '../modules/utils.nf'
include './helper.nf'

if (params.help) { help(); exit 0 }
prep_params(params, workflow)

// validate inputs

// set up channels
  nofile = Channel.fromPath("NO_FILE")
  design = Channel
    .fromPath(params.design, checkIfExists: true)
    .ifEmpty { exit 1, "sample design table missing: ${params.design}" }
  reads = get_reads(design)

include fqc from '../modules/fastq.nf'
include {fqd; fqz; fqv; mqc; upd} from '../modules/barn.nf'

def summary = getSummary()
log.info showHeader(summary)

workflow {
  main:
    clean_reads = Channel.empty()
    if( params.local ) {
      reads | fqz | fqv
      clean_reads = fqv.out
    } else {
      reads | fqd
      clean_reads = fqd.out
    }
    clean_reads | fqc
    mqc(fqc.out.zip.collect())
    upd(design, mqc.out)
  publish:
    clean_reads to: "${params.fq_dir}/${params.name}", mode:'copy', overwrite: true
    upd.out to: "${params.meta_dir}", mode:'copy', overwrite: true
}

workflow.onComplete {
  mf = workflow.manifest
  name = params.name ?: mf.name; author = mf.author; url = mf.homePage;

  c_reset = params.monochrome_logs ? '' : "\033[0m";
  c_purple = params.monochrome_logs ? '' : "\033[0;35m";
  c_green = params.monochrome_logs ? '' : "\033[0;32m";
  c_red = params.monochrome_logs ? '' : "\033[0;31m";

  log.info "[${c_purple}${name}${c_reset}] ${c_green}samples processed${c_reset}"

  if (workflow.stats.ignoredCount > 0 && workflow.success) {
    log.info "- ${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}"
    log.info "- ${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}"
    log.info "- ${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}"
  }

  if (workflow.success) {
    log.info "[${c_purple}${name}${c_reset}] ${c_green}Pipeline completed successfully${c_reset}"
  } else {
    log.info "[${c_purple}${name}${c_reset}] ${c_red}Pipeline completed with errors${c_reset}"
  }
}


