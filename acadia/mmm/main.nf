#!/usr/bin/env nextflow
nextflow.preview.dsl = 2
include '../modules/utils.nf'
include './helper.nf'

if (params.help) { help(); exit 0 }
prep_params(params, workflow)

// set up channels
  seqdb = Channel.fromPath([params.seqdb, params.seqdb_idx], checkIfExists: true)
    .ifEmpty { exit 1, "no seqdb found: ${params.seqdb}" }
    .toList()
  fimo_bg = Channel.fromPath(params.fimo_bg, checkIfExists: true)
    .ifEmpty { exit 1, "no fimo_bg file found: ${params.fimo_bg}" }
  mtf = Channel.fromPath(params.mtf, checkIfExists: true)
    .ifEmpty { exit 1, "no motif file found: ${params.mtf}" }
  lsts = Channel.fromPath(params.lst, checkIfExists: true)
    .ifEmpty { exit 1, "no seq list found: ${params.lst}" }
    .splitCsv(header:true, sep:"\t")
    .map { row -> [ row.lid, file("${params.lstdir}/${row.lid}.txt", checkIfExists:true) ]}

include '../modules/mmm.nf'

def summary = getSummary()
log.info showHeader(summary)

workflow {
  main:
    seqret(lsts.combine(seqdb))
    fimo(seqret.out.combine(mtf).combine(fimo_bg))
    mg_fimo(fimo.out.sum.collect({it[1]}))
  publish:
    fimo.out to: "${params.outdir}/11_fimo", mode:'link', overwrite: true
    mg_fimo.out to: "${params.outdir}", mode:'link', overwrite: true
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


