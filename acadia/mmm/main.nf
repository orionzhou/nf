#!/usr/bin/env nextflow
nextflow.preview.dsl = 2
include {show_header; check_host; yml_summary; has_ext} from '../modules/utils.nf'
include {help; prep_params; summary} from './helper.nf'

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
  mtfs = Channel.fromPath(params.mtf_lst, checkIfExists: true)
    .ifEmpty { exit 1, "no motif list found: ${params.mtf}" }
    .splitCsv(header:false)
    .map { row -> [ row[0] ] }
  lsts = Channel.fromPath(params.lst, checkIfExists: true)
    .ifEmpty { exit 1, "no seq list found: ${params.lst}" }
    .splitCsv(header:true, sep:"\t")
    .map { row -> [ row.lid, file("${params.lstdir}/${row.lid}.txt", checkIfExists:true) ]}
  lst_pairs = Channel.fromPath(params.lst, checkIfExists: true)
    .ifEmpty { exit 1, "no seq list found: ${params.lst}" }
    .splitCsv(header:true, sep:"\t")
    .map { row -> [ row.lid, row.clid ] }
  bg_lsts = Channel.fromPath(params.bg_lst, checkIfExists: true)
    .ifEmpty { exit 1, "no bg list found: ${params.lst}" }
    .splitCsv(header:true, sep:"\t")
    .map { row -> [ row.lid, file("${params.bg_lstdir}/${row.lid}.txt", checkIfExists:true) ]}


include {mmk; mmd} from '../modules/mmm.nf'

def sum = summary()
log.info show_header(sum)

workflow {
  main:
    //mmk(seqdb, mtfs, mtf, fimo_bg)
    mmd(seqdb, lsts, bg_lsts, lst_pairs)
  //publish:
    //mmk.out.fimo to: "${params.outdir}/11_fimo_raw", mode:'copy', overwrite: true
    //mmk.out.fimo2 to: "${params.outdir}/12_fimo_sum", mode:'copy', overwrite: true
    //mmk.out.mg_fimo to: "${params.outdir}", mode:'copy', overwrite: true
    //mmd.out.seqs to: "${params.outdir}/21_seqs", mode:'copy', overwrite: true
    //mmd.out.dreme to: "${params.outdir}/22_dreme", mode:'copy', overwrite: true
    //mmd.out.mg_dreme to: "${params.outdir}", mode:'copy', overwrite: true
}

workflow.onComplete {
  mf = workflow.manifest
  name = params.name ?: mf.name; author = mf.author; url = mf.homePage;

  c_reset = params.monochrome_logs ? '' : "\033[0m";
  c_purple = params.monochrome_logs ? '' : "\033[0;35m";
  c_green = params.monochrome_logs ? '' : "\033[0;32m";
  c_red = params.monochrome_logs ? '' : "\033[0;31m";

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


