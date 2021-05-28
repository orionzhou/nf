#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include {show_header; check_host; yml_summary; has_ext} from '../modules/utils.nf'
include {help; prep_params; summary} from './helper.nf'

if (params.help) { help(); exit 0 }

prep_params(params, workflow)

// Stage config files
  ch_out_doc = file("$baseDir/docs/output.md", checkIfExists: true)

// validate inputs
  comps = Channel.fromPath(params.comps, checkIfExists: true)
      .ifEmpty { exit 1, "comparison file not found: ${params.comps}" }
      .splitCsv(header:false, sep:',')
      .map {r -> [r[0], r[1]]}
  qrys = comps.map{r -> r[0]}.unique().map{r->[r,params.genomes[r]]}
  tgts = comps.map{r -> r[1]}.unique().map{r->[r,params.genomes[r]]}
  // qry seq
    qry_fas = qrys
      .map {r -> [r[0], file(r[1].fasta, checkIfExists:true)]}
    qry_gbed = qrys
      .map {r -> [r[0], file(r[1].genome_bed, checkIfExists:true)]}
    qry_sizes = qrys
      .map {r -> [r[0], file(r[1].genome_sizes, checkIfExists:true)]}
    qry_gap = qrys
      .map {r -> [r[0], file(r[1].gap_bed, checkIfExists:true)]}
    qry_2bit = qrys
      .map {r -> [r[0], file(r[1].blat, checkIfExists:true)]}
    qry_gatk = qrys
      .map {r -> [r[0], file(r[1].gatk, checkIfExists:true)]}
  // qry annotation
    qry_gff = qrys
      .map {r -> [r[0], file(r[1].gff, checkIfExists:true)]}
    qry_pgff = qrys
      .map {r -> [r[0], file(r[1].pgff, checkIfExists:true)]}
    qry_fna = qrys
      .map {r -> [r[0], file(r[1].fna, checkIfExists:true)]}
    qry_faa = qrys
      .map {r -> [r[0], file(r[1].faa, checkIfExists:true)]}
    qry_pfna = qrys
      .map {r -> [r[0], file(r[1].pfna, checkIfExists:true)]}
    qry_pfaa = qrys
      .map {r -> [r[0], file(r[1].pfaa, checkIfExists:true)]}
    qry_pbed = qrys
      .map {r -> [r[0], file(r[1].pbed, checkIfExists:true)]}
  // tgt seq
    tgt_fas = tgts
      .map {r -> [r[0], file(r[1].fasta, checkIfExists:true)]}
    tgt_gbed = tgts
      .map {r -> [r[0], file(r[1].genome_bed, checkIfExists:true)]}
    tgt_sizes = tgts
      .map {r -> [r[0], file(r[1].genome_sizes, checkIfExists:true)]}
    tgt_gap = tgts
      .map {r -> [r[0], file(r[1].gap_bed, checkIfExists:true)]}
    tgt_2bit = tgts
      .map {r -> [r[0], file(r[1].blat, checkIfExists:true)]}
    tgt_gatk = tgts
      .map {r -> [r[0], file(r[1].gatk, checkIfExists:true)]}
  // tgt annotation
    tgt_gff = tgts
      .map {r -> [r[0], file(r[1].gff, checkIfExists:true)]}
    tgt_pgff = tgts
      .map {r -> [r[0], file(r[1].pgff, checkIfExists:true)]}
    tgt_fna = tgts
      .map {r -> [r[0], file(r[1].fna, checkIfExists:true)]}
    tgt_faa = tgts
      .map {r -> [r[0], file(r[1].faa, checkIfExists:true)]}
    tgt_pfna = tgts
      .map {r -> [r[0], file(r[1].pfna, checkIfExists:true)]}
    tgt_pfaa = tgts
      .map {r -> [r[0], file(r[1].pfaa, checkIfExists:true)]}
    tgt_pbed = tgts
      .map {r -> [r[0], file(r[1].pbed, checkIfExists:true)]}
  tgt_blastp = tgts
    .map {r -> [r[0], file(r[1].blastp, checkIfExists:true)]}

include {wgc} from '../modules/wgc.nf'
include {syn} from '../modules/syntelog.nf'

def sum = summary()
log.info show_header(sum)
check_host()

workflow {
  main:
    //version()
    wgc(comps,
        qry_fas, qry_gbed, qry_sizes, qry_gap, qry_2bit,
        qry_gff, qry_pgff, qry_gatk,
        tgt_fas, tgt_gbed, tgt_sizes, tgt_gap, tgt_2bit,
        tgt_gff, tgt_pgff, tgt_gatk)
    syn(comps, qry_pfaa, qry_pbed,
        tgt_pfaa, tgt_pbed, tgt_blastp)
    //outdoc(ch_out_doc)
  //publish:
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


