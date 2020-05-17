include '../modules/utils.nf'

def help() {
  log.info showHeader()
  log.info """

  Usage:

  The typical command for running the pipeline is as follows:

  nextflow run $nfc/barn

  """.stripIndent()
}

def prep_params(params, workflow) {
  if (!params.name) params.name = workflow.runName
  design_dir = params.local ? params.local_dir : params.sra_dir
  params.design = "${design_dir}/${params.name}.tsv"
}

def getSummary() {
  def summary = [:]
  if (workflow.revision) summary['Pipeline Release'] = workflow.revision
  summary['Run Name'] = params.name ?: workflow.runName
  design_dir = params.local ? params.local_dir : params.sra_dir
  params.design = "${design_dir}/${params.name}.tsv"
  summary['Design'] = params.design
  summary['Config Profile'] = workflow.profile
  if (params.config_profile_description) summary['Config Description'] = params.config_profile_description
  if (params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
  if (params.config_profile_url)         summary['Config URL']         = params.config_profile_url
  if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
  }
  return summary
}

def get_reads(design) {
  reads = Channel.empty()
  if (params.local) {
    reads = design
      .splitCsv(header:true, sep:"\t")
      .map {row -> [row.SampleID, row.paired.toBoolean(),
        row.interleaved.toBoolean(), file(row.r0), file(row.r1), file(row.r2),
        hasExtension(row.r0, ".gz"), hasExtension(row.r1, ".gz")]}
  } else {
    reads = design
      .splitCsv(header:true, sep:"\t")
      .map {row -> [row.SampleID, row.paired.toBoolean()]}
  }
  return reads
}


