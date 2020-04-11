include '../modules/utils.nf'

def help() {
  log.info showHeader()
  log.info """

  Usage:

  The typical command for running the pipeline is as follows:

  nextflow run $nfc/mmm

  """.stripIndent()
}

def prep_params(params, workflow) {
  if (!params.name) params.name = workflow.runName
}

def getSummary() {
  def summary = [:]
  if (workflow.revision) summary['Pipeline Release'] = workflow.revision
  summary['Run Name'] = params.name ?: workflow.runName
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
  reads = Channel.create()
  if (params.source == 'local') {
    reads = design
      .splitCsv(header:true, sep:"\t")
      .map {row -> [row.SampleID, row.paired.toBoolean(), false, 'r0','r1','r2']}
  } else if (params.source == 'sra') {
    reads = design
      .splitCsv(header:true, sep:"\t")
      .map {row -> [row.SampleID, row.paired.toBoolean(),
        row.interleaved.toBoolean(), row.r0, row.r1, row.r2]}
  }
  return reads
}


