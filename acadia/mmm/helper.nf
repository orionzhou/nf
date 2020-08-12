include {show_header; prep_params_genome} from '../modules/utils.nf'

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
  prep_params_genome(params)
}

def summary() {
  def summary = [:]
  if (workflow.revision) summary['Pipeline Release'] = workflow.revision
  summary['Run Name'] = params.name ?: workflow.runName
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

