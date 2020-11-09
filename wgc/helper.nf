include {show_header; } from '../modules/utils.nf'

def help() {
  log.info show_header()

  log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

      nextflow run $nf/wgc -profile msi

    Mandatory arguments:
      -qry [str]                 query genome
      -tgt [str]                 target genome
      -profile [str]             Configuration profile to use. Can use multiple (comma separated)
                                 Available: msi, msi_mesabi, msi_mangi

  """.stripIndent()
}

def prep_params(params, workflow) {
  if (!params.name) params.name = workflow.runName
  //if (!params.qry || !params.genomes.containsKey(params.qry))
    //exit 1, "The query genome '${params.qry}' is not available. Current available genomes are ${params.genomes.keySet().join(", ")}"
  //if (!params.tgt || !params.genomes.containsKey(params.tgt))
    //exit 1, "The target genome '${params.tgt}' is not available. Current available genomes are ${params.genomes.keySet().join(", ")}"
}

def summary() {
  def summary = [:]
  if (workflow.revision) summary['Pipeline Release'] = workflow.revision
  summary['Run Name'] = params.name ?: workflow.runName
  summary['Comparisons'] = "${params.comps}"
  summary['Max Resources'] = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
  if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
  summary['Output dir'] = params.outdir
  summary['Launch dir'] = workflow.launchDir
  summary['Working dir'] = workflow.workDir
  summary['Script dir'] = workflow.projectDir
  summary['User'] = workflow.userName
  if (workflow.profile == 'awsbatch') {
    summary['AWS Region']     = params.awsregion
    summary['AWS Queue']      = params.awsqueue
  }
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



