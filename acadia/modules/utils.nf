def showHeader(summary) {
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

def hasExtension(it, extension) {
  it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

def prep_params_genome(params) {
  if (params.genomes && params.genome && !params.genomes.containsKey(params.genome))
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
  // genome parameters
  params.star_index = params.genome ? params.genomes[ params.genome ].star ?: false : false
  params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
  params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
  params.gff = params.genome ? params.genomes[ params.genome ].gff ?: false : false
  params.bed12 = params.genome ? params.genomes[ params.genome ].bed12 ?: false : false
  params.win11 = params.genome ? params.genomes[ params.genome ].win11 ?: false : false
  params.win56 = params.genome ? params.genomes[ params.genome ].win56 ?: false : false
  params.win127 = params.genome ? params.genomes[ params.genome ].win127 ?: false : false
  params.transcript_fasta = params.genome ? params.genomes[ params.genome ].transcript_fasta ?: false : false
  params.rcfg = params.genome ? params.genomes[ params.genome ].rcfg ?: false : false
  params.tx2gene = params.genome ? params.genomes[ params.genome ].tx2gene ?: false : false
  params.hisat2_index = params.genome ? params.genomes[ params.genome ].hisat2 ?: false : false
  params.hisat2_ss = params.genome ? params.genomes[ params.genome ].hisat2_ss ?: false : false
  params.salmon_index = params.genome ? params.genomes[ params.genome ].salmon ?: false : false
}

def create_workflow_summary(summary) {
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

def checkHostname() {
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

def get_ch_reads(design) {
  ch_reads_se = design
    .splitCsv(header:true, sep:"\t")
    .filter { it.paired != 'TRUE' }
    .map { row -> [ row.SampleID, false, [
      file("${params.seqdir}/${params.name}/${row.SampleID}.fq.gz", checkIfExists:true)
      ] ] }
  ch_reads_pe = design
    .splitCsv(header:true, sep:"\t")
    .filter { it.paired == 'TRUE' }
    .map { row -> [ row.SampleID, true, [
      file("${params.seqdir}/${params.name}/${row.SampleID}_1.fq.gz", checkIfExists:true),
      file("${params.seqdir}/${params.name}/${row.SampleID}_2.fq.gz", checkIfExists:true)
      ] ] }
  return ch_reads_se.concat(ch_reads_pe)
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



