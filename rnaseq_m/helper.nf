include {show_header; prep_params_genome} from '../modules/utils.nf'

def help() {
  log.info showHeader()
  log.info """

  Usage:

  The typical command for running the pipeline is as follows:

  nextflow run rnaseq -params-file genomes.yml --genome Zmays_B73v5 -profile conda,msi

  Mandatory arguments:
    -design                       tab-separated sample design table
    -profile                      Configuration profile to use. Can use multiple (comma separated)
                                  Available: conda, docker, singularity, awsbatch, test and more.

  Generic:
    --source                      ['sra', 'local', 's3', 'mixed']
    --paired                      ['SE', 'PE', 'mixed']
    --stranded                    ['no', 'forward', 'rerverse']

  References:                     If not specified in the configuration file or you wish to overwrite any of the references.
    --genome                      Name of iGenomes reference
    --star_index                  Path to STAR index
    --hisat2_index                Path to HiSAT2 index
    --salmon_index                Path to Salmon index
    --fasta                       Path to genome fasta file
    --transcript_fasta            Path to transcript fasta file
    --gtf                         Path to GTF file
    --gff                         Path to GFF3 file
    --bed12                       Path to bed12 file
    --rcfg                        Path to rconfig file
    --saveReference               Save the generated reference files to the results directory
    --gencode                     Use fc_group_features_type = 'gene_type' and pass '--gencode' flag to Salmon

  Trimming:
    --skip_trimming               Skip Trim Galore step
    --clip_r1 [int]               Instructs Trim Galore to remove bp from the 5' end of read 1 (or single-end reads)
    --clip_r2 [int]               Instructs Trim Galore to remove bp from the 5' end of read 2 (paired-end reads only)
    --three_prime_clip_r1 [int]   Instructs Trim Galore to remove bp from the 3' end of read 1 AFTER adapter/quality trimming has been performed
    --three_prime_clip_r2 [int]   Instructs Trim Galore to remove bp from the 3' end of read 2 AFTER adapter/quality trimming has been performed
    --trim_nextseq [int]          Instructs Trim Galore to apply the --nextseq=X option, to trim based on quality after removing poly-G tails
    --pico                        Sets trimming and standedness settings for the SMARTer Stranded Total RNA-Seq Kit - Pico Input kit. Equivalent to: --forwardStranded --clip_r1 3 --three_prime_clip_r2 3
    --save_trimmed                 Save trimmed FastQ file intermediates

  Ribosomal RNA removal:
    --removeRiboRNA               Removes ribosomal RNA using SortMeRNA
    --save_nonrRNA_reads          Save FastQ file intermediates after removing rRNA
    --rRNA_database_manifest      Path to file that contains file paths for rRNA databases, optional

  Alignment:
    --aligner                     Specifies the aligner to use (available are: 'hisat2', 'star')
    --stringtie_ignore_gtf        Perform reference-guided de novo assembly of transcripts using StringTie i.e. dont restrict to those in GTF file
    --seq_center                  Add sequencing center in @RG line of output BAM header
    --saveBAM                     Save the BAM files from the aligment step - not done by default
    --save_unmapped               Save unaligned reads from either STAR, HISAT2 or Salmon to extra output files.
    --skipAlignment               Skip alignment altogether (usually in favor of pseudoalignment)
    --percent_aln_skip            Percentage alignment below which samples are removed from further processing. Default: 5%

  Read Counting:
    --fc_extra_attributes         Define which extra parameters should also be included in featureCounts (default: 'gene_name')
    --fc_group_features           Define the attribute type used to group features. (default: 'gene_id')
    --fc_count_type               Define the type used to assign reads. (default: 'exon')
    --fc_group_features_type      Define the type attribute used to group features based on the group attribute (default: 'gene_biotype')
    --mapQuality                  Minimum read mapping quality (default: 20)

  QC:
    --skipQC                      Skip all QC steps apart from MultiQC
    --skipFastQC                  Skip FastQC
    --skipPreseq                  Skip Preseq
    --skipDupRadar                Skip dupRadar (and Picard MarkDuplicates)
    --skipQualimap                Skip Qualimap
    --skipBiotypeQC               Skip Biotype QC
    --skipRseQC                   Skip RSeQC
    --skipEdgeR                   Skip edgeR MDS plot and heatmap
    --skipMultiQC                 Skip MultiQC

  Other options:
    --sampleLevel                 Used to turn off the edgeR MDS and heatmap. Set automatically when running on fewer than 3 samples
    --outdir                      The output directory where the results will be saved
    -w/--work-dir                 The temporary directory where intermediate data will be saved
    --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
    --email_on_fail               Same as --email, except only send mail if the workflow is not successful
    --max_multiqc_email_size      Threshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB)
    -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic

  AWSBatch options:
    --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
    --awsregion                   The AWS Region for your AWS Batch job to run on
  """.stripIndent()
}

def prep_params(params, workflow) {
  if (!params.design) exit 1, "no design / meta file specified"
  if (!params.name) params.name = workflow.runName
  params.genome_list = params.genome.split(",").collect()
  for (genome in params.genome_list ) {
    if (!params.genomes.containsKey(genome))
      exit 1, "The provided genome '${genome}' is not available. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
  }
  if (params.read_type != 'illumina' && params.read_type != 'nanopore')
    exit 1, "Invalid read type: ${params.read_type}. Valid options: 'illumina', 'nanopore'"
  if (params.trimmer != 'no' && params.trimmer != 'trim_galore')
    exit 1, "Invalid trimmer option: ${params.trimmer}. Valid options: 'no', 'trim_galore'"
  // Preset trimming options
  if (params.pico) {
    params.clip_r1 = 3
    params.clip_r2 = 0
    params.three_prime_clip_r1 = 0
    params.three_prime_clip_r2 = 3
    params.stranded = 'forward'
  }
  // aligner
  if (params.aligner != 'star' && params.aligner != 'hisat2' && params.aligner != 'minimap2')
    exit 1, "Invalid aligner option: ${params.aligner}. Valid options: 'star', 'hisat2', 'minimap2'"
  params.biotype = params.gencode ? "gene_type" : params.fc_group_features_type

  // uppmax
  if (workflow.profile == 'uppmax' || workflow.profile == 'uppmax-devel')
    if (!params.project) exit 1, "No UPPMAX project ID found! Use --project"
  // awsbatch
  if (workflow.profile == 'awsbatch') {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    // Check outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
    // Prevent trace files to be stored on S3 since S3 does not support rolling files.
    if (workflow.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
  }
}

def summary() {
  def summary = [:]
  if (workflow.revision) summary['Pipeline Release'] = workflow.revision
  summary['Run Name'] = params.name ?: workflow.runName
  summary['Design'] = params.design
  summary['Sequence'] = "source [$params.source]; type [$params.read_type]; paired [$params.paired]; stranded [$params.stranded]"
  summary['Genome'] = params.genome
  summary['Trimmer'] = params.trimmer
  //summary['Trimming'] = "5'R1: $params.clip_r1 / 5'R2: $params.clip_r2 / 3'R1: $params.three_prime_clip_r1 / 3'R2: $params.three_prime_clip_r2 / NextSeq Trim: $params.trim_nextseq"
  //summary['Remove rRNA'] = params.removeRiboRNA
  if (params.pico) summary['Library Prep'] = "SMARTer Stranded Total RNA-Seq Kit - Pico Input"
  summary['Aligner'] = params.aligner
  summary['Run StringTie'] = params.run_stringtie
  summary['Run Salmon'] = params.run_salmon
  if (params.gencode) summary['GENCODE'] = params.gencode
  if (params.stringtie_ignore_gtf) summary['StringTie Ignore GTF'] = params.stringtie_ignore_gtf
  // if (params.fc_group_features_type) summary['Biotype GTF field'] = params.biotype
  summary['Save prefs'] = "Raw FastQ: " + (params.save_fastq ? "T" : "F")
  summary['Save prefs'] += " / Trimmed FastQ: " + (params.save_trimmed ? 'T' : 'F')
  summary['Save prefs'] += " / Final BAM: " + (params.saveBAM ? 'T' : 'F')
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
    summary['MultiQC maxsize']   = params.max_multiqc_email_size
  }
  return summary
}





