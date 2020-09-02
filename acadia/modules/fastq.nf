include {has_ext} from '../modules/utils.nf'
process fqd {
  label 'low_memory'
  tag "$id"
  publishDir "${params.outdir}/00_fastq", mode:'link', overwrite: true,
    saveAs: { fn -> params.save_fastq ? "$fn" : null }

  input:
  tuple id, paired, acc

  output:
  tuple id, paired, path("${id}_R?.fq.gz")

  script:
  def mem = task.memory
  if( paired )
    """
    fasterq-dump --split-files -e ${task.cpus} -m ${task.memory.toGiga()}GB \
        -O ./ -t $TMPDIR ${acc}
    pigz -p ${task.cpus} --fast -c ${acc}_1.fastq > ${id}_R1.fq.gz
    pigz -p ${task.cpus} --fast -c ${acc}_2.fastq > ${id}_R2.fq.gz
    """
  else
    """
    fasterq-dump --split-files -e ${task.cpus} -m ${task.memory.toGiga()}GB \
        -O ./ -t $TMPDIR ${acc}
    pigz -p ${task.cpus} --fast -c ${acc}.fastq > ${id}_R0.fq.gz
    """
}

process fqz {
  label 'low_memory'
  tag "$id"

  input:
  tuple id, paired, path(r0), path(r1), path(r2), gz0, gz1

  output:
  tuple id, paired, path("${id}_R0.fq.gz"), path("${id}_R1.fq.gz"), path("${id}_R2.fq.gz")

  script:
  if( paired && gz1 )
    """
    ln -f $r1 ${id}_R1.fq.gz
    ln -f $r2 ${id}_R2.fq.gz
    touch ${id}_R0.fq.gz
    """
  else if( paired && !gz1 )
    """
    pigz -p ${task.cpus} -c $r1 > ${id}_R1.fq.gz
    pigz -p ${task.cpus} -c $r2 > ${id}_R2.fq.gz
    touch ${id}_R0.fq.gz
    """
  else if( !paired && gz0 )
    """
    ln -f $r0 ${id}_R0.fq.gz
    touch ${id}_R1.fq.gz ${id}_R2.fq.gz
    #cp -fL $r0 ${id}_R0.fq.gz
    """
  else
    """
    pigz -p ${task.cpus} -c $r0 > ${id}_R0.fq.gz
    touch ${id}_R1.fq.gz ${id}_R2.fq.gz
    """
}

process fqv {
  label 'low_memory'
  tag "$id"
  publishDir "${params.outdir}/00_fastq", mode:'link', overwrite: true,
    saveAs: { fn -> params.save_fastq ? "$fn" : null }

  input:
  tuple id, paired, path(r0), path(r1), path(r2)

  output:
  tuple id, paired, path("${id}_R?.fq.gz", includeInputs: true)

  script:
  if( paired && params.interleaved )
    """
    rm $r1 $r2
    zcat $r0 |\
      deinterleave_fastq.sh $r1 $r2 ${task.cpus} compress
    rm $r0
    """
  else if ( paired && !params.interleaved )
    """
    rm $r0
    """
  else
    """
    rm $r1 $r2
    """
}


process fqc {
  label 'low_memory'
  tag "$id"
  publishDir "${params.outdir}/01_fastqc", mode:'link', overwrite:'true',
    saveAs: { fn ->
      fn.endsWith(".zip") ? "zips/$fn" : "$fn"
    }

  when:
  !params.skip_qc && !params.skip_fastqc

  input:
  tuple id, paired, path(reads)

  output:
  path "${id}_*_fastqc.zip", emit: zip
  path "${id}_*_fastqc.html", emit: html

  script:
  """
  fastqc --quiet --threads ${task.cpus} --noextract --format fastq -o . $reads
  """
}

process trim {
  label 'low_memory'
  tag "$id"
  publishDir "${params.outdir}/02_trim_galore", mode:'link', overwrite:'true',
    saveAs: { fn ->
      if (fn.endsWith(".html")) "fastqc/$fn"
      else if (fn.endsWith(".zip")) "fastqc/zips/$fn"
      else if (fn.endsWith("trimming_report.txt")) "logs/$fn"
      else params.save_trimmed ? fn : null
    }

  input:
  tuple id, paired, path(reads)

  output:
  tuple id, paired, path("*fq.gz"), emit: reads
  path "*trimming_report.txt", emit: log
  path "*_fastqc.zip", emit: fastqc_zip
  path "*_fastqc.html", emit: fastqc_html

  script:
  c_r1 = params.clip_r1 > 0 ? "--clip_r1 ${clip_r1}" : ''
  c_r2 = params.clip_r2 > 0 ? "--clip_r2 ${clip_r2}" : ''
  tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${three_prime_clip_r1}" : ''
  tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${three_prime_clip_r2}" : ''
  nextseq = params.trim_nextseq > 0 ? "--nextseq ${params.trim_nextseq}" : ''
  pair_str = paired ? "--paired" : ""
  clip_str = paired ? "$c_r1 $c_r2" : c_r1
  tpc_str = paired ? "$tpc_r1 $tpc_r2" : tpc_r1
  """
  trim_galore $pair_str --fastqc --gzip $clip_str $tpc_str $nextseq $reads
  """
}

process mqc {
  label 'low_memory'
  tag "${params.name}"

  input:
  path(fqcs)

  output:
  path("mqc.txt")

  script:
  """
  multiqc -f -o . ${fqcs}
  mv multiqc_data/multiqc_fastqc.txt mqc.txt
  """
}

process upd {
  label 'low_memory'
  tag "${params.name}"
  conda "$NXF_CONDA_CACHEDIR/r"
  publishDir "${params.outdir}", mode:'copy', overwrite:'true'
  publishDir "${params.qcdir}/${params.genome}/${params.name}", mode:'copy', overwrite:'true'

  input:
  path(design)
  path(mqc)

  output:
  path("00.meta.tsv")

  script:
  """
  mv $design in.tsv
  samplelist_addstat.R in.tsv $mqc 00.meta.tsv
  """
}

def get_reads(design) {
  reads = Channel.empty()
  if( params.source == 'local' ) {
    reads_se = design
      .splitCsv(header:true, sep:"\t")
      .filter { it.paired != 'TRUE' }
      .map {row -> [row.SampleID, row.paired.toBoolean(),
        file(row.r0), file('f1'), file('f2'), has_ext(row.r0, ".gz"), null]}
    reads_pe = design
      .splitCsv(header:true, sep:"\t")
      .filter { it.paired == 'TRUE' }
      .map {row -> [row.SampleID, row.paired.toBoolean(),
        file('f0'), file(row.r1), file(row.r2), null, has_ext(row.r1, ".gz")]}
    reads = reads_se.concat(reads_pe)
  } else if (params.source == 'sra') {
    if (params.lib in ['chipseq','dapseq']) {
      reads = design
        .splitCsv(header:true, sep:"\t")
        .map {row -> [row.SampleID, row.paired.toBoolean(), row.acc]}
    } else {
      reads = design
        .splitCsv(header:true, sep:"\t")
        .map {row -> [row.SampleID, row.paired.toBoolean(), row.SampleID]}
    }
  } else {
    exit 1, "unknown sequence source: ${params.source}"
  }
  return reads
}

def get_reads_o1(design) {
  reads_se = design
    .splitCsv(header:true, sep:"\t")
    .filter { it.paired != 'TRUE' }
    .map { row -> [ row.SampleID, false, [
      file("${params.seqdir}/${params.name}/${row.SampleID}.fq.gz", checkIfExists:true)
      ] ] }
  reads_pe = design
    .splitCsv(header:true, sep:"\t")
    .filter { it.paired == 'TRUE' }
    .map { row -> [ row.SampleID, true, [
      file("${params.seqdir}/${params.name}/${row.SampleID}_1.fq.gz", checkIfExists:true),
      file("${params.seqdir}/${params.name}/${row.SampleID}_2.fq.gz", checkIfExists:true)
      ] ] }
  return reads_se.concat(reads_pe)
}
def get_reads_o2(design) {
  reads_se = design
    .splitCsv(header:true, sep:"\t")
    .filter { it.paired != 'TRUE' }
    .map { row -> [ row.SampleID, false, [
      file("${params.seqdir}/${params.name}/${row.SampleID}.fq.gz", checkIfExists:true)
      ] ] }
  reads_pe = design
    .splitCsv(header:true, sep:"\t")
    .filter { it.paired == 'TRUE' }
    .map { row -> [ row.SampleID, true, [
      //file("${params.seqdir}/${params.name}/${row.SampleID}_{1,R1}.fq.gz", checkIfExists:true),
      //file("${params.seqdir}/${params.name}/${row.SampleID}_{2,R2}.fq.gz", checkIfExists:true)
      file("${params.seqdir}/${params.name}/${row.SampleID}_1.fq.gz", checkIfExists:true),
      file("${params.seqdir}/${params.name}/${row.SampleID}_2.fq.gz", checkIfExists:true)
      ] ] }
  return reads_se.concat(reads_pe)
}

workflow fq {
  take: design
  main:
    raw_read_list = get_reads(design)
    raw_reads = Channel.empty()
    if( params.source == 'local' ) {
      raw_read_list | fqz | fqv
      raw_reads = fqv.out
    } else {
      raw_read_list | fqd
      raw_reads = fqd.out
    }
    raw_reads | (fqc & trim)
    mqc(fqc.out.zip.collect())
    upd(design, mqc.out)
  emit:
    raw_reads = raw_reads
    trim_reads = trim.out.reads
    trim_log = trim.out.log
    raw_fqc = fqc.out.zip
    trim_fqc = trim.out.fastqc_zip
    readlist = upd.out
}



