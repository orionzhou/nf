include {has_ext} from '../modules/utils.nf'
process fqu {
  label 'low_memory'
  tag "$id"

  input:
  tuple val(id), val(source), val(paired), val(r0), val(r1), val(r2)

  output:
  tuple val(id), val(paired), path("${id}_R?.fq.gz")

  script:
  mem = task.memory
  opt = params.interleaved ? "--interleaved" : ""
  reads = (paired == 'SE' || source == 'sra') ? "--r0 $r0" : "--r1 $r1 --r2 $r2"
  """
  $baseDir/bin/nf.fastq.py $id --source $source --paired $paired $opt \
    --cpu ${task.cpus} --mem ${task.memory.toGiga()} --tmp ${NXF_TEMP} \
    $reads
  """
}

process fqd {
  label 'low_memory'
  tag "$id"
  publishDir "${params.outdir}/00_fastq", mode:'copy', overwrite: true,
    saveAs: { fn -> params.save_fastq ? "$fn" : null }

  input:
  tuple val(id), val(paired), val(acc)

  output:
  tuple val(id), val(paired), path("${id}_R?.fq.gz")

  script:
  def mem = task.memory
  if( paired == 'PE' )
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
  tuple val(id), val(paired), path(r0), path(r1), path(r2), val(gz0), val(gz1)

  output:
  tuple val(id), val(paired), path("${id}_R0.fq.gz"), path("${id}_R1.fq.gz"), path("${id}_R2.fq.gz")

  script:
  if( paired =='PE' && gz1 )
    """
    ln -f $r1 ${id}_R1.fq.gz
    ln -f $r2 ${id}_R2.fq.gz
    touch ${id}_R0.fq.gz
    """
  else if( paired == 'PE' && !gz1 )
    """
    pigz -p ${task.cpus} -c $r1 > ${id}_R1.fq.gz
    pigz -p ${task.cpus} -c $r2 > ${id}_R2.fq.gz
    touch ${id}_R0.fq.gz
    """
  else if( !paired == 'PE' && gz0 )
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
  publishDir "${params.outdir}/00_fastq", mode:'copy', overwrite: true,
    saveAs: { fn -> params.save_fastq ? "$fn" : null }

  input:
  tuple val(id), val(paired), path(r0), path(r1), path(r2)

  output:
  tuple val(id), val(paired), path("${id}_R?.fq.gz", includeInputs: true)

  script:
  if( paired =='PE' && params.interleaved )
    """
    rm $r1 $r2
    zcat $r0 |\
      deinterleave_fastq.sh $r1 $r2 ${task.cpus} compress
    rm $r0
    """
  else if ( paired == 'PE' && !params.interleaved )
    """
    rm $r0
    """
  else
    """
    rm $r1 $r2
    """
}

process fqs {
  label 'low_memory'
  tag "$id"

  input:
  tuple val(id), val(paired), val(r0), val(r1), val(r2)

  output:
  tuple val(id), val(paired), path("${id}_R?.fq.gz", includeInputs: true)

  script:
  if( paired == 'PE')
    """
    s3cmd get $r1 ${id}_R1.fq.gz
    s3cmd get $r2 ${id}_R2.fq.gz
    """
  else
    """
    s3cmd get $r0 ${id}_R0.fq.gz
    """
}


process fqc {
  label 'low_memory'
  tag "$id"
  publishDir "${params.outdir}/01_fastqc", mode:'copy', overwrite:'true',
    saveAs: { fn ->
      fn.endsWith(".zip") ? "zips/$fn" : "$fn"
    }

  when:
  !params.skip_qc && !params.skip_fastqc

  input:
  tuple val(id), val(paired), path(reads)

  output:
  path "${id}_*_fastqc.zip", emit: zip
  path "${id}_*_fastqc.html", emit: html

  script:
  """
  fastqc --quiet --threads ${task.cpus} --noextract --format fastq -o . $reads
  """
}

process trim_galore {
  label 'low_memory'
  tag "$id"
  publishDir "${params.outdir}/02_trim_galore", mode:'copy', overwrite:'true',
    saveAs: { fn ->
      if (fn.endsWith(".html")) "fastqc/$fn"
      else if (fn.endsWith(".zip")) "fastqc/zips/$fn"
      else if (fn.endsWith("trimming_report.txt")) "logs/$fn"
      else params.save_trimmed ? fn : null
    }

  input:
  tuple val(id), val(paired), path(reads)

  output:
  tuple val(id), val(paired), path("*fq.gz"), emit: reads
  path "*trimming_report.txt", emit: log
  path "*_fastqc.zip", emit: fastqc_zip
  path "*_fastqc.html", emit: fastqc_html

  script:
  c_r1 = params.clip_r1 > 0 ? "--clip_r1 ${clip_r1}" : ''
  c_r2 = params.clip_r2 > 0 ? "--clip_r2 ${clip_r2}" : ''
  tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${three_prime_clip_r1}" : ''
  tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${three_prime_clip_r2}" : ''
  nextseq = params.trim_nextseq > 0 ? "--nextseq ${params.trim_nextseq}" : ''
  pair_str = paired == 'PE' ? "--paired" : ""
  clip_str = paired == 'PE' ? "$c_r1 $c_r2" : c_r1
  tpc_str = paired == 'PE' ? "$tpc_r1 $tpc_r2" : tpc_r1
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
  //conda "$NXF_CONDA_CACHEDIR/r"
  publishDir "${params.outdir}/50_final", mode:'copy', overwrite:'true'
  //publishDir "${params.qcdir}/${params.genome}/${params.name}", mode:'copy', overwrite:'true'

  input:
  path(design)
  path(mqc)
  val(paired)

  output:
  path("00.meta.tsv")

  script:
  """
  mv $design in.tsv
  $baseDir/bin/samplelist_addstat.R in.tsv $mqc --paired $paired 00.meta.tsv
  """
}

def get_reads(design) {
  reads = design
    .splitCsv(header:true, sep:"\t")
    .map { row -> [ row.SampleID,
            (params.source=='mixed') ? row.source : params.source,
            (params.paired=='mixed') ? row.paired : params.paired,
            row.r0, row.r1, row.r2
            ]
    }
  return reads
}

def get_reads2(design) {
  reads = Channel.empty()
  if( params.source == 'local' ) {
    reads_se = design
      .splitCsv(header:true, sep:"\t")
      .filter { it.paired != 'TRUE' }
      .map {row -> [row.SampleID, row.paired.toBoolean(),
        file(row.r0), file('f1', checkIfExists:false), file('f2', checkIfExists:false), has_ext(row.r0, ".gz"), null]}
    reads_pe = design
      .splitCsv(header:true, sep:"\t")
      .filter { it.paired == 'TRUE' }
      .map {row -> [row.SampleID, row.paired.toBoolean(),
        file('f0', checkIfExists:false), file(row.r1), file(row.r2), null, has_ext(row.r1, ".gz")]}
    reads = reads_se.concat(reads_pe)
  } else if( params.source == 's3' ) {
    reads = design
      .splitCsv(header:true, sep:"\t")
      .filter { it.paired != 'TRUE' }
      .map {row -> [row.SampleID, row.paired.toBoolean(), row.r0, row.r1, row.r2]}
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

workflow fq {
  take: design
  main:
    raw_read_list = get_reads(design)
    raw_read_list | fqu | fqc
    trim_reads = Channel.empty(); trim_log = Channel.empty(); trim_fqc = Channel.empty()
    if (params.trimmer == "no") {
      trim_reads = fqu.out
    } else if (params.trimmer == 'trim_galore'){
      fqu.out | trim_galore
      trim_reads = trim_galore.out.reads
      trim_log = trim_galore.out.log
      trim_fqc = trim_galore.out.fastqc_zip
    }
    mqc(fqc.out.zip.collect())
    upd(design, mqc.out, params.paired)
  emit:
    raw_reads = fqu.out
    trim_reads = trim_reads
    trim_log = trim_log
    raw_fqc = fqc.out.zip
    trim_fqc = trim_fqc
    readlist = upd.out
}



