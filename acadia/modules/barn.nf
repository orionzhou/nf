
process fqd {
  label 'low_memory'
  tag "$id"

  input:
  tuple id, paired, interleaved, path(r0), path(r1), path(r2)

  output:
  tuple id, paired, path(reads)

  script:
  def mem = task.memory
  def reads = []
  """
  fasterq-dump --split-files -e {task.cpus} -m {task.memory} \
      -O ./ -t $TMPDIR ${id}
  """
  if( paired ) {
    """
    pigz -p {task.cpus} --fast -c ${id}_1.fastq > ${id}_R1.fq.gz
    pigz -p {task.cpus} --fast -c ${id}_2.fastq > ${id}_R2.fq.gz
    """
    reads = ["${id}_R1.fq.gz", "${id}_R2.fq.gz"]
  } else {
    """
    pigz -p {task.cpus} --fast -c ${id}.fastq > ${id}_R0.fq.gz
    """
    reads = ["${id}_R0.fq.gz"]
  }
}

process fqz {
  label 'low_memory'
  tag "$id"

  input:
  tuple id, paired, interleaved, path(r0), path(r1), path(r2)

  output:
  tuple id, paired, interleaved, path("${id}_R0.fq.gz"), path("${id}_R1.fq.gz"), path("${id}_R2.fq.gz")

  script:
  if( paired ) {
    if( hasExtention(r1, ".gz") ) {
      try {
        """
        ln -f $r1 ${id}_R1.fq.gz
        ln -f $r2 ${id}_R2.fq.gz
        """
      } catch (all) {
        """
        cp -fL $r1 ${id}_R1.fq.gz
        cp -fL $r2 ${id}_R2.fq.gz
        """
      }
    } else {
      """
      pigz -p {task.cpus} -c $r1 > ${id}_R1.fq.gz
      pigz -p {task.cpus} -c $r2 > ${id}_R2.fq.gz
      """
    }
    """
    touch ${id}_R0.fq.gz
    """
  } else {
    if( hasExtention(r1, ".gz") ) {
      try {
        """
        ln -f $r0 ${id}_R0.fq.gz
        """
      } catch (all) {
        """
        cp -fL $r0 ${id}_R0.fq.gz
        """
      }
    } else {
      """
      pigz -p {task.cpus} -c $r0 > ${id}_R0.fq.gz
      """
    }
    """
    touch ${id}_R1.fq.gz {$id}_R1.fq.gz
    """
  }
}

process fqv {
  label 'low_memory'
  tag "$id"

  input:
  tuple id, paired, interleaved, path(r0), path(r1), path(r2)

  output:
  tuple id, paired, path(reads)

  script:
  def reads = []
  if( interleaved ) {
    """
    zcat $r0 |\
      deinterleave_fastq.sh $r1 $r2 {task.cpus} compress
    """
  }
  if( paired )
    reads = [r1, r2]
  else
    reads = [r0]
}

process fqc {
  label 'low_memory'
  tag "$id"

  input:
  tuple id, paired, path(reads)

  output:
  path(zips)

  script:
  """
  fastqc --threads {task.cpus} --noextract --format fastq -o . ${reads}
  """
  def zips = []
  if( paired ) {
    zips = ["${id}_R1_fastqc.zip", "${id}_R2_fastqc.zip"]
  } else
    zips = ["${id}_R0_fastqc.zip"]
  }
}

process mqc {
  label 'low_memory'
  tag ""

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
  tag ""
  conda '/home/springer/zhoux379/software/miniconda3/envs/r'

  input:
  path(design)
  path(mqc)

  output:
  path("${params.name}.tsv")

  script:
  """
  samplelist_addstat.R $design $mqc ${params.name}.tsv
  """
}

