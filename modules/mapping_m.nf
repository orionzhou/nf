process hs2 {
  label 'high_memory'
  tag "$id"
  publishDir "${params.outdir}/11_hisat2", mode:'copy', overwrite:'true',
    saveAs: { fn ->
      if (fn.indexOf(".hisat2_summary.txt") > 0) "logs/$fn"
      //else if (!params.saveBAM && fn == 'where_are_my_files.txt') fn
      //else if (params.saveBAM && fn != 'where_are_my_files.txt') fn
      else null
    }

  input:
  tuple val(sid), val(paired), path(reads), val(genome), val(pre), path(hs2_indices)

  output:
  tuple val(id), path("${id}.bam"), emit: bam
  path "${id}.hisat2_summary.txt", emit: log
  path "unmapped.hisat2*" optional true

  script:
  id = "${sid}-${genome}"
  seq_center = params.seq_center ? 'CN:$params.seq_center' : ''
  rg = "--rg-id ${id} --rg ${seq_center.replaceAll('\\s','_')} SM:${id}"
  opt_splice = params.lib == 'chipseq' ? '--no-spliced-alignment' : ''
  def strandness = ''
  if (params.lib == 'rnaseq' && params.stranded == 'forward') {
    strandness = paired == 'SE' ? '--rna-strandness F' : '--rna-strandness FR'
  } else if (params.lib == 'rnaseq' && params.stranded == 'reverse') {
    strandness = paired == 'SE' ? '--rna-strandness R' : '--rna-strandness RF'
  }
  input = paired == 'PE' ? "-1 ${reads[0]} -2 ${reads[1]}" : "-U ${reads}"
  extra = ""//params.cage ? "--sp 1,0.1 --score-min L,0,-2" : ""
  // unaligned = params.save_unmapped ? !paired ? "--un-gz unmapped.hisat2.gz" : "--un-conc-gz unmapped.hisat2.gz" : ''
  // --dta --no-mixed --no-discordant
  // --known-splicesite-infile $splicesites \\
  // ${unaligned} \\
  """
  hisat2 -x ${pre} \\
     ${input} $opt_splice $strandness $extra \\
     -p ${task.cpus} --met-stderr --new-summary \\
     --summary-file ${id}.hisat2_summary.txt $rg \\
     | samtools view -bSh -o ${id}.bam -
  """
}

process star {
  label 'high_memory'
  tag "$id"
  publishDir "${params.outdir}/11_star", mode:'copy', overwrite:'true',
    saveAs: { fn ->
      if (fn.indexOf(".bam") == -1) "logs/$fn"
      else if (params.save_unmapped && fn != "where_are_my_files.txt" && 'Unmapped' in fn) "unmapped/$fn"
      else if (!params.saveBAM && fn == 'where_are_my_files.txt') fn
      else if (params.saveBAM && fn != 'where_are_my_files.txt') fn
      else null
    }

  input:
  tuple val(sid), val(paired), path(reads), val(genome), val(pre), path(index), path(gtf)

  output:
  tuple val(id), path('*.bam'), emit: bam
  path "*.out", emit: log
  path "*SJ.out.tab"
  path "*Log.out"
  path "*Unmapped*" optional true

  script:
  id = "${sid}-${genome}"
  def avail_mem = "--limitBAMsortRAM ${task.memory.toBytes() - 100000000}"
  seq_center = params.seq_center ? 'CN:$params.seq_center' : ''
  unaligned = params.save_unmapped ? "--outReadsUnmapped Fastx" : ''
  extra = """\\
    --outFilterType BySJout \\
    --outSAMmapqUnique 60 \\
    --outFilterMultimapNmax 10 \\
    --outFilterMismatchNmax 999 \\
    --outFilterMismatchNoverLmax 1 \\
    --outFilterMismatchNoverReadLmax 1.0 \\
    --outFilterMatchNminOverLread 0 \\
    --outFilterMatchNmin 0 \\
    --outFilterScoreMinOverLread 0 \\
    --alignSJoverhangMin 8 \\
    --alignSJDBoverhangMin 1 \\
    --alignIntronMin 20 \\
    --alignIntronMax 1000000 \\
    --alignMatesGapMax 1000000 \\
    --outSAMunmapped Within KeepPairs \\
    --outSAMattrRGline ID:$id $seq_center 'SM:$id'
    """.stripIndent()
  """
  STAR --genomeDir $pre \\
    --sjdbGTFfile $gtf \\
    --runThreadN ${task.cpus} \\
    --readFilesIn $reads  \\
    --readFilesCommand zcat \\
    --twopassMode Basic \\
    --outSAMtype BAM Unsorted $avail_mem \\
    --outFileNamePrefix ${id}_ \\
    --runDirPerm All_RWX $unaligned \\
    --outStd Log $extra
  mv ${id}_Aligned.out.bam ${id}.bam
  """
    //--outWigType bedGraph \\
}

process bwa {
  label 'high_memory'
  tag "$id"
  publishDir "${params.outdir}/11_bwa", mode:'copy', overwrite:'true',
    saveAs: { fn ->
      if (fn.indexOf(".bam") == -1) "logs/$fn"
      else if (!params.saveBAM && fn == 'where_are_my_files.txt') fn
      else if (params.saveBAM && fn != 'where_are_my_files.txt') fn
      else null
    }

  input:
  tuple val(sid), val(paired), path(reads), val(genome), val(pre), path(index)

  output:
  tuple val(id), path("${id}.bam"), emit: bam
  path "${id}.log", emit: log optional true

  script:
  id = "${sid}-${genome}"
  rg = "'@RG\\tID:${id}\\tSM:${id}\\tPL:ILLUMINA'"
  """
  bwa mem -t ${task.cpus} $pre \\
    -R $rg $reads |\\
    samtools view -@ ${task.cpus} -b -h -O BAM -o ${id}.bam -
  """
}

process bwa_aln {
  label 'high_memory'
  tag "$id"
  publishDir "${params.outdir}/11_bwa_aln", mode:'copy', overwrite:'true',
    saveAs: { fn ->
      if (fn.indexOf(".bam") == -1) "logs/$fn"
      else if (!params.saveBAM && fn == 'where_are_my_files.txt') fn
      else if (params.saveBAM && fn != 'where_are_my_files.txt') fn
      else null
    }

  input:
  tuple val(sid), val(paired), path(reads), val(genome), val(pre), path(index)

  output:
  tuple val(id), val(genome), path("${id}.bam"), emit: bam
  path "${id}.log", emit: log optional true

  script:
  id = "${sid}-${genome}"
  rg = "'@RG\\tID:${id}\\tSM:${id}\\tPL:ILLUMINA'"
  if( paired == 'PE')
    """
    bwa aln -t ${task.cpus} $pre $reads[0] > out1.sai
    bwa aln -t ${task.cpus} $pre $reads[1] > out2.sai
    bwa sampe -r $rg $pre out1.sai out2.sai $reads |\\
      samtools view -@ ${task.cpus} -b -h -O BAM -o ${id}.bam -
    """
  else
    """
    bwa aln -t ${task.cpus} $pre $reads > out.sai
    bwa samse -r $rg $pre out.sai $reads |\\
      samtools view -@ ${task.cpus} -b -h -O BAM -o ${id}.bam -
    """
}

process bsmk {
  label 'high_memory'
  tag "$id"
  publishDir "${params.outdir}/11_bismark", mode:'copy', overwrite:'true',
    saveAs: { fn ->
      if (fn.indexOf(".bam") == -1) "logs/$fn"
      else if (!params.saveBAM && fn == 'where_are_my_files.txt') fn
      else if (params.saveBAM && fn != 'where_are_my_files.txt') fn
      else null
    }

  input:
  tuple val(sid), val(paired), path(reads), val(genome), val(pre), path(index)

  output:
  tuple val(id), path("${id}.bam"), emit: bam
  path "${id}.txt", emit: log

  script:
  id = "${sid}-${genome}"
  extra = "-n 1 --rg_tag --rg_id ${id} --rg_sample ${id}"
  """
  mkdir tmp
  bismark --multicore ${task.cpus.intdiv(2)} \\
    $pre ${reads} \\
    --temp_dir tmp \\
    ${extra} \\
    --output_dir .
  mv ${sid}*bismark*.bam ${id}.bam
  mv ${sid}*bismark*_report.txt ${id}.txt
  """
}

process minimap2 {
  label 'high_memory'
  tag "$id"
  publishDir "${params.outdir}/11_minimap2", mode:'copy', overwrite:'true',
    saveAs: { fn ->
      if (fn.indexOf(".bam") == -1) "logs/$fn"
      else if (!params.saveBAM && fn == 'where_are_my_files.txt') fn
      else if (params.saveBAM && fn != 'where_are_my_files.txt') fn
      else null
    }

  input:
  tuple val(sid), val(paired), path(reads), val(genome), val(pre), path(index)

  output:
  tuple val(id), path("${id}.bam"), emit: bam
  path "${id}.log", emit: log optional true

  script:
  id = "${sid}-${genome}"
  rg = "'@RG\\tID:${id}\\tSM:${id}\\tPL:ILLUMINA'"
  //minimap2 -ax splice -uf -k14 -G 90000 \\
  """
  minimap2 -ax splice -k14 -G 90000 \\
    ${pre} ${reads} -t ${task.cpus} -R $rg |\\
    samtools view -@ ${task.cpus} -bSh -O BAM -o ${id}.bam -
  """
}


workflow aln {
  take:
    reads
    genomes
  main:
    aln = Channel.empty()
    if (params.aligner == 'hisat2') {
    // [genome, db_pre, db_files]
      hs2_indices = genomes
        .map {r -> [r[0], r[1].hisat2]}
        .map {r -> [r[0], r[1].substring(r[1].lastIndexOf(File.separator)+1),
          file("${r[1]}*ht2*").toList()]}
      hs2(reads.combine(hs2_indices))
      aln = hs2.out
    } else if (params.aligner == 'star') {
      gff = genomes
        .map {r -> [r[0], file(r[1].gff, checkIfExists: true)]}
      gtf = genomes
        .map {r -> [r[0], file(r[1].gtf, checkIfExists: true)]}
      ch_index = genomes
        .map {r -> [r[0], r[1].star.replaceAll("\\/+\$", ''),
          file(r[1].gff, checkIfExists: true)]}
        .map {r -> [r[0], r[1].substring(r[1].lastIndexOf(File.separator)+1),
          file(r[1], checkIfExists: true), r[2]
          ]}
      star(reads.combine(ch_index))
      aln = star.out
    } else if (params.aligner in ['bwa','bwa_aln']) {
      bwa_index = genomes
        .map {r -> [r[0], r[1].bwa]}
        .map {r -> [r[0], r[1].substring(r[1].lastIndexOf(File.separator)+1),
          file("${r[1]}*").toList()]}
      if (params.aligner == 'bwa_aln') {
        bwa_aln(reads.combine(bwa_index))
        aln = bwa_aln.out
      } else {
        bwa(reads.combine(bwa_index))
        aln = bwa.out
      }
    } else if (params.aligner == 'minimap2') {
      ch_fasta = genomes
        .map {r -> [r[0], r[1].fasta]}
        .map {r -> [r[0], r[1].substring(r[1].lastIndexOf(File.separator)+1),
          file(r[1], checkIfExists: true)]}
      minimap2(reads.combine(ch_fasta))
      aln = minimap2.out
    }
  emit:
    bam = aln.bam
    l = aln.log
}
