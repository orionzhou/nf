process hs2 {
  label 'high_memory'
  tag "$id"
  publishDir "${params.outdir}/11_hisat2", mode:'copy', overwrite:'true',
    saveAs: { fn ->
      if (fn.indexOf(".hisat2_summary.txt") > 0) "logs/$fn"
      else if (!params.saveBAM && fn == 'where_are_my_files.txt') fn
      else if (params.saveBAM && fn != 'where_are_my_files.txt') fn
      else null
    }

  input:
  tuple val(id), val(paired), path(reads), val(pre), path(hs2_indices)

  output:
  tuple val(id), path("${id}.bam"), emit: bam
  path "${id}.hisat2_summary.txt", emit: log
  path "unmapped.hisat2*" optional true

  script:
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
  extra = params.cage ? "--sp 1,0.1 --score-min L,0,-2" : ""
  // unaligned = params.save_unmapped ? !paired ? "--un-gz unmapped.hisat2.gz" : "--un-conc-gz unmapped.hisat2.gz" : ''
  // --dta --no-mixed --no-discordant
  // --known-splicesite-infile $splicesites \\
  // ${unaligned} \\
  """
  hisat2 -x ${params.hisat2_index} \\
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
  tuple val(id), val(paired), path(reads), path(index), path(gtf)

  output:
  tuple val(id), path('*.bam'), emit: bam
  path "*.out", emit: log
  path "*SJ.out.tab"
  path "*Log.out"
  path "*Unmapped*" optional true

  script:
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
  STAR --genomeDir $index \\
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
  tuple val(id), val(paired), path(reads), val(pre), path(index)

  output:
  tuple val(id), path("${id}.bam"), emit: bam
  path "${id}.log", emit: log optional true

  script:
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
  tuple val(id), val(paired), path(reads), val(pre), path(index)

  output:
  tuple val(id), path("${id}.bam"), emit: bam
  path "${id}.log", emit: log optional true

  script:
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
  tuple val(id), val(paired), path(reads), val(pre), path(index)

  output:
  tuple val(id), path("${id}.bam"), emit: bam
  path "${id}.txt", emit: log

  script:
  extra = "-n 1 --rg_tag --rg_id ${id} --rg_sample ${id}"
  """
  mkdir tmp
  bismark --multicore ${task.cpus.intdiv(2)} \\
    $pre ${reads} \\
    --temp_dir tmp \\
    ${extra} \\
    --output_dir .
  mv ${id}*bismark*.bam ${id}.bam
  mv ${id}*bismark*_report.txt ${id}.txt
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
  tuple val(id), val(paired), path(reads), path(index)

  output:
  tuple val(id), path("${id}.bam"), emit: bam
  path "${id}.log", emit: log optional true

  script:
  rg = "'@RG\\tID:${id}\\tSM:${id}\\tPL:ILLUMINA'"
  //minimap2 -ax splice -uf -k14 -G 90000 \\
  """
  minimap2 -ax splice -k14 -G 90000 \\
    ${index} ${reads} -t ${task.cpus} -R $rg |\\
    samtools view -@ ${task.cpus} -bSh -O BAM -o ${id}.bam -
  """
}

workflow aln {
  take: reads
  main:
    aln = Channel.empty()
    if (params.aligner == 'hisat2') {
      hs2_indices = Channel
        .fromPath("${params.hisat2_index}*ht2*", checkIfExists: true)
        .ifEmpty { exit 1, "HISAT2 index not found: ${params.hisat2_index}" }
        .collect()
        .map {r -> [params.hisat2_index, r]}
        .map {r -> [r[0].substring(r[0].lastIndexOf(File.separator)+1), r[1]]}
      //splicesites = Channel
        //.fromPath(params.bed12, checkIfExists: true)
        //.ifEmpty { exit 1, "HISAT2 splice sites file not found: $splicesites" }
      hs2(reads.combine(hs2_indices))
      aln = hs2.out
    } else if (params.aligner == 'star') {
      gff = Channel.fromPath(params.gff, checkIfExists: true)
         .ifEmpty { exit 1, "GFF annotation file not found: ${params.gff}" }
      gtf = Channel.fromPath(params.gtf, checkIfExists: true)
         .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
      star_index = Channel
        .fromPath(params.star_index.replaceAll("\\/+\$", ""), checkIfExists: true)
        .ifEmpty { exit 1, "STAR index not found: ${params.star_index}" }
      star(reads.combine(star_index).combine(gtf))
      aln = star.out
    } else if (params.aligner in ['bwa','bwa_aln']) {
      lastPath = params.bwa_index.lastIndexOf(File.separator)
      bwa_dir =  params.bwa_index.substring(0,lastPath+1)
      bwa_base = params.bwa_index.substring(lastPath+1)
      bwa_index = Channel
        .fromPath("${bwa_dir}*", checkIfExists: true)
        .ifEmpty { exit 1, "bwa index not found: ${params.bwa_index}" }
        .collect()
        .map {r -> [params.bwa_index, r]}
        .map {r -> [r[0].substring(r[0].lastIndexOf(File.separator)+1), r[1]]}
      if (params.aligner == 'bwa_aln') {
        bwa_aln(reads.combine(bwa_index))
        aln = bwa_aln.out
      } else {
        bwa(reads.combine(bwa_index))
        aln = bwa.out
      }
    } else if (params.aligner == 'minimap2') {
      genome_fasta = Channel.fromPath(params.fasta, checkIfExists: true)
        .ifEmpty { exit 1, "Genome fasta file not found: ${params.fasta}" }
      minimap2(reads.combine(genome_fasta))
      aln = minimap2.out
    }
  emit:
    bam = aln.bam
    l = aln.log
}





