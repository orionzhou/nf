
process seqret {
  label 'low_memory'
  tag "$id"
  publishDir "${params.outdir}/21_seqs", mode:'copy', overwrite: true

  input:
  tuple id, path(lst), path(seqdb), path(seqdb_idx)

  output:
  tuple id, path("${id}.fas")

  script:
  """
  fasta.py extract --list $seqdb $lst > ${id}.fas
  """
}

process fimo_old {
  label 'low_memory'
  tag "$id"

  input:
  tuple id, path(seq), path(mtf), path(bg)

  output:
  tuple id, path("${id}.raw.tsv"), emit: raw
  tuple id, path("${id}.sum.tsv"), emit: sum

  script:
  """
  fimo --skip-matched-sequence --text --bfile $bg $mtf $seq > ${id}.raw.tsv
  grep ^M ${id}.raw.tsv |\
    bioawk -tH '{if(\$7>0) {print \$1"%"\$3, \$4-1, \$5, ".", \$7, \$6}}' |\
    sortBed | mergeBed -s -c 6,5 -o distinct,max |\
    bioawk -tH '{split(\$1,a,/%/); print a[1], 1, 2, \$3-\$2, \$5}' |\
    mergeBed -c 4,5 -o sum,sum |\
    cut -f1,4,5 > ${id}.sum.tsv
  """
}

process fimo {
  label 'low_memory'
  tag "$id"
  publishDir "${params.outdir}/11_fimo_raw", mode:'copy', overwrite: true

  input:
  tuple id, path(mtf), path(seq), path(seq_idx), path(bg)

  output:
  tuple id, path("${id}.tsv")

  script:
  """
  fimo --skip-matched-sequence --text --motif $id --bfile $bg $mtf $seq > ${id}.tsv
  """
}

process fimo2 {
  label 'low_memory'
  tag "$id"
  publishDir "${params.outdir}/12_fimo_sum", mode:'copy', overwrite: true

  input:
  tuple id, path("in.tsv")

  output:
  tuple id, path("${id}.tsv")

  script:
  """
  parse_fimo.py in.tsv ${id}.tsv
  """
}

process meme {
  label 'low_memory'
  tag "$id"

  input:
  tuple id, clid, path(seq), path(cseq)

  output:
  tuple id, path("${id}.txt")

  script:
  """
  meme -objfun de -dna -mod zoops -revcomp -time 72000 \\
    -p ${task.cpus} -oc out $seq -neg $cseq
  mv out/meme.txt ${id}.txt
  """
}

process dreme {
  label 'low_memory'
  tag "$id"
  publishDir "${params.outdir}/22_dreme", mode:'copy', overwrite: true

  input:
  tuple id, clid, path(seq), path(cseq)

  output:
  tuple id, path("${id}.txt")

  script:
  mink = 5
  maxk = 20
  """
  dreme -p $seq -n $cseq -dna -e 0.05 -t 36000 -mink $mink -maxk $maxk -oc out
  mv out/dreme.txt ${id}.txt
  """
}

process mg_fimo {
  label 'medium_memory'
  publishDir "${params.outdir}", mode:'copy', overwrite: true

  input:
  path(fis)

  output:
  path("fimo.tsv")

  script:
  """
  bioawk -t '{print substr(FILENAME,1,length(FILENAME)-4), \$1, \$2, \$3, \$4, \$5}' $fis > fimo.tsv
  """
}

process mg_dreme {
  label 'medium_memory'
  conda '/home/springer/zhoux379/software/miniconda3/envs/r'
  publishDir "${params.outdir}", mode:'copy', overwrite: true

  input:
  path(fis)

  output:
  path("dreme.rds")

  script:
  """
  merge.dreme.R -o dreme.rds $fis
  """
}

process mg_meme {
  label 'medium_memory'
  conda '/home/springer/zhoux379/software/miniconda3/envs/r'

  input:
  path(fis)

  output:
  path("meme.rds")

  script:
  """
  merge.meme.R -o meme.rds $fis
  """
}

workflow mmk {
  take:
    seqdb
    mtfs
    mtf
    fimo_bg
  main:
    fimo(mtfs.combine(mtf).combine(seqdb).combine(fimo_bg))
    fimo2(fimo.out)
    mg_fimo(fimo2.out.collect({it[1]}))
  emit:
    fimo = fimo.out
    fimo2 = fimo2.out
    mg_fimo = mg_fimo.out
}

workflow mmd {
  take:
    seqdb
    lsts
    bg_lsts
    lst_pairs
  main:
    seqret(bg_lsts.concat(lsts).combine(seqdb))
    dreme_in = lst_pairs.combine(seqret.out, by: 0)
      .map { row -> [ row[1], row[0], row[2] ] }
      .combine(seqret.out, by: 0)
      .map { row -> [ row[1], row[0], row[2], row[3] ] }
    dreme(dreme_in)
    mg_dreme(dreme.out.collect({it[1]}))
    //meme(dreme_in)
    //mg_meme(meme.out.collect({it[1]}))
  emit:
    seqs = seqret.out
    dreme = dreme.out
    mg_dreme = mg_dreme.out
    //meme = meme.out
    //mg_meme = mg_meme.out
}


