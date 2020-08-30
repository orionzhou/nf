
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
  publishDir "${params.outdir}", mode:'copy', overwrite: true,
    saveAs: { fn ->
      if (fn.indexOf(".tsv") > 0) "11_fimo/${fn.replaceAll('.raw','')}"
      else null
    }

  input:
  tuple id, path(mtf), path(seq), path(seq_idx), path(bg)

  output:
  tuple id, path("${id}.tsv")

  script:
  pval = 1e-4
  """
  fimo --motif $id --bfile $bg --thresh $pval $mtf $seq
  mv fimo_out/fimo.tsv ${id}.tsv
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
  publishDir "${params.outdir}", mode:'copy', overwrite: true,
    saveAs: { fn ->
      if (fn.indexOf(".dreme") > 0) "22a_dreme_mtf/$fn"
      else if (fn.indexOf(".tsv") > 0) "22b_dreme_kmer/${fn.replaceAll('.sum','')}"
      else null
    }

  input:
  tuple id, clid, path(seq), path(cseq)

  output:
  tuple id, path("${id}.dreme"), path("${id}.tsv")

  script:
  mink = 6
  maxk = 13
  pval = 1e-4
  runtime = 180000
  //fimo --bfile --motif-- --thresh $pval2 ${id}.dreme $seq || (mkdir fimo_out; touch fimo_out/fimo.tsv)
  //mv fimo_out/fimo.tsv ${id}.tsv
  """
  dreme -p $seq -n $cseq -dna -e $pval -t $runtime -mink $mink -maxk $maxk -oc out
  mv out/dreme.txt ${id}.dreme
  dreme.py 2tsv ${id}.dreme >${id}.tsv
  """
}

process mg_fimo {
  label 'medium_memory'
  conda '/home/springer/zhoux379/software/miniconda3/envs/r'
  publishDir "${params.outdir}", mode:'copy', overwrite: true

  input:
  path(fis)

  output:
  path("fimo.rds")

  script:
  //bioawk -t '{print substr(FILENAME,1,length(FILENAME)-4), \$1, \$2, \$3, \$4, \$5}' $fis > fimo.tsv
  """
  merge.fimo.R -o fimo.rds $fis
  """
}

process mg_dreme {
  label 'medium_memory'
  conda '/home/springer/zhoux379/software/miniconda3/envs/r'
  publishDir "${params.outdir}", mode:'copy', overwrite: true

  input:
  path(f_mtfs)
  path(f_tsvs)

  output:
  tuple path("dreme.rds"), path("dreme.meme"), path("dreme.txt"), path("dreme_kmer.tsv")

  script:
  """
  merge.dreme.R -o dreme.rds --meme dreme.meme --txt dreme.txt $f_mtfs
  merge.dreme.kmer.R -o dreme_kmer.tsv $f_tsvs
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
    mg_fimo(fimo.out.collect({it[1]}))
  emit:
    fimo = fimo.out
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
    mg_dreme(dreme.out.collect({it[1]}), dreme.out.collect({it[2]}))
  emit:
    seqs = seqret.out
    dreme = dreme.out
    mg_dreme = mg_dreme.out
}


