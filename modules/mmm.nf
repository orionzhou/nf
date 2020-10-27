
process seqret {
  label 'low_memory'
  tag "$id"
  publishDir "${params.outdir}/21_seqs", mode:'copy', overwrite: true

  input:
  tuple val(id), path(lst), path(seqdb), path(seqdb_idx)

  output:
  tuple val(id), path("${id}.fas")

  script:
  """
  fasta.py extract --list $seqdb $lst > ${id}.fas
  """
}

process fimo_old {
  label 'low_memory'
  tag "$id"

  input:
  tuple val(id), path(seq), path(mtf), path(bg)

  output:
  tuple val(id), path("${id}.raw.tsv"), emit: raw
  tuple val(id), path("${id}.sum.tsv"), emit: sum

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
  tuple val(id), path(mtf), path(seq), path(seq_idx), path(bg)

  output:
  tuple val(id), path("${id}.tsv")

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
  tuple val(id), val(clid), path(seq), path(cseq)

  output:
  tuple val(id), path("${id}.txt")

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
  tuple val(id), val(clid), path(seq), path(cseq)

  output:
  tuple val(id), path("${id}.dreme"), path("${id}.tsv")

  script:
  mink = 6
  maxk = 13
  pval = 1e-4
  runtime = 180000
  runtime = 324000
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
  //conda "$NXF_CONDA_CACHEDIR/r"
  publishDir "${params.outdir}", mode:'copy', overwrite: true

  input:
  path(fis)

  output:
  path("fimo.rds")

  script:
  //bioawk -t '{print substr(FILENAME,1,length(FILENAME)-4), \$1, \$2, \$3, \$4, \$5}' $fis > fimo.tsv
  """
  $baseDir/bin/mmm/merge.fimo.R -o fimo.rds $fis
  """
}

process mg_dreme {
  label 'medium_memory'
  //conda "$NXF_CONDA_CACHEDIR/r"
  publishDir "${params.outdir}", mode:'copy', overwrite: true

  input:
  path(f_mtfs)
  path(f_tsvs)

  output:
  tuple path("23.dreme.rds"), path("23.dreme.meme"), path("23.dreme.txt"), path("23.kmer.tsv"), path("23.kmer.motif.tsv"), path("23.kmer.txt")

  script:
  """
  $baseDir/bin/mmm/merge.dreme.R -o 23.dreme.rds --meme 23.dreme.meme --txt 23.dreme.txt $f_mtfs
  $baseDir/bin/mmm/merge.dreme.kmer.R --fok 23.kmer.tsv --fom 23.kmer.motif.tsv --fos 23.kmer.txt $f_tsvs
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

process ml1 {
  label 'medium_memory'
  publishDir "${params.outdir}", mode:'copy', overwrite: true,
    saveAs: { fn ->
      if (fn.indexOf(".rds") > 0) "41_ml/$fn"
      else null
    }
  tag "${id}"

  input:
  tuple val(id), path(fi)

  output:
  tuple val(id), path("${id}.rds")

  script:
  alg = 'rf'
  fold = 10
  nlevel = 4
  perm = 10
  """
  ml_classification.R --perm $perm --alg $alg --fold $fold --nlevel $nlevel --downsample --seed $perm --cpu ${task.cpus} $fi ${id}.rds
  """
}

process mg_ml {
  label 'medium_memory'
  publishDir "${params.outdir}", mode:'copy', overwrite: true

  input:
  path(fis)

  output:
  path("42.ml.rds")

  script:
  """
  $baseDir/bin/merge.stats.R --opt rds -o 42.ml.rds $fis
  """
}


workflow ml {
  take:
    data_lst
  main:
    data_lst | ml1
    mg_ml(ml1.out.collect({it[1]}))
  emit:
    mg_ml = mg_ml.out
}
