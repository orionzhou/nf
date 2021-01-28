
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

process getfasta1 {
  label 'low_memory'
  tag "$lid"
  publishDir "${params.dm_dir}/${params.dm_tag}/20_seqs", mode:'copy'

  input:
  tuple val(lid), val(bin), val(epi), path(fg)

  output:
  tuple val(lid), path("${lid}.fas") optional true

  script:
  """
  $baseDir/bin/kmer.py getfasta --bin '$bin' --epi $epi $fg ${lid}.fas
  """
}
process getfasta2 {
  label 'low_memory'
  tag "$lid"
  publishDir "${params.dm_dir}/${params.dm_tag}/20_seqs", mode:'copy'

  input:
  tuple val(lid), val(bin), val(epi), path(fg)

  output:
  tuple val(lid), path("${lid}.fas") optional true

  script:
  """
  $baseDir/bin/kmer.py getfasta --bin '$bin' --epi $epi $fg ${lid}.fas
  """
}

process nseq {
  label 'low_memory'
  publishDir "${params.dm_dir}/${params.dm_tag}", mode:'copy',
    saveAs: { fn ->
      if (fn.indexOf("nseq.txt") >= 0) "20.nseq.txt"
      else null
    }

  input:
  path(seqs)

  output:
  path("nseq.txt")

  script:
  """
  grep -c '>' $seqs > nseq.txt
  """
}

process dreme {
  label 'low_memory'
  tag "$id"
  publishDir "${params.dm_dir}/${params.dm_tag}/21_dreme_raw", mode:'copy', overwrite: true

  input:
  tuple val(id), val(ctrl), path(seq), path(cseq)

  output:
  tuple val(id), path("${id}")

  script:
  mink = 6
  maxk = 13
  minw = 6
  maxw = 20
  pval = 1e-2
  runtime = 324000
  runtime = 180000
  runtime = 72000
  //fimo --bfile --motif-- --thresh $pval2 ${id}.dreme $seq || (mkdir fimo_out; touch fimo_out/fimo.tsv)
  //mv fimo_out/fimo.tsv ${id}.tsv
  //dreme -p $seq -n $cseq -dna -e $pval -t $runtime -mink $mink -maxk $maxk -oc out
  //mv out/dreme.txt ${id}.dreme
  """
  streme --p $seq --n $cseq --dna --pvt $pval --time $runtime -minw $minw -maxw $maxw -oc $id
  """
}

process dreme2 {
  label 'low_memory'
  tag "$id"
  publishDir "${params.dm_dir}/${params.dm_tag}", mode:'copy', overwrite: true,
    saveAs: { fn ->
      if (fn.indexOf(".dreme") > 0) "22_motifs/$fn"
      else if (fn.indexOf(".meme") > 0) "22_motifs/$fn"
      else if (fn.indexOf(".bed") > 0) "22_fimo/${fn.replaceAll('.sum','')}"
      else null
    }

  input:
  tuple val(id), val(ctrl), path(seq), path(cseq), path("in")

  output:
  tuple val(id), path("${id}.dreme"), path("${id}.tsv"), path("${id}.bed")

  script:
  """
  $baseDir/bin/mmm/streme.py xml2tsv in/streme.xml > ${id}.tsv
  $baseDir/bin/mmm/streme.py addscore in/streme.txt in/streme.xml ${id}.dreme
  cat $seq $cseq > merge1.fas
  $baseDir/bin/mmm/fimo.py locate ${id}.dreme merge1.fas ${id}.bed
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
  publishDir "${params.dm_dir}/${params.dm_tag}", mode:'copy', overwrite: true,
    saveAs: { fn ->
      if (fn.indexOf("dreme.rds") >= 0) "23.dreme.rds"
      else if (fn.indexOf("dreme.meme") >= 0) "23.dreme.meme"
      else if (fn.indexOf("fimo.rds") >= 0) "23.fimo.rds"
      else if (fn.indexOf("kmer.tsv") >= 0) "23.kmer.tsv"
      else if (fn.indexOf("kmer.motif.tsv") >= 0) "23.kmer.motif.tsv"
      else null
    }

  input:
  path(f_mtfs)
  path(f_beds)

  output:
  tuple path("dreme.rds"), path("fimo.rds")

  script:
  """
  $baseDir/bin/mmm/merge.dreme.R -o dreme.rds --meme dreme.meme --txt dreme.txt $f_mtfs
  $baseDir/bin/mmm/merge.dreme.fimo.R -o fimo.rds $f_beds
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

workflow dm {
  take:
    dm_cfg
  main:
    // setup channels
      qrys = dm_cfg.map{r -> [r.lid, r.bin, r.epi,
        file("${params.dm_dir}/${params.dm_tag}/02_gene_lists/${r.cid}.txt", checkIfExists:true)
      ]}.unique()
      tgts = dm_cfg.map{r -> [r.clid, r.bin, r.epi,
        file("${params.dm_dir}/${params.dm_tag}/02_gene_lists/${r.ccid}.txt", checkIfExists:true)
      ]}.unique()
    qrys | getfasta1
    tgts | getfasta2
    getfasta1.out.concat(getfasta2.out).collect({it[1]}) | nseq
    dreme_in = dm_cfg
      .map {r -> [r.lid, r.clid]}
      .combine(getfasta1.out, by:0)
      .combine(getfasta2.out.map {r -> [r[1], r[0]]}, by:1)
      .map {r -> [r[1],r[0],r[2],r[3]]}
    dreme_in | dreme
    dreme_in.combine(dreme.out, by:0) | dreme2
    mg_dreme(dreme2.out.collect({it[1]}), dreme2.out.collect({it[3]}))
  emit:
    dreme = dreme.out
    mg_dreme = mg_dreme.out
}

process ml1 {
  label 'low_memory'
  publishDir "${params.ml_dir}/${params.ml_tag}/41_ml_in", mode:'copy', overwrite: true
  tag "${id}"

  input:
  tuple val(id), val(bin), val(epi), val(nfea), val(mod), path("module.tsv"), path("mtfs.meme")

  output:
  tuple val(id), path("${id}.tsv")

  script:
  //$baseDir/bin/kmer.py prepare_ml --bin '$bin' --epi $epi --nfea $nfea --mod $mod module.tsv mtf.tsv ${id}.tsv
  """
  $baseDir/bin/mmm/fimo.py prepare_ml --bin '$bin' --epi $epi --nfea $nfea --mod $mod --fmt wide module.tsv mtfs.meme ${id}.tsv
  """
}

process ml2 {
  label 'medium_memory'
  publishDir "${params.ml_dir}/${params.ml_tag}/42_ml_out", mode:'copy', overwrite: true
  tag "${id}"

  input:
  tuple val(id), path(fi)

  output:
  tuple val(id), path("${id}.rds")

  script:
  alg = 'rf'
  fold = 10
  nlevel = 3
  perm = 20
  """
  $baseDir/bin/mmm/ml_classification.R --perm $perm --alg $alg --fold $fold --nlevel $nlevel --downsample --seed $perm --cpu ${task.cpus} $fi ${id}.rds
  """
}

process mg_ml {
  label 'high_memory'
  publishDir "${params.ml_dir}/${params.ml_tag}", mode:'copy', overwrite: true
  publishDir "${params.ml_dir}", mode:'copy', overwrite: true,
    saveAs: { fn ->
      if (fn.indexOf(".rds") > 0) "${params.ml_tag}.rds"
      else null
    }

  input:
  path(fis)

  output:
  path("43.ml.rds")

  script:
  """
  $baseDir/bin/merge.stats.R --opt rds -o 43.ml.rds $fis
  """
}


workflow ml {
  take:
    ml_cfg
  main:
    ml_cfg | ml1 | ml2
    mg_ml(ml2.out.collect({it[1]}))
  emit:
    mg_ml = mg_ml.out
}
