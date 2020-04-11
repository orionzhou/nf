
process seqret {
  label 'low_memory'
  tag "$id"

  input:
  tuple id, path(lst), path(seqdb), path(seqdb_idx)

  output:
  tuple id, path("${id}.fas")

  script:
  """
  fasta.py extract --list $seqdb $lst > ${id}.fas
  """
}

process fimo {
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

process mg_fimo {
  label 'medium_memory'
  tag ""

  input:
  path(fis)

  output:
  path("fimo.tsv")

  script:
  """
  bioawk -t '{print substr(FILENAME,1,5), \$1, \$2, \$3}' $fis > fimo.tsv
  """
}


