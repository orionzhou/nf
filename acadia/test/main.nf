#!/usr/bin/env nextflow
nextflow.preview.dsl = 2

ch = Channel.from( [['r08','j'], ['r06','hj'], ['r05', 'hkk']])
ch2 = Channel.from( [['r04','j'], ['r05','hj'], ['r06', 'hkk']])

process mg {
  executor = 'local'
  input:
  val a

  output:
  path "00.raw.rds"

  script:
  b = a.sort(false)
  """
  echo $a
  echo $b
  touch 00.raw.rds
  """
}


workflow {
  main:
    x = ch2.join(ch, by:0, remainder: true)
      .toSortedList {entry -> entry[0]}
      .map {allPairs -> allPairs.collect{it[0]}}
    mg(x.collect())
}

chl = ch.toList()
workflow.onComplete {
  def nsam = chl.getVal().size()
  log.info "$nsam"
}


