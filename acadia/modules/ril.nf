process ril1 {
  label "process_medium"
  tag "${params.name}.$rid"

  input:
  tuple path(bams), path(bais), path(ref), path(sites), path(sites_idx), rid, region

  output:
  tuple rid, path("${rid}.vcf.gz"), path("${rid}.vcf.gz.tbi")

  when:
  params.ril

  script:
  """
  bcftools mpileup -f $ref -r $region -Ou $bams |\
    bcftools call -m -C alleles -T $sites -Oz -o ${rid}.vcf.gz
  bcftools index -t ${rid}.vcf.gz
  """
}

process ril2 {
  label "process_medium"
  tag "${params.name}"

  input:
  val rids
  path(vcfs)
  path(tbis)

  output:
  tuple path("ril.vcf.gz"), path("ril.vcf.gz.tbi")

  when:
  params.ril

  script:
  vcf_str = rids.sort().collect {"${it}.vcf.gz"}.join(' ')
  """
  echo $rids
  bcftools concat -n $vcf_str -Oz -o ril.vcf.gz
  bcftools index -t ril.vcf.gz
  """
}

process ril3 {
  label "process_medium"
  tag "${params.name}.$rid"
  conda "$NXF_CONDA_CACHEDIR/snpbinner"

  input:
  tuple path(vcf), path(tbi), rid, region

  output:
  path "${rid}.csv", emit: csv
  path "${rid}.txt", emit: txt

  when:
  params.ril

  script:
  """
  (bcftools query -l $vcf | tr "\\n" "\\t" |\
    sed "s/^/marker\\tposition\\t/; s/\\t\$/\\n/" &&
  bcftools query -i 'INFO/DP>=3' -r $region -f'%CHROM\\t%POS[\\t%GT]\\n' $vcf |\
    sed -e 's/\\.\\/\\./\\-/g; s/0\\/0/a/g; s/1\\/1/b/g; s/0\\/1/h/g' |\
    awk '{{OFS="\\t"}}; {{\$1 = \$1"_"\$2; print}}') > t1.tsv

  snpbinner crosspoints -i t1.tsv -o ${rid}.txt -r ${params.min_ratio}
  snpbinner bins -i ${rid}.txt -o ${rid}.csv -l ${params.min_bin_size}
  """
}



