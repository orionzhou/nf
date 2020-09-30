#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser(description = 'compile nextflow RNA-Seq result file')
parser$add_argument("meta", nargs=1, help="meta table (*.tsv)")
parser$add_argument("out", nargs=1, help="output file (*.rds)")
parser$add_argument("--bamstat", default='', help="bamstat")
parser$add_argument("--fcnt", default='', help="featurecounts")
parser$add_argument("--salmon_gcnt", default='', help="salmon gene counts")
parser$add_argument("--salmon_gtpm", default='', help="salmon gene TPM")
parser$add_argument("--salmon_tcnt", default='', help="salmon tx counts")
parser$add_argument("--salmon_ttpm", default='', help="salmon tx TPM")
parser$add_argument("--ase_gene", default='', help="gene ase")
parser$add_argument("--ase_snp", default='', help="snp ase")
parser$add_argument("--ril", default='', help="snpbinner (RIL) output")
args <- parser$parse_args()

require(tidyverse)

meta = args$meta
f_out = args$out

if( file.access(meta) == -1 )
    stop(sprintf("file ( %s ) cannot be accessed", meta))
th = read_tsv(meta)

bamstat = if (args$bamstat == '') NA else read_tsv(args$bamstat)
fcnt = if(args$fcnt == '') NA else readRDS(args$fcnt) %>% rename(SampleID=sid, ReadCount=cnt)

salmon = NA; salmon_tx = NA
salm = args$salmon_gcnt != '' && args$salmon_gtpm != '' && args$salmon_tcnt != '' && args$salmon_ttpm != ''
if (salm) {
    ti1 = readRDS(args$salmon_gcnt) %>% rename(ReadCount=val)
    ti2 = readRDS(args$salmon_gtpm) %>% rename(TPM=val)
    ti3 = readRDS(args$salmon_tcnt) %>% rename(ReadCount=val)
    ti4 = readRDS(args$salmon_ttpm) %>% rename(TPM=val)
    salmon = ti1 %>% inner_join(ti2, by=c("sid",'gid')) %>% rename(SampleID=sid)
    salmon_tx = ti3 %>% inner_join(ti4, by=c("sid",'gid')) %>%
        rename(SampleID=sid, tid=gid)
}

ase = args$ase_gene != '' && args$ase_snp != ''
ase_gene = NA; ase_snp = NA
if (ase) {
    ase_gene = readRDS(args$ase_gene) %>% rename(SampleID=sid)
    ase_snp = readRDS(args$ase_snp) %>% rename(SampleID=sid)
}
ril = if (args$ril == '') NA else readRDS(args$ril)

res = list(th=th, bamstat=bamstat, fcnt=fcnt,
           salmon = salmon, salmon_tx = salmon_tx,
           ase_gene = ase_gene, ase_snp = ase_snp, ril = ril)
saveRDS(res, file = f_out)

