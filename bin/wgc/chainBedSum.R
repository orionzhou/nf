#!/usr/bin/env Rscript
suppressPackageStartupMessages(require(argparse))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(readr))

parser <- ArgumentParser(description = 'chain BED summary')
parser$add_argument("fi", nargs=1, help="input chain BED")
parser$add_argument("fo", nargs=1, help="output summary tsv")
args <- parser$parse_args()

fi = args$fi
fo = args$fo
if( file.access(fi) == -1 )
    stop(sprintf("Input file ( %s ) cannot be accessed", fi))

cols8 = c('tchrom','tstart','tend','srd','qchrom','qstart','qend','cid')
cols9 = c(cols8, 'mm')
ti = read_tsv(fi, col_names = F)
if(ncol(ti) == 8) {
    colnames(ti) = cols8
} else if(ncol(ti) == 9) {
    colnames(ti) = cols9
}

tis = ti %>% mutate(alnlen = tend - tstart) %>%
    group_by(cid) %>%
    summarise(tchrom=tchrom[1], tstart=min(tstart), tend=max(tend),
              qchrom=qchrom[1], qstart=min(qstart), qend=max(qend),
              srd=srd[1],
              tsize = tend-tstart, qsize = qend - qstart,
              alnlen = sum(alnlen), n_block = n(), mismatch = sum(mm)) %>%
    ungroup() %>%
    mutate(tstart = tstart + 1, qstart = qstart + 1)

cat(sprintf("%d chains read\n", nrow(tis)))
to = tis
write_tsv(to, fo, col_names = T)


