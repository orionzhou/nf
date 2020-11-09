#!/usr/bin/env Rscript
suppressPackageStartupMessages(require(argparse))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(readr))
#suppressPackageStartupMessages(require(BBmisc))

parser <- ArgumentParser(description = 'wgc5 chain post processing')
parser$add_argument("fi", nargs=1, help="input chain BED")
parser$add_argument("fo", nargs=1, help="output chain BED")
parser$add_argument("fvi", nargs=1, help="input variant BED")
parser$add_argument("fvo", nargs=1, help="output variant BED")
parser$add_argument("--minsize", type="integer", default=50000,
                    help="output chain [default: %(default)s]")
args <- parser$parse_args()

fi = args$fi
fo = args$fo
minsize = args$minsize
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
              alnlen = sum(alnlen), block = n()) %>%
    filter(alnlen >= minsize)
to = ti %>% filter(cid %in% tis$cid)

cat(sprintf("%d out of %d chains passed filtering\n", nrow(tis), length(ti$cid)))
write_tsv(to, fo, col_names = F)

fvi = args$fvi
fvo = args$fvo
tv = read_tsv(fvi, col_names = F) %>%
    filter(X8 %in% tis$cid) %>%
    arrange(X1, X2, X3)
write_tsv(tv, fvo, col_names = F)

