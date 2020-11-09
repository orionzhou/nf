#!/usr/bin/env Rscript
suppressPackageStartupMessages(require(argparse))

parser <- ArgumentParser(description = 'prepare for vntcall from 8-col chainBed file')
parser$add_argument("fi", nargs=1, help="input chain BED")
parser$add_argument("fo", nargs=1, help="output interval BED")
parser$add_argument("--maxsize", type="integer", default=0,
                    help="max indel size [default: %(default)s]")
args <- parser$parse_args()

require(dplyr)
require(readr)
fi = args$fi
fo = args$fo
maxsize = ifelse(args$maxsize == 0, Inf, args$maxsize)
if( file.access(fi) == -1 )
    stop(sprintf("Input file ( %s ) cannot be accessed", fi))

cols8 = c('tchrom','tstart','tend','srd','qchrom','qstart','qend','cid')
cols9 = c(cols8, 'mm')
ti = read_tsv(fi, col_names = F)
if(ncol(ti) == 8) {
    colnames(ti) = cols8
} else if(ncol(ti) == 9) {
    colnames(ti) = cols9
    ti = ti[,cols8]
}
tc = ti %>% mutate(alnlen = tend - tstart) %>%
    group_by(cid) %>%
    summarise(tchrom=tchrom[1], tstart=min(tstart), tend=max(tend),
              qchrom=qchrom[1], qstart=min(qstart), qend=max(qend),
              srd=srd[1], alnlen = sum(alnlen), block = n())

ta = ti %>% arrange(cid, tchrom, tstart)
ta1 = ta[-nrow(ta),] %>%
    transmute(cid1 = cid, tchrom = tchrom, qchrom = qchrom, srd = srd, ts1 = tstart, te1 = tend, qs1 = qstart, qe1 = qend)
ta2 = ta[-1,] %>%
    transmute(cid2 = cid, ts2 = tstart, te2 = tend, qs2 = qstart, qe2 = qend)
ta3 = ta1 %>% bind_cols(ta2) %>%
    filter(cid1 == cid2)
stopifnot(nrow(ta3) + nrow(tc) == nrow(ta))
ta4 = ta3 %>% mutate(cid = cid1) %>% select(-cid1,-cid2) %>%
    mutate(tstart = te1, tend = ts2,
           qstart = ifelse(srd == '-', qe2, qe1),
           qend = ifelse(srd == '-', qs1, qs2),
           td = tend - tstart, qd = qend - qstart) %>%
    select(tchrom, tstart, tend, srd, qchrom, qstart, qend, cid, td, qd)

tp1 = ti %>% mutate(type = 'aln')
tp2 = ta4 %>% filter(td + qd <= maxsize) %>%
    select(-td, -qd) %>% mutate(type = 'indel')
tp = tp1 %>% bind_rows(tp2) %>% arrange(tchrom, tstart, tend)

write_tsv(tp, fo, col_names = F)
