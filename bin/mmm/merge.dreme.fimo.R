#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

ps <- ArgumentParser(description = 'merge fimo.py outputs')
ps$add_argument("fi", nargs='+', help = "fimo.py output file(s)")
ps$add_argument("-o", dest = 'fo', metavar = 'output',
                nargs=1, default="out.rds",
                help = "output file [default: %(default)s]")
args <- ps$parse_args()

fis = args$fi
fo = args$fo

require(tidyverse)

nfile = length(fis)
ti = tibble(fi = fis) %>%
    mutate(fname = map_chr(fi, basename)) %>%
    mutate(lid = str_replace(fname, '[\\.]gz$', '')) %>%
    mutate(lid = str_replace(lid, '[\\.]\\w+$', '')) %>%
    select(-fname)

read_fimo <- function(fi) {
    #{{{
    ti = read_tsv(fi)
    if(nrow(ti) == 0) {
        tibble()
    } else {
        if (is.na(ti$motif_alt_id[1])) ti = ti %>% mutate(ti, motif_alt_id = motif_id)
        ti = ti %>%
            select(mid=2,sid=3,start=4,end=5,srd=6,score=7,pval=8) %>%
            filter(score > 0) %>%
            separate(sid, c('sid','sstart','send'), sep='-') %>%
            select(-sstart, -send) %>%
            mutate(chrom = str_c(mid, sid, sep='_'))
        #
        ti2 = genome_cluster(ti, by=c('chrom','start','end'), cluster_column_name='cid') %>%
            group_by(cid, chrom, mid, sid) %>%
            summarise(start=min(start), end=max(end), score=max(score), pval=min(pval)) %>% ungroup() %>%
            mutate(size=end-start+1, start2 = 1, end2 = 2)
        #
        ti3 = genome_cluster(ti2, by=c("chrom",'start2','end2'), cluster_column_name='cid') %>%
            select(cid, mid, sid, size, score, pval) %>%
            group_by(cid, mid, sid) %>%
            summarise(nseq=n(), size=sum(size), score=sum(score), pval=min(pval)) %>% ungroup() %>%
            select(-cid)
            ti3
    }
    #}}}
}
read_fimo_py <- function(fi) {
    #{{{
    ti = read_tsv(fi, col_names=c('id','start','end'))
    if(nrow(ti) == 0) {
        tibble()
    } else {
        ti %>%
            separate('id', c('mid','gid'), sep='%') %>%
            separate('mid', c('mid','mtf'), sep='-') %>%
            mutate(pos = round((start+1+end)/2)) %>%
            select(mid, gid, pos)
    }
    #}}}
}

to = ti %>%
    mutate(x = map(fi, read_fimo_py)) %>% select(-fi) %>%
    unnest(x) %>%
    mutate(mid = str_c(lid,mid,sep='_')) %>%
    group_by(lid, mid) %>% nest() %>% ungroup() %>%
    rename(loc=data)

saveRDS(to, fo)

