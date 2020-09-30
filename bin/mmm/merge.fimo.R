#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser(description = 'merge fimo outputs')
parser$add_argument("fi", nargs='+', help = "fimo output file(s)")
parser$add_argument("-o", dest = 'fo', metavar = 'output',
                    nargs=1, default="fimo.rds",
                    help = "output file [default: %(default)s]")
args <- parser$parse_args()

fis = args$fi
fo = args$fo

require(tidyverse)
require(universalmotif)
require(tidygenomics)

nfile = length(fis)
ti = tibble(fi = fis) %>%
    mutate(fname = map_chr(fi, basename)) %>%
    mutate(lid = str_replace(fname, '[\\.]gz$', '')) %>%
    mutate(lid = str_replace(lid, '[\\.]\\w+$', '')) %>%
    select(-fname)

read_meme0 <- function(fi) {
    #{{{
    r = read_meme(fi)
    if (length(r) == 1)
        list(r)
    else
        r
    #}}}
}
rename_motif <- function(mtf, mid) {
    #{{{
    mtf['altname'] = mtf['name']
    mtf['name'] = mid
    mtf
    #}}}
}

fi = '/home/springer/zhoux379/projects/stress/nf/raw/22b_dreme_fimo/l0421.tsv'
fi = '/home/springer/zhoux379/projects/stress/nf/raw/11_fimo/M08424_2.00.tsv'
read_fimo <- function(fi) {
    #{{{
    ti = read_tsv(fi, comment = "#")
    if(nrow(ti) == 0) {
        list(raw=tibble(), cnt=tibble())
    } else {
        if (is.na(ti$motif_alt_id[1])) ti = ti %>% mutate(ti, motif_alt_id = motif_id)
        ti = ti %>%
            select(mid=2,sid=3,start=4,end=5,srd=6,score=7,pval=8) %>%
            filter(score > 0) %>%
            separate(sid, c('sid','sstart','send'), sep='-') %>%
            select(-sstart, -send)
        ti1 = ti %>%
            mutate(chrom = str_c(mid, sid, sep='_'))
        #
        ti2 = genome_cluster(ti1, by=c('chrom','start','end'), cluster_column_name='cid') %>%
            group_by(cid, chrom, mid, sid) %>%
            summarise(start=min(start), end=max(end), score=max(score), pval=min(pval)) %>% ungroup() %>%
            mutate(size=end-start+1, start2 = 1, end2 = 2)
        #
        ti3 = genome_cluster(ti2, by=c("chrom",'start2','end2'), cluster_column_name='cid') %>%
            select(cid, mid, sid, size, score, pval) %>%
            group_by(cid, mid, sid) %>%
            summarise(nseq=n(), size=sum(size), score=sum(score), pval=min(pval)) %>% ungroup() %>%
            select(-cid)
        list(raw=ti, cnt=ti3)
    }
    #}}}
}

to = ti %>% mutate(r = map(fi, read_fimo)) %>%
    mutate(raw = map(r, 'raw'), cnt = map(r, 'cnt')) %>%
    select(lid, fi, raw, cnt)

saveRDS(to, fo)

