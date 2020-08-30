#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser(description = 'merge dreme outputs')
parser$add_argument("fi", nargs='+', help = "dreme output file(s)")
parser$add_argument("-o", dest = 'fo', metavar = 'output',
                    nargs=1, default="out.rds",
                    help = "output file [default: %(default)s]")
parser$add_argument("--meme", dest='fm', metavar = 'meme',
                    nargs=1, default="out.meme",
                    help = "merged motifs in meme format [default: %(default)s]")
parser$add_argument("--txt", dest='ft', metavar = 'list',
                    nargs=1, default="out.txt",
                    help = "motif ID list [default: %(default)s]")
args <- parser$parse_args()

fis = args$fi
fo = args$fo
fm = args$fm
ft = args$ft

require(tidyverse)
require(universalmotif)
require(tidygenomics)

#fis = c(
#'/home/springer/zhoux379/projects/stress/nf/raw/22a_dreme_mtf/l0421.dreme',
#'/home/springer/zhoux379/projects/stress/nf/raw/22a_dreme_mtf/l0422.dreme',
#'/home/springer/zhoux379/projects/stress/nf/raw/22b_dreme_fimo/l0421.tsv',
#'/home/springer/zhoux379/projects/stress/nf/raw/22b_dreme_fimo/l0422.tsv'
#)

nfile = length(fis)
ti = tibble(fi = fis) %>%
    mutate(opt = ifelse(str_detect(fi, '\\.dreme$'), 'dreme','fimo')) %>%
    mutate(fname = map_chr(fi, basename)) %>%
    mutate(lid = str_replace(fname, '[\\.]gz$', '')) %>%
    mutate(lid = str_replace(lid, '[\\.]\\w+$', '')) %>%
    select(-fname) %>%
    spread(opt, fi)

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

to1 = ti %>%
    mutate(motif = map(dreme, read_meme0)) %>%
    select(lid, motif) %>% unnest(motif) %>%
    mutate(mid = map_chr(motif, 'altname')) %>%
    mutate(mid = sprintf("%s_%03d", lid, as.integer(str_replace(mid,'DREME-','')))) %>%
    mutate(mtf = map2(motif, mid, rename_motif)) %>%
    select(lid, mid, mtf)

#to2 = ti %>% mutate(hits = map(fimo, read_fimo)) %>%
    #select(lid, hits) %>% unnest(hits) %>%
    #group_by(lid, mid) %>% nest() %>% ungroup() %>%
    #rename(hits = data) %>%
    #mutate(mid = sprintf("%s_%03d", lid, as.integer(str_replace(mid,'DREME-',''))))

to = to1# %>% inner_join(to2, by=c('lid','mid'))
saveRDS(to, fo)

mtfs = to$mtf
write_meme(mtfs, fm, overwrite=T)

tol = to %>% select(mid, lid)
write_tsv(tol, ft)
