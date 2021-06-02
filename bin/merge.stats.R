#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

ps <- ArgumentParser(description = 'merge stats tables')
ps$add_argument("fi", nargs='+', help = "input stats file(s)")
ps$add_argument("-o", dest = 'fo', metavar = 'output',
                 nargs=1, default="stats.tsv",
                 help = "output file [default: %(default)s]")
ps$add_argument("--opt", metavar = 'option',
                 nargs=1, default="bam_stat",
                 help = "stats file format [default: %(default)s]")
args <- ps$parse_args()

fis = args$fi
fo = args$fo
opt = args$opt

require(tidyverse)

nfile = length(fis)
ti = tibble(fi = fis) %>%
    mutate(fname = map_chr(fi, basename)) %>%
    mutate(sid = str_replace(fname, '[\\.]gz$', '')) %>%
    mutate(sid = str_replace(sid, '[\\.]\\w+$', ''))

read_featurecount <- function(fi)
    read_tsv(fi, skip = 1) %>%
        select(gid = 1, cnt = 7)

read_mmquant <- function(fi)
    read_tsv(fi, col_types = 'ci') %>%
        select(gid = 1, cnt = 2)

read_salmon <- function(fi)
    read_csv(fi, skip = 1) %>%
        select(gid = 1, val = 2)

read_snpbinner <- function(fi) {
    #{{{
    ti = read_csv(fi, col_names=F)
    ti2 = tibble(start=as.integer(ti[1,-1]), end=as.integer(ti[2,-1]))
    ti2b = t(column_to_rownames(ti[-c(1:3),], "X1")) %>% as_tibble()
    ti3 = ti2 %>% bind_cols(ti2b) %>%
        gather(SampleID, gt, -start, -end)
    ti3
    #}}}
}

bin2tib <- function(binstr) {
    #{{{
    parts = str_split(str_remove(binstr, ',$'), ',')[[1]]
    sid = parts[1]
    ps = parts[2:length(parts)]
    np = length(parts)
    gts = ps[seq(2, np-1, by=2)]
    starts = as.integer(ps[seq(1, np-2, by=2)])
    ends = as.integer(ps[seq(3, np, by=2)])
    tibble(sid=sid, start=starts, end=ends, gt=gts)
    #}}}
}

if (opt == 'bam_stat') {
    #{{{
    to = ti %>%
        mutate(data = map(fi, read_tsv, col_names=c('type','cnt'), col_types='cd')) %>%
        select(sid, data) %>%
        unnest(data) %>%
        spread(type, cnt)
    write_tsv(to, fo)
    #}}}
} else if (opt == 'featurecounts') {
    #{{{
    to = ti %>%
        mutate(data = map(fi, read_featurecount)) %>%
        select(sid, data) %>%
        mutate(sid = str_replace(sid, "_gene\\.featureCounts$", '')) %>%
        unnest(data)
    saveRDS(to, file=fo)
    #}}}
} else if (opt == 'salmon') {
    #{{{
    to = ti %>%
        mutate(sid = str_replace(sid, "_salmon\\w+$", '')) %>%
        mutate(data = map(fi, read_salmon)) %>%
        select(sid, data) %>%
        unnest(data)
    saveRDS(to, file=fo)
    #}}}
} else if (opt == 'mmquant') {
    #{{{
    to = ti %>%
        mutate(data = map(fi, read_mmquant)) %>%
        select(sid, data) %>%
        unnest(data)
    saveRDS(to, file=fo)
    #}}}
} else if (opt == 'snpbinner') {
    #{{{
    to = ti %>% rename(rid=sid) %>%
        mutate(data = map(fi, read_snpbinner)) %>%
        select(rid, data) %>%
        unnest(data)
    to2 = ti %>% rename(rid=sid) %>%
        mutate(fi = str_replace(fi, ".csv", ".txt")) %>%
        mutate(data = map(fi, read_tsv, col_names = F)) %>%
        select(rid, data) %>% unnest(data) %>%
        mutate(res = map(X1, bin2tib)) %>%
        select(rid, res) %>% unnest(res)
    saveRDS(list(cp=to2, bin=to), file=fo)
    #}}}
} else if (opt == 'ase_gene') {
    #{{{
    to = ti %>%
        mutate(data = map(fi, read_featurecount)) %>%
        select(sid, data) %>%
        unnest(data) %>%
        separate(sid, c('sid','allele'), sep='[\\.]', extra='merge') %>%
        mutate(allele=str_c("allele", allele, sep='')) %>%
        spread(allele, cnt)
    saveRDS(to, file=fo)
    #}}}
} else if (opt == 'ase_snp') {
    #{{{
    to = ti %>%
        mutate(data = map(fi, read_tsv)) %>%
        select(sid, data) %>%
        unnest(data) %>%
        mutate(allele1 = ifelse(gt=='1|0', altsupport, refsupport)) %>%
        mutate(allele2 = ifelse(gt=='1|0', refsupport, altsupport)) %>%
        select(sid,chr,pos,ref,alt,gt,allele1,allele2)
    saveRDS(to, file=fo)
    #}}}
} else if (opt == 'rds') {
    to = ti %>%
        mutate(data = map(fi, readRDS)) %>%
        select(sid, data) %>%
        unnest(data)
    saveRDS(to, file=fo)
} else {
    stop("unknown option", opt, "\n")
}
