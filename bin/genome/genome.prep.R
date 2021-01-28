#!/usr/bin/env Rscript
suppressPackageStartupMessages(require(argparse))

ps <- ArgumentParser(description = 'save genome stats to .rds')
ps$add_argument("genome", default='Zmays_B73', help="genome")
ps$add_argument("fo", help="output (.rds) file")
ps$add_argument("--dirg", default='~/projects/genome/data',
                    help="genome directory [default: %(default)s]")
args <- ps$parse_args()

genome = args$genome
dirg = args$dirg

#require(devtools)
#load_all('~/git/rmaize')
require(tidyverse)

flattern_gcoord_prepare <- function(t_size, gap=2e7) {
    #{{{
    offsets = c(0, cumsum(t_size$size + gap)[-nrow(t_size)])
    t_size %>% mutate(offset = offsets) %>%
        mutate(start=1+offset, end=size+offset, pos=(start+end)/2)
    #}}}
}

dirw = file.path(dirg, genome)
fi = file.path(dirw, '15_intervals', '01.chrom.bed')
bed.chrom = read_tsv(fi, col_names = c("chrom", "start", "end")) %>%
    mutate(start = start + 1, end = as.numeric(end))
size = bed.chrom %>% select(chrom, size=end)
chrom = flattern_gcoord_prepare(size, gap=0)

fi = file.path(dirw, '15_intervals', '11.gap.bed')
bed.gap = read_tsv(fi, col_names = c("chrom", "start", "end")) %>%
    mutate(start = start + 1, end = as.numeric(end))

fi = file.path(dirw, '50_annotation', '15.tsv')
gene.loc = read_tsv(fi, col_types = 'ccccciic')
fi = file.path(dirw, '50_annotation', '15.desc.tsv')
gene.desc = read_tsv(fi)

gene.size.rna = gene.loc %>%
    group_by(gid, tid) %>%
    summarise(ttype = ttype[1], chrom=chrom[1],
              start = min(start), end = max(end), srd = srd[1],
              size.rna = end-start+1) %>% ungroup()
gene.size.exon = gene.loc %>% filter(etype == 'exon') %>%
    group_by(gid, tid) %>%
    summarise(size.exon = sum(end - start + 1)) %>% ungroup()
gene.size.cds = gene.loc %>% filter(etype == 'CDS') %>%
    group_by(gid, tid) %>%
    summarise(size.cds = sum(end - start + 1)) %>% ungroup()
gene = gene.size.rna %>%
    left_join(gene.size.exon, by=c('gid','tid')) %>%
    left_join(gene.size.cds, by=c('gid','tid')) %>%
    left_join(gene.desc, by=c('gid'='id')) %>%
    replace_na(list(size.cds=0,note1='',note2=''))
cat(sprintf("%5d genes for %s\n", nrow(gene), genome))

gidx = gene %>% filter(ttype=='mRNA') %>%
    arrange(chrom,start,end) %>%
    select(gid,tid,chrom,start,end) %>%
    mutate(idx=1:n()) %>%
    group_by(chrom) %>% mutate(cidx=1:n()) %>% ungroup()
gidx.chrom = gidx %>% group_by(chrom) %>%
    summarise(start=min(idx), end=max(idx), pos=(start+end)/2) %>% ungroup()

res = list(chrom=chrom, gap=bed.gap, gene=gene, gene.loc=gene.loc,
           gidx=gidx, gidx.chrom=gidx.chrom)
fo = file.path(dirw, '55.rds')
saveRDS(res, file = args$fo)
