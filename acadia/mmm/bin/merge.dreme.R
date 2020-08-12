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

nfile = length(fis)
ti = tibble(fi = fis) %>%
    mutate(fname = map_chr(fi, basename)) %>%
    mutate(sid = str_replace(fname, '[\\.]gz$', '')) %>%
    mutate(sid = str_replace(sid, '[\\.]\\w+$', ''))

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

to = ti %>%
    mutate(motif = map(fi, read_meme0)) %>%
    select(lid=sid, motif) %>% unnest(motif) %>%
    group_by(lid) %>% mutate(i = 1:n()) %>% ungroup() %>%
    mutate(mid = sprintf("%s_%03d", lid, i)) %>%
    mutate(mtf = map2(motif, mid, rename_motif)) %>%
    select(lid, i, mid, mtf)
saveRDS(to, fo)

mtfs = to$mtf
write_meme(mtfs, fm, overwrite=T)

write(to$mid, ft)
