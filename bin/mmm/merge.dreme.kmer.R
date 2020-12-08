#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

ps <- ArgumentParser(description = 'merge parsed dreme kmer outputs')
ps$add_argument("fi", nargs='+', help = "dreme kmer output file(s)")
ps$add_argument("--fok", default="kmer.tsv",
                help = "output kmer file [default: %(default)s]")
ps$add_argument("--fom", default="kmer.motif.tsv",
                help = "output kmer motif file [default: %(default)s]")
ps$add_argument("--fos", default="kmer.txt",
                help = "output kmer sequence file [default: %(default)s]")
args <- ps$parse_args()

fis = args$fi
fok = args$fok
fom = args$fom
fos = args$fos

require(tidyverse)
#require(tidygenomics)

nfile = length(fis)
ti = tibble(fi = fis) %>%
    mutate(fname = map_chr(fi, basename)) %>%
    mutate(lid = str_replace(fname, '[\\.]gz$', '')) %>%
    mutate(lid = str_replace(lid, '[\\.]\\w+$', '')) %>%
    select(-fname)

rename_mid <- function(lid, mid) sprintf("%s_%03d", lid, as.numeric(str_replace_all(mid,'.*REME-','')))
tk = ti %>% mutate(r = map(fi, read_tsv, col_types='cccciiiidd')) %>%
    select(lid, r) %>% unnest(r) %>%
    mutate(mid = map2_chr(lid, mid, rename_mid))
tom = tk %>% filter(!is.na(re)) %>% select(-re)
tok = tk %>% filter(is.na(re)) %>% select(-re) %>%
    group_by(lid, mid) %>% mutate(kid=sprintf("%s_%d", mid, 1:n())) %>%
    ungroup() %>%
    select(lid,mid,kid,everything())
tos = tok %>% distinct(seq) %>% arrange(seq) %>% pull(seq)

write_tsv(tom, fom, na='')
write_tsv(tok, fok, na='')
write(tos, fos)

