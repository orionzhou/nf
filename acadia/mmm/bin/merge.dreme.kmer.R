#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser(description = 'merge parsed dreme kmer outputs')
parser$add_argument("fi", nargs='+', help = "dreme kmer output file(s)")
parser$add_argument("-o", dest = 'fo', metavar = 'output',
                    nargs=1, default="dreme_kmer.tsv",
                    help = "output file [default: %(default)s]")
args <- parser$parse_args()

fis = args$fi
fo = args$fo

require(tidyverse)
require(tidygenomics)

nfile = length(fis)
ti = tibble(fi = fis) %>%
    mutate(fname = map_chr(fi, basename)) %>%
    mutate(lid = str_replace(fname, '[\\.]gz$', '')) %>%
    mutate(lid = str_replace(lid, '[\\.]\\w+$', '')) %>%
    select(-fname)

to = ti %>% mutate(r = map(fi, read_tsv, col_types='cccciidd')) %>%
    select(lid, r) %>% unnest(r)

write_tsv(to, fo, na='')

