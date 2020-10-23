#!/usr/bin/env Rscript
suppressPackageStartupMessages(require(argparse))

ps <- ArgumentParser(description = 'create bioconductor txDB object')
ps$add_argument("gff", help="input gff file")
ps$add_argument("sizes", help="input chrom size file")
ps$add_argument("fo", help="output (.sqlite) file")
ps$add_argument("--organism", default='unknown',
                help="organism [default: %(default)s]")
args <- ps$parse_args()

f_gff = args$gff
f_size = args$sizes
fo = args$fo

require(tidyverse)
require(GenomicFeatures)

chromInfo = read_tsv(f_size, col_names=c('chrom','length'))
txdb = makeTxDbFromGFF(f_gff, format='gff3',
                       organism=args$organism, chrominfo=chromInfo)

saveDb(txdb, file=fo)

