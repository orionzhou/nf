#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser(description = 'convert CTSS BED to bigwig using CAGEr')
parser$add_argument("fi", nargs=1, help = "CTSS BED file")
parser$add_argument("-o1", dest = 'fo1', metavar = 'output bigwig (+)',
                    nargs=1, default="CTSS.plus.bw",
                    help = "bigwig (+) file [default: %(default)s]")
parser$add_argument("-o2", dest = 'fo2', metavar = 'output bigwig (-)',
                    nargs=1, default="CTSS.minus.bw",
                    help = "bigwig (-) file [default: %(default)s]")
parser$add_argument("-genome", dest = 'genome', metavar = 'genome',
                    nargs=1, default="BSgenome.Zmays.B73",
                    help = "BSgenome [default: %(default)s]")
args <- parser$parse_args()

fi = args$fi
fo1 = args$fo1
fo2 = args$fo2
gen = args$genome

require(CAGEr)
#require(BSgenome.Zmays.B73)
x = CAGEr:::import.bedCTSS(fi)

x2 = as.data.frame(x)
colnames(x2)[1]='chr'
x3 = as(x2, 'CAGEexp')
genomeName(x3) = gen
CAGEr:::exportCTSStoBedGraph(x3, values = "raw", format = "BigWig")

system(sprintf("mv score.CTSS.raw.plus.bw %s", fo1))
system(sprintf("mv score.CTSS.raw.minus.bw %s", fo2))
