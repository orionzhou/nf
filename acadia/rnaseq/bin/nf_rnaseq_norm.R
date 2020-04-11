#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser(description = 'Normalize raw read count matrix to CPM and FPKM by DESeq2 and edgeR')
parser$add_argument("fi", help="input file (*.rds)")
parser$add_argument("fo", help="output file (*.rds)")
parser$add_argument("--meta", default='none',
                    help="new meta table [default: %(default)s]")
parser$add_argument("--rcfg", default='none',
                    help="genome configuration file (*.rds) [default: %(default)s]")
args <- parser$parse_args()

fi = args$fi
fo = args$fo
meta = args$meta

require(DESeq2)
require(edgeR)
require(devtools)
load_all("~/git/rmaize")

res = readRDS(fi)
th=res$th; fcnt=res$fcnt

if(meta != 'none' & file.access(meta) != -1 )
    th = read_tsv(meta)

size.gene = if( args$rcfg == '') F : readRDS(args$rcfg)$gene %>% select(gid, size=size.exon)

t_rc = fcnt %>% filter(SampleID %in% th$SampleID)

rn = readcount_norm(t_rc, size.gene)
tl = rn$tl; tm = rn$tm

ths = th %>% distinct(Tissue, Genotype, Treatment, Replicate) %>%
    dplyr::count(Tissue, Genotype, Treatment) %>% dplyr::rename(n_rep=n) %>%
    mutate(nSampleID = sprintf("s%03d", 1:length(Tissue)))
t_map = th %>% inner_join(ths, by = c("Tissue", "Genotype", "Treatment")) %>%
    select(SampleID, nSampleID)
th_m = ths %>% select(SampleID=nSampleID, Tissue, Genotype, Treatment, n_rep)

t_rc_m = t_rc %>% inner_join(t_map, by = 'SampleID') %>%
    mutate(SampleID = nSampleID) %>%
    group_by(gid, SampleID) %>%
    summarise(ReadCount = sum(ReadCount)) %>%
    ungroup()
rn2 = readcount_norm(t_rc_m, size.gene)
tm_m = rn2$tm

x = res
x$th = th; x$tl = tl; x$tm = tm; x$th_m = th_m; x$tm_m = tm_m
saveRDS(x, file = fo)

