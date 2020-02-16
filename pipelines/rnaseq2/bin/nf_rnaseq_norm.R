#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser(description = 'Normalize raw read count matrix to CPM and FPKM by DESeq2 and edgeR')
parser$add_argument("fi", nargs=1, help="input file (*.rds)")
parser$add_argument("fo", nargs=1, help="output file (*.rds)")
parser$add_argument("--meta", default='none',
                    help="new meta table [default: %(default)s]")
parser$add_argument("--config", default='none',
                    help="genome configuration file (*.rds) [default: %(default)s]")
args <- parser$parse_args()

fi = args$fi
fo = args$fo
meta = args$meta
f_cfg = args$config

require(DESeq2)
require(edgeR)
require(devtools)
load_all("~/git/rmaize")

res = readRDS(fi)
th=res$th; fcnts=res$fcnts; salmon=res$salmon; salmonT=res$salmon_transcript

if(meta != 'none' & file.access(meta) != -1 )
    th = read_tsv(meta)

size.gene = F
if( f_cfg != 'none')
    size.gene = readRDS(f_cfg)$gene %>% select(gid, size=size.exon)

t_rc = fcnts %>% filter(SampleID %in% th$SampleID)

res = readcount_norm(t_rc, size.gene)
tl = res$tl; tm = res$tm

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
resm = readcount_norm(t_rc_m, size.gene)
tm_m = resm$tm

r = list(th=th, fcnts=fcnts, salmon=salmon, salmon_transcript=salmonT,
    tl=tl, tm=tm, th_m=th_m, tm_m=tm_m)
saveRDS(r, file = fo)

