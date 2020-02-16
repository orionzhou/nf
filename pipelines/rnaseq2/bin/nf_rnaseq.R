#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser(description = 'make raw read count matrix using nextflow output')
parser$add_argument("yid", default='none', help="project/study ID [default: %(default)s]")
parser$add_argument("out", nargs=1, help="output file (*.rds)")
parser$add_argument("--metadir", default='/home/springer/zhoux379/projects/barn/data/15_read_list',
                    help="meta dir [default: %(default)s]")
parser$add_argument("--rawdir", default='/home/springer/zhoux379/projects/nf/raw',
                    help="nextflow raw dir [default: %(default)s]")
parser$add_argument("--config", default='none',
                    help="genome configuration file (*.rds) [default: %(default)s]")
args <- parser$parse_args()

yid = args$yid
f_out = args$out
dir_meta = args$metadir
dir_raw = args$rawdir
f_cfg = args$config

require(devtools)
load_all("~/git/rmaize")

f_sam = sprintf("%s/%s.tsv", dir_meta, yid)
if( file.access(f_sam) == -1 )
    stop(sprintf("file ( %s ) cannot be accessed", f_sam))
th = read_tsv(f_sam)

size.gene = F
if( f_cfg != 'none')
    size.gene = readRDS(f_cfg)$gene %>% select(gid, size=size.exon)

# featureCounts
fi1 = sprintf("%s/%s/featureCounts/merged_gene_counts.txt", dir_raw, yid)
if( file.access(fi1) == -1 )
    stop(sprintf("file ( %s ) cannot be accessed", fi1))
fcnts = read_tsv(fi1) %>% select(-2) %>%
    rename(gid = 1) %>% gather(SampleID, ReadCount, -gid) %>%
    filter(SampleID %in% th$SampleID)

# salmon
fi1 = sprintf("%s/%s/salmon/salmon_merged_gene_counts.csv", dir_raw, yid)
if( file.access(fi1) == -1 )
    stop(sprintf("file ( %s ) cannot be accessed", fi1))
ti1 = read_csv(fi1) %>% rename(gid = 1) %>%
    gather(SampleID, ReadCount, -gid)
fi1 = sprintf("%s/%s/salmon/salmon_merged_gene_tpm.csv", dir_raw, yid)
if( file.access(fi1) == -1 )
    stop(sprintf("file ( %s ) cannot be accessed", fi1))
ti2 = read_csv(fi1) %>% rename(gid = 1) %>%
    gather(SampleID, TPM, -gid)

fi1 = sprintf("%s/%s/salmon/salmon_merged_transcript_counts.csv", dir_raw, yid)
if( file.access(fi1) == -1 )
    stop(sprintf("file ( %s ) cannot be accessed", fi1))
ti3 = read_csv(fi1) %>% rename(tid = 1) %>%
    gather(SampleID, ReadCount, -tid)
fi1 = sprintf("%s/%s/salmon/salmon_merged_transcript_tpm.csv", dir_raw, yid)
if( file.access(fi1) == -1 )
    stop(sprintf("file ( %s ) cannot be accessed", fi1))
ti4 = read_csv(fi1) %>% rename(tid = 1) %>%
    gather(SampleID, TPM, -tid)

salmon = ti1 %>% inner_join(ti2, by=c("SampleID",'gid'))
salmon_transcript = ti3 %>% inner_join(ti4, by=c("SampleID",'tid'))

res = list(th=th, fcnts=fcnts, salmon=salmon, salmon_transcript=salmon_transcript)
saveRDS(res, file = f_out)

