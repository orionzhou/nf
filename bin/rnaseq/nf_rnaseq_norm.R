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

require(tidyverse)
fi = args$fi
fo = args$fo
meta = args$meta

require(DESeq2)
require(edgeR)
readcount_norm <- function(t_rc, t_gs = F) {
    #{{{ normalize
    smMap = t_rc %>% distinct(SampleID) %>% dplyr::rename(oSampleID=SampleID) %>%
        mutate(nSampleID = str_replace_all(oSampleID, '[^a-zA-z0-9_]', '.'))
    tm = t_rc
    if (! identical(smMap$oSampleID, smMap$nSampleID))
        tm = tm %>% inner_join(smMap, by=c('SampleID'='oSampleID')) %>%
            select(-SampleID) %>% select(SampleID=nSampleID, everything())
    tw = tm %>%
        select(SampleID, gid, ReadCount) %>%
        spread(SampleID, ReadCount) %>%
        replace(., is.na(.), 0)
    tm = tw %>% gather(SampleID, ReadCount, -gid)
    gids = tw$gid
    twd = data.frame(tw[,-1])
    rownames(twd) = tw$gid
    # nRC with DESeq2
    require(DESeq2)
    th = tm %>% distinct(SampleID) %>% arrange(SampleID)
    th2 = th %>% mutate(sid = SampleID, SampleID = factor(SampleID))
    thd = column_to_rownames(as.data.frame(th2), var = 'sid')
    stopifnot(identical(rownames(thd), colnames(twd)))
    dds = DESeqDataSetFromMatrix(countData=twd, colData=thd, design = ~ 1)
    dds = estimateSizeFactors(dds)
    sf = sizeFactors(dds)
    t_sf = tibble(SampleID = names(sf), sizeFactor = as.numeric(sf))
    t_nrc = counts(dds, normalized = T) %>% as_tibble() %>%
        mutate(gid = names(dds)) %>% gather(SampleID, nRC, -gid)
    # rCPM and CPM with edgeR
    require(edgeR)
    y = DGEList(counts = twd)
    y = calcNormFactors(y, method = 'TMM')
    t_nf = y$samples %>% as_tibble() %>%
        mutate(SampleID = rownames(y$samples)) %>%
        select(SampleID, libSize=lib.size, normFactor=norm.factors)
    t_cpm1 = cpm(y) %>% as_tibble() %>% mutate(gid = rownames(cpm(y))) %>%
        select(gid, everything()) %>%
        gather(SampleID, CPM, -gid)
    t_cpm2 = cpm(y, normalized = F) %>% as_tibble() %>%
        mutate(gid = rownames(cpm(y))) %>%
        gather(SampleID, rCPM, -gid)
    t_cpm = t_cpm1 %>% inner_join(t_cpm2, by=c('SampleID','gid'))
    # rFPKM & FPKM
    if(is.list(t_gs)) {
        t_tpm = tm %>% left_join(t_gs, by='gid') %>%
            mutate(RPK = ReadCount / (size/1000)) %>%
            group_by(SampleID) %>% mutate(rTPM = RPK / sum(RPK) * 1e6) %>%
            ungroup() %>% inner_join(t_nf, by = 'SampleID') %>%
            mutate(TPM = rTPM / normFactor) %>%
            select(SampleID, gid, TPM, rTPM)
        t_cpm = t_cpm %>% left_join(t_gs, by = 'gid') %>%
            mutate(FPKM = CPM / (size / 1000)) %>%
            mutate(rFPKM = rCPM / (size / 1000)) %>%
            select(-size) %>%
            inner_join(t_tpm, by=c('SampleID','gid'))
    }
    #
    tl = th %>% inner_join(t_sf, by = 'SampleID') %>%
        inner_join(t_nf, by = 'SampleID')
    tm = tm %>%
        left_join(t_nrc, by = c('SampleID','gid')) %>%
        left_join(t_cpm, by = c('SampleID','gid'))
    stopifnot(nrow(tm) == length(gids) * nrow(th))
    if (! identical(smMap$oSampleID, smMap$nSampleID)) {
        tl = tl %>% inner_join(smMap, by=c('SampleID'='nSampleID')) %>%
            dplyr::select(-SampleID) %>% dplyr::select(SampleID=oSampleID, everything())
        tm = tm %>% inner_join(smMap, by=c('SampleID'='nSampleID')) %>%
            dplyr::select(-SampleID) %>% dplyr::select(SampleID=oSampleID, everything())
    }
    list(tl = tl, tm = tm, dds=dds)
    #}}}
}

res = readRDS(fi)
th=res$th; fcnt=res$fcnt

if(meta != 'none' & file.access(meta) != -1 )
    th = read_tsv(meta)

size.gene = F
if(args$rcfg != 'none')
    size.gene = readRDS(args$rcfg)$gene %>% select(gid, size=size.exon)

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

