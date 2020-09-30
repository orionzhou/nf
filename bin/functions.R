require(devtools)
load_all("~/git/rmaize")
require(ape)
require(ggtree)
require(ggforce)
require(Rtsne)
dirp = '~/projects/nf'
dird = file.path(dirp, 'data')
#dirc = file.path(dird, 'cache')
dirr = file.path(dird, 'raw')
f_cfg = '~/projects/master.xlsx'
t_cfg = read_xlsx(f_cfg, sheet='barn', col_names=T) %>%
    #filter(libtype %in% c("chipseq",'dapseq','atacseq')) %>%
    select(yid,libtype,author,study,genotype,tissue,n,ref) %>%
    replace_na(list(ref='Zmays_B73')) %>%
    mutate(lgd = sprintf("%s %s", str_to_title(author),study))

get_fastq <- function(sid, paired, read1, yid) {
    #{{{
    diri = file.path("/scratch.global/zhoux379/barn/data/fastq", yid)
    if(paired) {
        suf = ifelse(read1, "_1.fq.gz", "_2.fq.gz")
        sprintf("%s/%s%s", diri, sid, suf)
    } else {
        ifelse(read1, sprintf("%s/%s.fq.gz", diri, sid), "")
    }
    #}}}
}

read_samplelist <- function(yid, diri = '~/projects/barn/data/15_read_list')
    read_tsv(sprintf("%s/%s.tsv", diri, yid))
read_design <- function(yid, diri = '~/projects/nf/data/design')
    read_csv(sprintf("%s/%s.csv", diri, yid))



