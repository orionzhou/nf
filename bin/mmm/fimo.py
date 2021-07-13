#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import random

from jcvi.apps.base import sh, popen, mkdir
from jcvi.formats.base import read_block, is_number, get_number

off = 4050
range_dict = {
    'TSS:-500': [[1500,2000]],
    'TSS:+500': [[2000,2500]],
    'TSS:-/+500':[[1500,2500]],
    'TSS:-1k':[[1000,2000]],
    'TSS:+1k':[[2000,3000]],
    'TSS:-/+1k':[[1000,3000]],
    'TSS:-2k':[[0,2000]],
    'TSS:+2k':[[2000,4000]],
    'TSS:-/+2k':[[0,4000]],
    'TTS:-500':[[off+1500,off+2000]],
    'TTS:+500':[[off+2000,off+2500]],
    'TTS:-/+500':[[off+1500,off+2500]],
    'TTS:-1k':[[off+1000,off+2000]],
    'TTS:+1k':[[off+2000,off+3000]],
    'TTS:-/+1k':[[off+1000,off+3000]],
    'TTS:-2k':[[off+0,off+2000]],
    'TTS:+2k':[[off+2000,off+4000]],
    'TTS:-/+2k':[[off+0,off+4000]],
    'TSS:-/+2k,TTS:-500':[[0,4000],[off+1500,off+2000]],
    'TSS:-/+2k,TTS:-1k':[[0,4000],[off+1000,off+2000]],
    'TSS:-/+2k,TTS:-2k':[[0,4000],[off+0,off+2000]],
    'TSS:-/+2k,TTS:-/+500':[[0,4000],[off+1500,off+2500]],
    'TSS:-/+2k,TTS:-/+1k':[[0,4000],[off+1000,off+3000]],
    'TSS:-/+2k,TTS:-/+2k':[[0,4000],[off+0,off+4000]]
}

def read_meme(fi):
    mtfs = []
    fhi = open(fi, 'r')
    for head,content in read_block(fhi, 'MOTIF'):
        ps = head.split(' ')
        pre, mid = ps[:2]
        score = ''
        if len(ps) >= 3:
            score = ps[2]
        #mtf = mid.split("-")[1]
        if is_number(score):
            score = float(score)
        width = len(content)-2
        mtfs.append([mid,width,score])
        #print(mid,'\t',width)
    return mtfs

def get_motif_width(fi):
    fhi = open(fi, 'r')
    for line in fhi:
        ps = line.strip().split("\t")
        return int(ps[2])-int(ps[1])

def read_motif(fm, mtf_str):
    if fm == 'none': return []
    mtfs0 = read_meme(fm)
    mtfs1 = [x[0] for x in mtfs0]
    print(f'{len(mtfs1)} motifs read')
    mtfs = []
    if mtf_str == 'all':
        mtfs = mtfs0
    elif mtf_str.startswith("top"):
        n = int(mtf_str.replace("top",''))
        mtfs = mtfs0[:n]
    else:
        mids = set(mtf_str.split(","))
        mtfs = [ [x,wd,score] for x,wd,score in mtfs0 if x in mids ]
    if len(mtfs) == 0:
        print(f"zero motifs found in {fm}")
        sys.exit(1)
    return mtfs

def read_gene_list(fg):
    glst = []
    for line in open(fg,'r'):
        line = line.rstrip("\n")
        if not line: continue
        ps = line.split()
        gid, status = ps[0], 0
        if gid == 'gid': continue
        if len(ps) >= 2: status = ps[1]
        glst.append([gid,status])
    return glst

def locate(args):
    fi, fo = args.fi, args.fo
    seq = args.seq
    mtfs = read_motif(fi, args.motif)
    #
    mtf_str = " ".join([f'--motif {mid}' for mid,wd,score in mtfs])
    pre = f"tmp.lc{random.randrange(1000)}"
    #
    sh(f'fimo --bfile --motif-- {mtf_str} --thresh {args.eval} --skip-matched-sequence --text {fi} {seq} > {pre}_0.txt')
    for mid,wd,score in mtfs:
        sh(f'grep -P "^{mid}\t" {pre}_0.txt > {pre}_0a.txt')
        #
        score_thresh = score
        if not score:
            xh = popen(f'cut -f7 {pre}_0a.txt | sed \'1d\' | sort -k1,1nr | head')
            max_score = float(xh.readline().decode("utf-8").strip())
            score_thresh = max_score * args.score_thresh
        #
        sh("bioawk -tH '{if($7>%f) {print $1\"%%\"$3, $4-1, $5}}' %s_0a.txt > %s_1.bed" % (score_thresh, pre, pre))
        hwd = round(wd * args.motif_frac)
        if os.stat(f"{pre}_1.bed").st_size == 0:
            sh(f'touch {pre}_4_{mid}.bed')
        else:
            sh(f'sortBed -i {pre}_1.bed | mergeBed > {pre}_2.bed')
            sh(f'bedtools makewindows -w {wd} -b {pre}_2.bed > {pre}_3.bed')
            sh(f'bed.py filter --minsize {hwd} {pre}_3.bed > {pre}_4_{mid}.bed')
    sh(f'cat {pre}_4_*.bed > {fo}')
    if not args.debug:
        sh(f'rm -rf {pre}_*')

def filter(args):
    fi, fg, fo = args.fi, args.fg, args.fo
    epi, bin = args.epi, args.bin
    glst = read_gene_list(fg)
    #
    pre = f"tmp.fl{random.randrange(1000)}"
    fht = open(f"{pre}_1.bed", 'w')
    locs = range_dict[bin]
    for gid, status in glst:
        for b, e in locs:
            fht.write(f'{gid}\t{b}\t{e}\n')
    fht.close()
    #
    sh(f'intersectBed -u -f 0.5 -a {fi} -b {pre}_1.bed > {pre}_2.bed')
    if epi == 'raw':
        sh(f"cp {pre}_2.bed {fo}")
    else:
        fe = ''
        if epi == 'umr':
            fe = args.umr
        elif epi == 'acrL':
            fe = args.acrL
        elif epi == 'acrE':
            fe = args.acrE
        else:
            print(f"unknown {epi} option")
            sys.exit(1)
        sh(f'intersectBed -u -f 0.5 -a {pre}_2.bed -b {fe} > {fo}')
    sh(f"rm {pre}_*")

def bed2wide(args):
    fg, fb, fo = args.fg, args.fb, args.fo
    mod, fm, nfea = args.mod, args.motif, args.nfea
    glst = read_gene_list(fg)
    mtfs = read_motif(fm, nfea)
    #
    mc = dict()
    for line in open(fb,'r'):
        gid,b,e,mid = line.strip().split("\t")[:4]
        if gid not in mc: mc[gid] = dict()
        if mid not in mc[gid]: mc[gid][mid] = 0
        mc[gid][mid] += 1
    #
    fho = open(fo,'w')
    mid_str = "\t".join([x for x,wd,score in mtfs])
    fho.write(f"gid\tstatus\t{mid_str}\n")
    for gid,status in glst:
        lstr = f"{gid}\t{status}"
        for mid,wd,score in mtfs:
            if gid not in mc or mid not in mc[gid]:
                lstr += "\t0"
            else:
                cnt = 1 if mc[gid][mid]>1 and mod=='zoops' else mc[gid][mid]
                lstr += f"\t{cnt}"
        fho.write(f"{lstr}\n")
    fho.close()

def prepare_ml(args):
    fg, fm, fo = args.fg, args.fm, args.fo
    db = args.db
    bin, epi, nfea, mod = args.bin, args.epi, args.nfea, args.mod
    fmt = args.fmt
    pre = f"tmp.pm{random.randrange(1000)}"
    sh(f"sed '1d' {fg} |cut -f1 > {pre}_0.txt")
    sh(f"fasta.py extract --list {db} {pre}_0.txt {pre}_1.fas")
    debug = '--debug' if args.debug else ''
    sh(f"fimo.py locate {debug} --motif {nfea} {fm} --score_thresh {args.score_thresh} --motif_frac {args.motif_frac} {pre}_1.fas {pre}_2.bed")
    sh("bioawk -t '{split($1,a,/%%/); print a[2], $2, $3, a[1]}' %s_2.bed | sort -k1,1 -k2,2n -k3,3n > %s_3.bed" % (pre, pre))
    sh(f"fimo.py filter {pre}_3.bed --bin '{bin}' --epi {epi} {fg} {pre}_4.bed")
    if fmt == 'long':
        sh(f"cp {pre}_4.bed {fo}")
    else:
        sh(f"fimo.py bed2wide --mod {mod} --motif {fm} --nfea {nfea} {fg} {pre}_4.bed {fo}")
    if not args.debug:
        sh(f"rm {pre}_*")

if __name__ == "__main__":
    import argparse
    ps = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = "fimo utilities"
    )
    sp = ps.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("locate",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            help = "run fimo to find given motifs in input sequences")
    sp1.add_argument('fi', help = 'input motif file (.meme/.dreme/.streme)')
    sp1.add_argument('seq', help = 'sequence file')
    sp1.add_argument('fo', help = 'output file')
    sp1.add_argument('--motif', default='all', help = 'motif ID / option')
    sp1.add_argument('--score_thresh', type=float, default=0.7, help = 'minimum (relative) score threshold')
    sp1.add_argument('--eval', type=float, default=1e-4, help = 'minimum (relative) score threshold')
    sp1.add_argument('--motif_frac', type=float, default=0.8, help = 'fraction of motif to be counted')
    sp1.add_argument('--debug', action='store_true', help = 'do not delete temporary files')
    sp1.set_defaults(func = locate)

    sp1 = sp.add_parser("filter",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            help = "filter BED file using window size / epigenetic marks")
    sp1.add_argument('fi', help = 'input BED file')
    sp1.add_argument('fg', help = 'gene list / status table')
    sp1.add_argument('fo', help = 'output BED file')
    sp1.add_argument('--bin', default='TSS:-/+2k,TTS:-/+2k', help = 'window option')
    sp1.add_argument('--epi', default='raw', choices=['raw','umr','acrE','acrL'], help = 'epigenetic filter')
    sp1.add_argument('--umr', default='%s/data/21_seq/15.UMR.bed' % os.environ['st'], help = 'UMR bed file')
    sp1.add_argument('--acr', default='%s/data/21_seq/15.ACR.bed' % os.environ['st'], help = 'ACR bed file')
    sp1.set_defaults(func = filter)

    sp1 = sp.add_parser("bed2wide",
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        help = "convert BED file to machine learing tables")
    sp1.add_argument('fg', help = 'gene list / status table')
    sp1.add_argument('fb', help = 'motif BED file')
    sp1.add_argument('fo', help = 'output table for ML')
    sp1.add_argument('--mod', default='zoops', choices=['zoops','anr'], help = 'encoding option')
    sp1.add_argument('--motif', default='none', help = 'input motif file (.meme/.dreme/.streme)')
    sp1.add_argument('--nfea', default='all', help = 'motif ID / option')
    sp1.set_defaults(func = bed2wide)

    sp1 = sp.add_parser("prepare_ml",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            help = "pipeline to find motifs and output in BED / ML input table")
    sp1.add_argument('fg', help = 'gene list / status table')
    sp1.add_argument('fm', help = 'input motif (.meme/.dreme/.streme) file')
    sp1.add_argument('fo', help = 'output BED / TSV file')
    sp1.add_argument('--db', default='%s/data/21_seq/02.fas' % os.environ['st'], help = 'database (fasta file) to search')
    sp1.add_argument('--bin', default='TSS:-/+2k', help = 'window option')
    sp1.add_argument('--epi', default='raw', choices=['raw','umr','acrE','acrL'], help = 'epigenetic filter')
    sp1.add_argument('--nfea', default='all', help = 'motif ID / option')
    sp1.add_argument('--mod', default='zoops', choices=['zoops','anr'], help = 'encoding option')
    sp1.add_argument('--umr', default='%s/data/21_seq/15.UMR.bed' % os.environ['st'], help = 'UMR bed file')
    sp1.add_argument('--acr', default='%s/data/21_seq/15.ACR.bed' % os.environ['st'], help = 'ACR bed file')
    sp1.add_argument('--score_thresh', type=float, default=0.7, help = 'minimum (relative) fimo score threshold')
    sp1.add_argument('--motif_frac', type=float, default=0.8, help = 'fraction of motif to be counted')
    sp1.add_argument('--fmt', default='wide', choices=['long','wide'], help = 'output format')
    sp1.add_argument('--debug', action='store_true', help = 'do not delete temporary files')
    sp1.set_defaults(func = prepare_ml)

    args = ps.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        ps.print_help()


