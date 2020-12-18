#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import random

from jcvi.apps.base import sh, mkdir
from jcvi.formats.base import read_block

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
        pre, mid, mid2 = head.split(' ')
        #mtf = mid.split("-")[1]
        width = len(content)-2
        mtfs.append([mid,width])
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
        mtfs = [ [x,wd] for x,wd in mtfs0 if x in mids ]
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
    pre = f"tmp.lc{random.randrange(1000)}"
    for mid,wd in mtfs:
        sh(f'fimo --motif {mid} --thresh {args.thresh} --oc fimo {fi} {seq}')
        cmd = (f'grep \'^{mid}\' fimo/fimo.tsv | '
            'bioawk -tH \'{if($9<0.05) {print $1"%"$3, $4-1, $5}}\' > zzzz.bed')
        sh(cmd)
        sh(f"mv zzzz.bed {pre}_1.bed")
        hwd = round(wd * args.motif_frac)
        if os.stat(f"{pre}_1.bed").st_size == 0:
            sh(f'touch {pre}_4_{mid}.bed')
        else:
            sh(f'sortBed -i {pre}_1.bed | mergeBed > {pre}_2.bed')
            sh(f'bedtools makewindows -w {wd} -b {pre}_2.bed > {pre}_3.bed')
            sh(f'bed.py filter --minsize {hwd} {pre}_3.bed > {pre}_4_{mid}.bed')
    sh(f'cat {pre}_4_*.bed > {fo}')
    sh(f'rm -rf fimo {pre}_*')

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
    mid_str = "\t".join([x for x,wd in mtfs])
    fho.write(f"gid\tstatus\t{mid_str}\n")
    for gid,status in glst:
        lstr = f"{gid}\t{status}"
        for mid,wd in mtfs:
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
    sh(f"fasta.py extract --list {db} {pre}_0.txt > {pre}_1.fas")
    sh(f"fimo.py locate --motif {nfea} {fm} --thresh {args.thresh} --motif_frac {args.motif_frac} {pre}_1.fas {pre}_2.bed")
    sh(f'mv {pre}_2.bed z0.bed')
    sh("bioawk -t '{split($1,a,/%/); print a[2], $2, $3, a[1]}' z0.bed | sort -k1,1 -k2,2n -k3,3n > z.bed")
    sh(f'mv z.bed {pre}_3.bed')
    sh(f"fimo.py filter {pre}_3.bed --bin '{bin}' --epi {epi} {fg} {pre}_4.bed")
    if fmt == 'long':
        sh(f"cp {pre}_4.bed {fo}")
    else:
        sh(f"fimo.py bed2wide --mod {mod} --motif {fm} --nfea {nfea} {fg} {pre}_4.bed {fo}")
    sh(f"rm {pre}_* z0.bed")

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
    sp1.add_argument('--thresh', type=float, default=0.05, help = 'q value threshold')
    sp1.add_argument('--motif_frac', type=float, default=0.8, help = 'fraction of motif to be counted')
    sp1.set_defaults(func = locate)

    sp1 = sp.add_parser("filter",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            help = "filter BED file using window size / epigenetic marks")
    sp1.add_argument('fi', help = 'input BED file')
    sp1.add_argument('fg', help = 'gene list / status table')
    sp1.add_argument('fo', help = 'output BED file')
    sp1.add_argument('--bin', default='TSS:-/+2k,TTS:-/+2k', help = 'window option')
    sp1.add_argument('--epi', default='raw', choices=['raw','umr','acrE','acrL'], help = 'epigenetic filter')
    sp1.add_argument('--umr', default='/home/springer/zhoux379/projects/stress/data/21_seq/15.umr.bed', help = 'UMR bed file')
    sp1.add_argument('--acrE', default='/home/springer/zhoux379/projects/stress/data/21_seq/15.acrE.bed', help = 'acrE bed file')
    sp1.add_argument('--acrL', default='/home/springer/zhoux379/projects/stress/data/21_seq/15.acrL.bed', help = 'acrL bed file')
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
    sp1.add_argument('--db', default='/home/springer/zhoux379/projects/stress/data/21_seq/02.fas', help = 'database (fasta file) to search')
    sp1.add_argument('--bin', default='TSS:-/+2k', help = 'window option')
    sp1.add_argument('--epi', default='raw', choices=['raw','umr','acrE','acrL'], help = 'epigenetic filter')
    sp1.add_argument('--nfea', default='all', help = 'motif ID / option')
    sp1.add_argument('--mod', default='zoops', choices=['zoops','anr'], help = 'encoding option')
    sp1.add_argument('--umr', default='/home/springer/zhoux379/projects/stress/data/21_seq/15.umr.bed', help = 'UMR bed file')
    sp1.add_argument('--acrE', default='/home/springer/zhoux379/projects/stress/data/21_seq/15.acrE.bed', help = 'acrE bed file')
    sp1.add_argument('--acrL', default='/home/springer/zhoux379/projects/stress/data/21_seq/15.acrL.bed', help = 'acrL bed file')
    sp1.add_argument('--thresh', type=float, default=0.05, help = 'q value threshold')
    sp1.add_argument('--motif_frac', type=float, default=0.8, help = 'fraction of motif to be counted')
    sp1.add_argument('--fmt', default='wide', choices=['long','wide'], help = 'output format')
    sp1.set_defaults(func = prepare_ml)

    args = ps.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        ps.print_help()


