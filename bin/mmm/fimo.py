#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys

from jcvi.apps.base import sh, mkdir
from jcvi.formats.base import read_block

def read_meme(fi):
    mtfs = dict()
    fhi = open(fi, 'r')
    for head,content in read_block(fhi, 'MOTIF'):
        pre, mid, mid2 = head.split(' ')
        mtf = mid.split("-")[1]
        width = len(mtf)
        mtfs[mid] = width
    return mtfs

def get_motif_width(fi):
    fhi = open(fi, 'r')
    for line in fhi:
        ps = line.strip().split("\t")
        return int(ps[2])-int(ps[1])

def fimo(args):
    fi, fo = args.fi, args.fo
    seq = args.seq
    mtfs = read_meme(fi)
    print(f'{len(mtfs)} motifs read')
    if args.motif and args.motif != 'all':
        if args.motif in mtfs:
            mtfs = {args.motif: mtfs[args.motif]}
        else:
            print("{args.motif} not in motif file {args.fi}")
            sys.exit(1)
    for mid, wd in mtfs.items():
        sh(f'fimo --motif {mid} --thresh {args.thresh} --oc fimo {fi} {seq}')
        cmd = (f'grep \'^{mid}\' fimo/fimo.tsv | '
            'bioawk -tH \'{if($9<0.05) {print $1"%"$3, $4-1, $5}}\' > x1.bed')
        sh(cmd)
        hwd = round(wd * args.motif_frac)
        if os.stat("x1.bed").st_size == 0:
            sh(f'touch fimo-{mid}.bed')
        else:
            ###wd = get_motif_width("x1.bed")
            #print(wd)
            sh(f'sortBed -i x1.bed | mergeBed > x2.bed')
            sh(f'bedtools makewindows -w {wd} -b x2.bed > x3.bed')
            sh(f'bed.py filter --minsize {hwd} x3.bed > fimo-{mid}.bed')
    sh(f'cat fimo-*.bed > {fo}')
    sh('rm -rf fimo x1.bed x2.bed x3.bed')

if __name__ == "__main__":
    import argparse
    ps = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = "run fimo and parse output"
    )

    ps.add_argument('fi', help = 'input motif file (.meme/.dreme/.streme)')
    ps.add_argument('seq', help = 'sequence file')
    ps.add_argument('fo', help = 'output file')
    ps.add_argument('--motif', default='all', help = 'motif ID')
    ps.add_argument('--thresh', type=float, default=0.05, help = 'q value threshold')
    ps.add_argument('--motif_frac', type=float, default=0.8, help = 'fraction of motif to be counted')

    args = ps.parse_args()
    fimo(args)

