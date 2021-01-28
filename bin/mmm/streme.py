#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import os.path as op
import random

import xml.etree.ElementTree as ET
from jcvi.apps.base import sh, popen, mkdir
from jcvi.formats.base import read_until, read_block, is_number

def read_streme_xml(fi):
    tree = ET.parse(fi)
    root = tree.getroot()
    mtfs = []
    for mtf in root[1]:
        d = mtf.attrib
        mtfs.append(d)
    return mtfs

def xml2tsv(args):
    mtfs = read_streme_xml(args.fi)
    print("mid","alt","wd","score_thresh","pos_train","neg_train","pos_test","neg_test",sep="\t")
    for d in mtfs:
        print(d['id'], d['alt'], d['width'], d['score_threshold'],
              d['train_pos_count'], d['train_neg_count'],
              d['test_pos_count'], d['test_neg_count'], sep="\t")

def add_score(args):
    mtfs = read_streme_xml(args.fx)
    fhi = open(args.fi, 'r')
    fho = open(args.fo, 'w')
    #
    start = "MOTIF"
    while 1:
        pos = fhi.tell()
        line = fhi.readline()
        if not line:
            break
        if line.startswith(start):
            fhi.seek(pos)
            break
        else:
            fho.write(line)
    #
    i = 0
    for head, content in read_block(fhi, 'MOTIF'):
        pre, mid, alt = head.split(' ')
        assert mtfs[i]['id'] == mid, 'motifs not in sync'
        new_head = f"{pre} {mid} {mtfs[i]['score_threshold']}"
        fho.write(new_head + "\n")
        for line in content:
            if line.startswith("0") or line.startswith("1"):
                fho.write(" " + line + "\n")
            else:
                fho.write(line + "\n")
        i += 1
    fhi.close()
    fho.close()

def pipe(args):
    fi, out = args.fi, args.out
    neg_str = f"-n {args.neg}" if args.neg else ""
    pre = f"tmp.sr{random.randrange(1000)}"
    sec = args.time * 3600
    cmd = f"streme --p {fi} {neg_str} --dna --pvt {args.pval} --time {sec} -minw {args.minw} -maxw {args.maxw} -oc {pre}"
    sh(cmd)
    sh(f"streme.py xml2tsv {pre}/streme.xml > {out}.tsv")
    sh(f"streme.py addscore {pre}/streme.txt {pre}/streme.xml {out}.dreme")
    sh(f"fimo.py locate {out}.dreme {fi} {out}.bed")
    sh(f"rm -rf {pre}")

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = 'STREME utilities'
    )
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("addscore", help = "add score_thresh to STREME output")
    sp1.add_argument('fi', help = 'STREME motif file')
    sp1.add_argument('fx', help = 'STREME xml file')
    sp1.add_argument('fo', help = 'output motif file')
    sp1.set_defaults(func = add_score)

    sp1 = sp.add_parser("xml2tsv", help = "convert STREME xml output to tsv")
    sp1.add_argument('fi', help = 'STREME xml file')
    #sp1.add_argument('fo', help = 'output tsv file')
    sp1.set_defaults(func = xml2tsv)

    sp1 = sp.add_parser("pipe", help = "run STREME pipeline",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    sp1.add_argument('fi', help = 'input sequences (positive)')
    sp1.add_argument('out', help = 'output prefix')
    sp1.add_argument('--minw', default=6, help = 'minimum motif width')
    sp1.add_argument('--maxw', default=20, help = 'maximum motif width')
    sp1.add_argument('--neg', default=None, help = 'control (negtive) sequences')
    sp1.add_argument('--pval', default=1e-2, help = 'p value')
    sp1.add_argument('--time', default=20, help = 'maximum running time in hours')
    sp1.set_defaults(func = pipe)

    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()
