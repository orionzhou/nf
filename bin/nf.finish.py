#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import logging
# import numpy as np
# import pandas as pd
from subprocess import Popen, PIPE, run

# from jinja2 import Template
from jcvi.apps.base import sh, mkdir
from jcvi.formats.base import must_open

def read_samplelist(fs):
    assert op.isfile(fs), "samplelist not found: %s" % fs
    cvts = dict(SampleID=str,Tissue=str,Genotype=str)
    sl = pd.read_csv(fs, sep="\t", header=0, converters=cvts)
    return sl

def main(args):
    os.chdir(op.join(args.nfdir, args.yid))
    dry = args.dry
    if op.islink("results"):
        print("results already a link, skip moving")
    elif op.isdir("results"):
        diro = op.join(args.s3dir, args.yid)
        if op.isdir(diro):
            if dry:
                print(f"mv {diro} {diro}.bak")
            else:
                sh(f"mv {diro} {diro}.bak")
        if dry:
            print(f"mv results {diro}")
            print(f"ln -sf {diro} results")
        else:
            sh(f"mv results {diro}")
            sh(f"ln -sf {diro} results")
    if args.upload:
        diro = op.join(args.s3dir, args.yid)
        if dry:
            print(f"s3.py up zhoup-nfo/{args.yid}")
        else:
            sh(f"s3.py up zhoup-nfo/{args.yid}")

if __name__ == "__main__":
    import argparse
    ps = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = 'nextflow post-processing'
    )

    ps.add_argument('yid', help = 'study/project ID')
    libs = ['rnaseq','smrnaseq','chipseq','atacseq','methylseq','dnaseq']
    #ps.add_argument('--lib', default='rnaseq', choices=libs, help = 'library type')
    ps.add_argument('--nfdir', default='/home/springer/zhoux379/projects/barn/nf', help = 'nextflow launch dir')
    ps.add_argument('--s3dir', default='/home/springer/zhoux379/projects/s3/zhoup-nfo', help = 's3 dir')
    ps.add_argument('--dry', action='store_true', help = 'dry run')
    ps.add_argument('--upload', action='store_true', help = 'upload to amazon s3?')

    args = ps.parse_args()
    main(args)

