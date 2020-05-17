#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import logging
import numpy as np
import pandas as pd
from subprocess import Popen, PIPE, run

from jinja2 import Template
from jcvi.apps.base import sh, mkdir
from jcvi.formats.base import must_open

def read_samplelist(fs):
    assert op.isfile(fs), "samplelist not found: %s" % fs
    cvts = dict(SampleID=str,Tissue=str,Genotype=str)
    sl = pd.read_csv(fs, sep="\t", header=0, converters=cvts)
    return sl

def main(args):
    os.chdir(op.join(args.workdir, args.lib))
    for yid in args.yid.split(","):
        print("cleanup %s workdir" % yid)
        sh("rm -rf %s" % yid)
    if args.fq:
        os.chdir(args.seqdir)
        for yid in args.yid.split(","):
            print("cleanup %s seqdir" % yid)
            sh("rm -rf %s" % yid)

if __name__ == "__main__":
    import argparse
    ps = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = 'cleanup nextflow work directories'
    )

    ps.add_argument('yid', help = 'study/project ID(s) separated by comma')
    libs = ['rnaseq','smrnaseq','chipseq','atacseq','methylseq','dnaseq']
    ps.add_argument('--lib', default='rnaseq', choices=libs, help = 'library type')
    ps.add_argument('--workdir', default='/scratch.global/zhoux379/nf/work', help = 'nextflow work dir')
    ps.add_argument('--seqdir', default='/scratch.global/zhoux379/barn/data/fastq', help = 'seq dir')
    ps.add_argument('--fq', action='store_true', help = 'remove fastq')

    args = ps.parse_args()
    main(args)

