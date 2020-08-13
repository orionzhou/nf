#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import logging

from jcvi.apps.base import sh, mkdir
from jcvi.formats.base import must_open

def init(args):
    assert op.isfile(fs), "samplelist not found: %s" % fs
    cvts = dict(SampleID=str,Tissue=str,Genotype=str)
    sl = pd.read_csv(fs, sep="\t", header=0, converters=cvts)
    return sl

if __name__ == "__main__":
    import argparse
    ps = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = 'help utilities to initialize/clean conda environments'
    )
    sp = ps.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("init",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            help = "initialize conda environments using *.yml configurations")
    sp1.add_argument('--envdir', default='~/git/nf/configs/environments', help = 'config folder')
    sp1.set_defaults(func = init)

    args = ps.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        ps.print_help()

