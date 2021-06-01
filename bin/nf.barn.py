#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import logging
import numpy as np
import pandas as pd

from jcvi.apps.base import sh, mkdir
from jcvi.formats.base import must_open

def main(args):
    yid = args.project
    fi = f'{args.excel}/{yid}.xlsx'
    if op.isfile(fi):
        sh(f"excel.py tsv {fi} tmp.tsv")
        cvts = dict(SampleID=str,Tissue=str,Genotype=str)
        sl = pd.read_csv('tmp.tsv', sep="\t", header=0, converters=cvts)
        for i in range(len(sl)):
            sid, paired = sl['SampleID'][i], sl['paired'][i]
            r0, r1, r2= sl['r0'][i], sl['r1'][i], sl['r2'][i]
            if paired == 'SE':
                f0 = op.join(args.barn, yid, r0)
                if not op.isfile(f0):
                    logging.warning(f"{f0} not found")
                sl['r0'][i] = f0
            else:
                f1 = op.join(args.barn, yid, r1)
                f2 = op.join(args.barn, yid, r2)
                if not op.isfile(f1):
                    logging.warning(f"{f1} not found")
                if not op.isfile(f2):
                    logging.warning(f"{f2} not found")
                sl['r1'][i] = f1
                sl['r2'][i] = f2
        sl.to_csv(args.output, sep="\t", header=True, index=False)
        sh("rm tmp.tsv")
    else:
        sys.exit(f"cannot read from excel [{fi}]")

if __name__ == "__main__":
    import argparse
    ps = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = "prepare sequence (fastq) meta table to run nextflow"
    )

    ps.add_argument('project', help = 'project ID')
    ps.add_argument('--output', '-o', default='design.tsv', help = 'output file')
    ps.add_argument('--excel', default=f'{os.environ["proj"]}/barn/data/05_input', help='input excel')
    ps.add_argument('--barn', default=f'{os.environ["s3"]}/zhoup-barn', help='s3 barn directory')

    args = ps.parse_args()
    main(args)

