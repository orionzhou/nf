#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import logging

from jcvi.apps.base import sh, mkdir
from jcvi.formats.base import must_open

def main(args):
    sid = args.sid
    if (args.source == 'sra'):
        acc = args.r0
        sh("fasterq-dump --split-files -e %d -m %dGB -O ./ -t %s %s" % (args.cpu, args.mem, args.tmp, acc))
        if (args.paired == 'PE'):
            sh("pigz -p %d --fast -c %s_1.fastq > %s_R1.fq.gz" % (args.cpu, acc, sid))
            sh("pigz -p %d --fast -c %s_2.fastq > %s_R2.fq.gz" % (args.cpu, acc, sid))
        else:
            sh("pigz -p %d --fast -c %s.fastq > %s_R0.fq.gz" % (args.cpu, acc, sid))
    elif (args.source == 'local'):
        if (args.interleaved):
            if (args.r0.endswith(".gz")):
                sh("zcat $r0 | deinterleave_fastq.sh %s_R1.fq.gz %s_R2.fq.gz %d compress" % (args.r0, sid, sid, args.cpu))
            else:
                sh("cat $r0 | deinterleave_fastq.sh %s_R1.fq.gz %s_R2.fq.gz %d compress" % (args.r0, sid, sid, args.cpu))
        else:
            if (args.paired == 'PE' and args.r1.endswith(".gz")):
                sh("ln -f %s %s_R1.fq.gz" % (args.r1, sid))
                sh("ln -f %s %s_R2.fq.gz" % (args.r2, sid))
            elif (args.paired == 'PE' and not args.r1.endswith(".gz")):
                sh("pigz -p %d -c %s > %s_R1.fq.gz" % (args.cpu, args.r1, sid))
                sh("pigz -p %d -c %s > %s_R2.fq.gz" % (args.cpu, args.r2, sid))
            elif (args.paired == 'SE' and args.r0.endswith(".gz")):
                sh("ln -f %s %s_R0.fq.gz" % (args.r0, sid))
            else:
                sh("pigz -p %d -c %s > %s_R0.fq.gz" % (args.cpu, args.r0, sid))
    elif (args.source == 's3'):
        if (args.paired == 'PE'):
            sh("s3cmd get %s %s_R1.fq.gz" % (args.r1, sid))
            sh("s3cmd get %s %s_R2.fq.gz" % (args.r2, sid))
        else:
            sh("s3cmd get %s %s_R0.fq.gz" % (args.r0, sid))

if __name__ == "__main__":
    import argparse
    ps = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = "prepare sequence (fastq) files to run workflow"
    )

    allowed_sources = ['local','sra','s3']
    ps.add_argument('sid', help = 'sample id (for output)')
    ps.add_argument('--source', default='sra', choices=allowed_sources, help='sequence source')
    ps.add_argument('--paired', default='SE', choices=['SE','PE'], help='single-end or paired-end?')
    ps.add_argument('--interleaved', action='store_true', help='interleaved format?')
    ps.add_argument('--r0', default='r0', help='read 0')
    ps.add_argument('--r1', default='r1', help='read 1')
    ps.add_argument('--r2', default='r2', help='read 2')
    ps.add_argument('--cpu', type=int, default=1, help='number processors/threads to use')
    ps.add_argument('--mem', type=int, default=20, help='size of memeroy (in GBs) to use')
    ps.add_argument('--tmp', default=os.environ['tmp'], help='temporary directory to use')

    args = ps.parse_args()
    main(args)

