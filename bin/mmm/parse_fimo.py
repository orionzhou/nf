#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys

from jcvi.apps.base import sh, mkdir

def parse_fimo(args):
    fi, fo = args.fi, args.fo
    if os.stat(fi).st_size > 0:
        cmd = """sed 1d %s | grep '^[a-zA-Z]' |\
            bioawk -tH '{if($7 > 0) {print $3, $4-1, $5, ".", $7, $6, $8}}' |\
            sortBed | mergeBed -s -c 5,6,7 -o max,distinct,min |\
            bioawk -tH '{split($1,a,/%%/); print a[1], 1, 2, $3-$2, $4, $6}' |\
            mergeBed -c 4,4,5,6 -o count,sum,sum,min |\
            cut -f1,4-7 > %s""" % (fi, fo)
        sh(cmd)
    else:
        sh("touch %s" % fo)

if __name__ == "__main__":
    import argparse
    ps = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = "parse FIMI output"
    )

    ps.add_argument('fi', help = 'input file')
    ps.add_argument('fo', help = 'output file')

    args = ps.parse_args()
    parse_fimo(args)

