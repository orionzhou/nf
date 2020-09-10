#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import logging

from jcvi.apps.base import sh, mkdir
from jcvi.formats.base import must_open

def init(args):
    for env in args.envs.split(","):
        fe = "%s/%s.yml" % (args.envdir, env)
        assert op.isfile(fe), "env file found: %s" % fe
        cmd = "conda env create -n %s -f %s" % (env, fe)
        sh(cmd)
        cmd = "conda update -n %s --all" % (env)
        sh(cmd)

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
    sp1.add_argument('envs', default="genome", help = 'comma separated list of environment names')
    sp1.add_argument('--envdir', default="%s/configs/environments" % os.environ["nf"], help = 'config folder')
    sp1.add_argument('--update', action="store_true", help = 'auto-update?')
    sp1.set_defaults(func = init)

    args = ps.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        ps.print_help()

