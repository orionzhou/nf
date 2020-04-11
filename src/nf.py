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

def nf_start(args):
    for yid in args.yid.split(","):
        print("working on %s" % yid)
        nf_start_one(yid, args)

def nf_start_one_2(yid, args):
    fi = "%s/%s.tsv" % (args.metadir, yid)
    sl = read_samplelist(fi)
    cols = sl.columns.values.tolist()
    cnts = sl["paired"].value_counts(sort=True, dropna=True)
    paired, cnt = cnts.index[0], cnts.values[0]
    if cnts.size > 1:
        cnts2 = cnts.sort_index()
        print("    warning: mix of single-end [%d] and paired-end [%d] reads" % (cnts2.values[0], cnts2.values[1]))
        if args.pair != 'auto':
            paired = True if args.pair == 'paired' else False
            cnt = cnts2.values[1] if args.pair == 'paired' else cnts2.values[0]
    paired_str = "paired" if paired else "single"
    print("    proceeding with %d %s-end reads" % (cnt, paired_str))

    if args.lib == 'smrnaseq' and paired:
        print("    %s reads not supported for lib[%s]" % (paired_str, args.lib))
        sys.exit(1)
    pre = yid
    if args.pair != 'auto':
        pre1 = 'pe' if args.pair == 'paired' else 'se'
        pre = "%s_%s" % (pre, pre1)

    fo = "%s/%s.tsv" % (args.dsgdir, pre)
    fho = must_open(fo, 'w')
    if paired:
        fho.write("SampleID\tGenotype\tfastq1\tfastq2\n")
    else:
        fho.write("SampleID\tGenotype\tfastq\n")
    for i in range(len(sl)):
        sid, paired0 = sl['SampleID'][i], sl['paired'][i]
        gt = sl['Genotype'][i]
        if paired0 != paired: continue
        if args.lib in ['chipseq', 'dapseq']:
            print("libtype[%s] not supported yet" % lib)
            sys.exit(1)
        if paired:
            f1 = "%s/%s/%s_1.fq.gz" % (args.seqdir, yid, sid)
            f2 = "%s/%s/%s_2.fq.gz" % (args.seqdir, yid, sid)
            assert op.isfile(f1), "fastq not found: %s" % f1
            assert op.isfile(f2), "fastq not found: %s" % f2
            fho.write("\t".join([sid, gt, f1, f2]) + '\n')
        else:
            f1 = "%s/%s/%s.fq.gz" % (args.seqdir, yid, sid)
            assert op.isfile(f1), "fastq not found: %s" % f1
            fho.write("\t".join([sid, gt, f1]) + '\n')
    fho.close()

    ft = "%s/tmpl.%s.nf" % (args.cfgdir, args.lib)
    fht = must_open(ft, 'r')
    tmp = Template(fht.read())
    msg = tmp.render(yid = pre, singleEnd = str(not paired).lower(),
                     ase = str(args.ase).lower(),
                     ril = str(args.ril).lower())
    fc = "%s/%s.nf" % (args.cfgdir, pre)
    fhc = must_open(fc, 'w')
    fhc.write(msg)
    fhc.close()

    os.chdir(args.rundir)
    if not op.isdir(pre):
        sh("rm -f %s" % pre)
        sh("mkdir %s" % pre)
    os.chdir(pre)
    fn = 'nextflow.config'
    if op.exists(fn):
        sh("rm -rf %s" % fn)
    sh("ln -sf %s %s" % (fc, fn))

def nf_start_one(yid, args):
    fi = "%s/%s.tsv" % (args.metadir, yid)
    sl = read_samplelist(fi)
    for i in range(len(sl)):
        sid, paired = sl['SampleID'][i], sl['paired'][i]
        if paired:
            f1 = "%s/%s/%s_1.fq.gz" % (args.seqdir, yid, sid)
            f2 = "%s/%s/%s_2.fq.gz" % (args.seqdir, yid, sid)
            assert op.isfile(f1), "fastq not found: %s" % f1
            assert op.isfile(f2), "fastq not found: %s" % f2
        else:
            f1 = "%s/%s/%s.fq.gz" % (args.seqdir, yid, sid)
            assert op.isfile(f1), "fastq not found: %s" % f1

    ft = "%s/tmpl.%s.nf" % (args.cfgdir, args.lib)
    fht = must_open(ft, 'r')
    tmp = Template(fht.read())
    srd = str(args.strand).lower()
    if srd != 'false': srd = "'%s'" % srd
    msg = tmp.render(yid = yid,
                     strand = srd,
                     ase = str(args.ase).lower(),
                     ril = str(args.ril).lower())
    fc = "%s/%s.nf" % (args.cfgdir, yid)
    fhc = must_open(fc, 'w')
    fhc.write(msg)
    fhc.close()

    os.chdir(args.rundir)
    if not op.isdir(yid):
        sh("rm -f %s" % yid)
        sh("mkdir %s" % yid)
    os.chdir(yid)
    fn = 'nextflow.config'
    if op.exists(fn):
        sh("rm -rf %s" % fn)
    sh("ln -sf %s %s" % (fc, fn))

    if not args.keep:
        sh("rm -rf %s/%s" % (args.rawdir, yid))

def nf_publish(args):
    for yid in args.yid.split(","):
        print("working on %s" % yid)
        nf_publish_one(yid, args)

def nf_publish_one(yid, args):
    diri1, fi1 = '', ''
    if args.lib in ['rnaseq']:
        diri1 = "%s/%s/MultiQC" % (args.rawdir, yid)
        fi1 = "%s/%s_multiqc_report.html" % (diri1, yid)
    elif args.lib in ['chipseq']:
        diri1 = "%s/%s/multiqc/broadPeak" % (args.rawdir, yid)
        fi1 = "%s/%s_multiqc_report.html" % (diri1, yid)
    else:
        print("libtype[%s] not supported yet" % lib)
        sys.exit(1)
    assert op.isdir(diri1), "multiqc dir not found: %s" % diri1
    assert op.isfile(fi1), "multiqc file not found: %s" % fi1
    diro1 = "%s/%s" % (args.mqcdir, yid)
    if op.exists(diro1):
        sh("rm -rf %s" % diro1)
    sh("cp -Rf %s %s" % (diri1, diro1))
    diro2 = "%s/%s" % (args.outdir, yid)
    if not op.isdir(diro2):
        sh("rm -f %s" % diro2)
        sh("mkdir -p %s" % diro2)
    fo1 = "%s/%s.html" % (args.webdir, yid)
    if op.exists(fo1):
        sh("rm -rf %s" % fo1)
    size_mqc = op.getsize(fi1) / 1000000
    print("    multiqc.html: %.1f Mb" % size_mqc)
    if size_mqc >= 100:
        print("    multiqc.html too large to publish")
    else:
        sh("cp -f %s %s" % (fi1, fo1))

def nf_clean(args):
    os.chdir(args.workdir)
    for yid in args.yid.split(","):
        print("cleanup %s workdir" % yid)
        sh("rm -rf %s" % yid)
    if args.fq:
        os.chdir(args.seqdir)
        for yid in args.yid.split(","):
            print("cleanup %s seqdir" % yid)
            sh("rm -rf %s" % yid)

def cpnf(args):
    tags = args.tags.split(",")
    print("will copy %d project results: %s" % (len(tags), ' '.join(tags)))
    for tag in tags:
        ds = op.join(args.srcdir, tag)
        dd = op.join(args.destdir, tag)
        if op.isdir(ds):
            print("%s : copying" % tag)
            sh("find %s -not \\( -name \"*.bam\" -o -name \"*.bigwig\" \\) | cpio -pdm %s" % (ds, dd))
        else:
            print("%s : not exist - skipped")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'nextflow helper utilities'
    )
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("start",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            help = "initiate one or more nextflow projects (creating design and config files)")
    sp1.add_argument('yid', help = 'study/project ID(s) separated by comma')
    libs = ['rnaseq','smrnaseq','chipseq','dapseq','atacseq','methylseq','dnaseq']
    sp1.add_argument('--lib', default='rnaseq', choices=libs, help = 'library type')
    sp1.add_argument('--metadir', default='/home/springer/zhoux379/projects/barn/data/15_read_list', help = 'meta dir')
    sp1.add_argument('--seqdir', default='/scratch.global/zhoux379/barn/data/fastq', help = 'seq dir')
    sp1.add_argument('--rundir', default='/home/springer/zhoux379/projects/nf/run', help = 'nextflow run dir')
    sp1.add_argument('--dsgdir', default='/home/springer/zhoux379/projects/nf/design', help = 'nextflow design/input dir')
    sp1.add_argument('--cfgdir', default='/home/springer/zhoux379/projects/nf/configs/projects', help = 'nextflow config dir')
    sp1.add_argument('--rawdir', default='/home/springer/zhoux379/projects/nf/raw', help = 'nextflow raw dir')
    sp1.add_argument('--strand', default='false', help = 'read strandedness')
    sp1.add_argument('--pair', default='auto', choices=['auto','single','paired'], help = 'force specify paired end option')
    sp1.add_argument('--ase', action='store_true', help = 'allele specific expression?')
    sp1.add_argument('--ril', action='store_true', help = 'genotype (RIL) samples?')
    sp1.add_argument('--keep', action='store_true', help = 'keep previous results?')
    sp1.set_defaults(func = nf_start)

    sp1 = sp.add_parser("publish",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            help = "publish nextflow results (multiqc files)")
    sp1.add_argument('yid', help = 'study/project ID(s) separated by comma')
    libs = ['rnaseq','smrnaseq','chipseq','dapseq','atacseq','methylseq','dnaseq']
    sp1.add_argument('--lib', default='rnaseq', choices=libs, help = 'library type')
    sp1.add_argument('--rawdir', default='/home/springer/zhoux379/projects/nf/raw', help = 'nextflow raw dir')
    sp1.add_argument('--outdir', default='/home/springer/zhoux379/projects/nf/out', help = 'output dir')
    sp1.add_argument('--mqcdir', default='/home/springer/zhoux379/projects/nf/multiqc', help = 'multiqc dir')
    sp1.add_argument('--webdir', default='/home/springer/zhoux379/git/orionzhou.github.io/public/multiqc', help = 'website dir')
    sp1.set_defaults(func = nf_publish)

    sp1 = sp.add_parser("clean",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            help = "clean nextflow work directory")
    sp1.add_argument('yid', help = 'study/project ID(s) separated by comma')
    sp1.add_argument('--workdir', default='/home/springer/zhoux379/projects/nf/work', help = 'nextflow work dir')
    sp1.add_argument('--seqdir', default='/scratch.global/zhoux379/barn/data/fastq', help = 'seq dir')
    sp1.add_argument('--fq', action='store_true', help = 'remove fastq')
    sp1.set_defaults(func = nf_clean)

    sp1 = sp.add_parser("cpnf",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            help = "copy nextflow result files")
    sp1.add_argument('tags', help = 'sub-dirs to copy')
    #sp1.add_argument('--srcdir', default='/scratch.global/zhoux379/nf/raw', help = 'source dir')
    sp1.add_argument('--srcdir', default='/home/springer/zhoux379/projects/nf/data/raw', help = 'source dir')
    sp1.add_argument('--destdir', default='/home/springer/zhoux379/projects/nf/data/out', help = 'destination dir')
    sp1.set_defaults(func = cpnf)

    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()

