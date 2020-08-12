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
    sl = pd.read_csv(fs, sep="\t", header=0, converters=cvts,
                     true_values=['1','Y','Yes','T','True'],
                     false_values=['0','N','No','F','False'])
    return sl

def check_fastq_0(args, yid):
    fi = "%s/%s.tsv" % (args.metadir, yid)
    sl = read_samplelist(fi)
    cdic = dict(single=0, pair=0)
    for i in range(len(sl)):
        sid, paired = sl['SampleID'][i], sl['paired'][i]
        if paired:
            f1 = "%s/%s/%s_1.fq.gz" % (args.seqdir, yid, sid)
            f1b = "%s/%s/%s_R1.fq.gz" % (args.seqdir, yid, sid)
            f2 = "%s/%s/%s_2.fq.gz" % (args.seqdir, yid, sid)
            f2b = "%s/%s/%s_R2.fq.gz" % (args.seqdir, yid, sid)
            assert op.isfile(f1) or op.isfile(f1b), "fastq not found: %s" % f1
            assert op.isfile(f2) or op.isfile(f2b), "fastq not found: %s" % f2
            cdic["pair"] += 1
        else:
            f1 = "%s/%s/%s.fq.gz" % (args.seqdir, yid, sid)
            f1b = "%s/%s/%s_R0.fq.gz" % (args.seqdir, yid, sid)
            assert op.isfile(f1) or op.isfile(f1b), "fastq not found: %s" % f1
            cdic["single"] += 1
    fcdic = { k: v for k, v in cdic.items() if v > 0}
    single_end = False
    if len(fcdic) > 1:
        print("mixed library: %d single + %d paired" % (fcdir['single'], fcdir['pair']))
    else:
        if list(fcdic.keys())[0] == 'single': single_end = True
    return single_end

def check_fastq(design):
    sl = read_samplelist(design)
    cdic = dict(single=0, pair=0)
    for i in range(len(sl)):
        sid, paired = sl['SampleID'][i], sl['paired'][i]
        if paired:
            cdic["pair"] += 1
        else:
            cdic["single"] += 1
    fcdic = { k: v for k, v in cdic.items() if v > 0}
    single_end = False
    if len(fcdic) > 1:
        print("warning: mixed library: %d single + %d paired" % (fcdir['single'], fcdir['pair']))
    else:
        if list(fcdic.keys())[0] == 'single': single_end = True
    return single_end

def nf_start(args):
    yid = args.yid
    metadir_x, metadir = args.metadir_sx, args.metadir_s
    if args.source == "local":
        metadir_x, metadir = args.metadir_lx, args.metadir_l
    xls = "%s/%s.xlsx" % (metadir_x, yid)
    design = "%s/%s.tsv" % (metadir, yid)
    sh("excel.py tsv %s %s" % (xls, design))
    single_end = check_fastq(design)

    ft = "%s/%s.config" % (args.cfgdir, args.lib)
    fht = must_open(ft, 'r')
    tmp = Template(fht.read())
    srd = str(args.strand).lower()
    if srd != 'false': srd = "'%s'" % srd

    aligner_map = dict(rnaseq='hisat2',chipseq='bwa',methylseq='bismark_hisat2')
    aligner = aligner_map[args.lib] if args.aligner == 'auto' else args.aligner
    msg = tmp.render(yid = yid,
                     design = design,
                     source = args.source,
                     interleaved = str(args.interleaved).lower(),
                     save_fastq = str(args.save_fastq).lower(),
                     save_trimmed = str(args.save_trimmed).lower(),
                     aligner = aligner,
                     saveBAM = str(args.saveBAM).lower(),
                     genome = args.genome,
                     skip_preseq = str(not args.preseq).lower(),

                     strand = srd,
                     ase = str(args.ase).lower(),
                     ril = str(args.ril).lower(),
                     salmon = str(args.salmon).lower(),
                     stringtie = str(args.stringtie).lower(),

                     single_end = str(single_end).lower(),
                     narrow_peak = str(args.narrow_peak).lower()
    )

    rundir = op.join(args.projdir, args.lib, 'nf', yid)
    workdir = op.join(args.workdir, args.lib, yid)
    rawdir = op.join(args.rawdir, args.lib, yid)
    mkdir(rundir, overwrite=True)
    mkdir(workdir, overwrite=True)
    mkdir(rawdir, overwrite=True)
    os.chdir(rundir)

    fc = "nextflow.config"
    fhc = must_open(fc, 'w')
    fhc.write(msg)
    fhc.close()

    sh("ln -sf %s/ work" % workdir)
    sh("ln -sf %s/ raw" % rawdir)

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
    ps = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = "initiate a nextflow project (creating design and config files)"
    )

    libs = ['rnaseq','smrnaseq','chipseq','dapseq','atacseq','methylseq','dnaseq']
    ps.add_argument('lib', choices = libs, help = 'library type')
    ps.add_argument('yid', help = 'study/project id')
    ps.add_argument('--projdir', default='/home/springer/zhoux379/projects', help = 'project dir')
    ps.add_argument('--cfgdir', default='/home/springer/zhoux379/git/nf/configs/templates', help = 'nextflow template config dir')
    ps.add_argument('--workdir', default='/scratch.global/zhoux379/nf/work', help = 'nextflow work dir')
    ps.add_argument('--rawdir', default='/scratch.global/zhoux379/nf/raw', help = 'nextflow raw output dir')
    ps.add_argument('--source', default='sra', choices=['local','sra','sra2'], help='sequence source')
    ps.add_argument('--interleaved', action='store_true', help='sequence source')
    ps.add_argument('--metadir_lx', default='/home/springer/zhoux379/projects/barn/data/06_local_list_excel', help = 'local meta excels')
    ps.add_argument('--metadir_l', default='/home/springer/zhoux379/projects/barn/data/07_local_list', help = 'local meta dir')
    ps.add_argument('--metadir_sx', default='/home/springer/zhoux379/projects/barn/data/08_sra_list_excel', help = 'sra meta excels')
    ps.add_argument('--metadir_s', default='/home/springer/zhoux379/projects/barn/data/09_sra_list', help = 'sra meta dir')
#    ps.add_argument('--seqdir', default='/scratch.global/zhoux379/barn/data/fastq', help = 'seq dir')
    ps.add_argument('--keep', action='store_true', help='keep previous results?')
    ps.add_argument('--save_fastq', action='store_true', help='save fastq files?')
    ps.add_argument('--save_trimmed', action='store_true', help='save trimmed fastq files?')
    ps.add_argument('--genome', default='Zmays_B73', help = 'reference genome')
    ps.add_argument('--aligner', default='auto', help='aligning software')
    ps.add_argument('--saveBAM', action='store_true', help='save bam files?')
    ps.add_argument('--preseq', action='store_true', help='run preseq?')

    g1 = ps.add_argument_group('rnaseq', 'rna-seq specific arguments')
    g1.add_argument('--strand', default='false', help = 'read strandedness')
    g1.add_argument('--ase', action='store_true', help='allele specific expression?')
    g1.add_argument('--ril', action='store_true', help='genotype (ril) samples?')
    g1.add_argument('--salmon', action='store_true', help='run salmon?')
    g1.add_argument('--stringtie', action='store_true', help='run stringtie?')

    g2 = ps.add_argument_group('chipseq', 'chip-seq specific arguments')
    #g2.add_argument('--pair', default='auto', choices=['auto','single','paired'], help = 'force specify paired end option')
    g2.add_argument('--narrow_peak', action='store_true', help = 'turn off broad peak calling mode in MACS2')

    args = ps.parse_args()
    nf_start(args)

