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

def check_fastq(design, paired, source):
    sl = read_samplelist(design)
    cdic = dict(SE=0, PE=0)
    if source == 'mixed': assert 'source' in sl
    if paired == 'mixed': assert 'paired' in sl
    for i in range(len(sl)):
        sid, r0, r1, r2 = sl['SampleID'][i], sl['r0'][i], sl['r1'][i], sl['r2'][i]
        source1 = source
        if 'source' in sl:
            if source == 'mixed':
                source1 = sl['source'][i]
            else:
                assert sl['source'][i] == source
        paired1 = paired
        if 'paired' in sl:
            if paired == 'mixed':
                paired1 = sl['paired'][i]
            else:
                assert sl['paired'][i] == paired
        if paired1 == 'SE':
            if source1 == 'local':
                assert op.isfile(r0)
            elif source1 == 'sra':
                assert r0.startswith("SRR") or r0.startswith("DRR")
            elif source1 == 's3':
                assert r0.startswith("s3://")
        else:
            if source1 == 'local':
                assert op.isfile(r1)
                assert op.isfile(r2)
            elif source1 == 'sra':
                assert r0.startswith("SRR")
            elif source1 == 's3':
                assert r1.startswith("s3://")
                assert r2.startswith("s3://")

def nf_start(args):
    yid = args.yid
    barn, genome = args.metadir, args.genome

    design = "design.tsv"
    ft = "%s/%s.config" % (args.cfgdir, args.lib)
    fht = must_open(ft, 'r')
    tmp = Template(fht.read())

    aligner_map = dict(rnaseq='hisat2',chipseq='bwa',methylseq='bismark_hisat2')
    aligner = aligner_map[args.lib] if args.aligner == 'auto' else args.aligner
    msg = tmp.render(yid = yid,
                     design = design,
                     source = args.source,
                     read_type = args.read_type,
                     paired = args.paired,
                     interleaved = str(args.interleaved).lower(),
                     save_fastq = str(args.save_fastq).lower(),
                     trimmer = args.trimmer,
                     save_trimmed = str(args.save_trimmed).lower(),
                     aligner = aligner,
                     saveBAM = str(args.saveBAM).lower(),
                     genome = args.genome,
                     skip_preseq = str(not args.preseq).lower(),
                     skip_markdup = str(not args.markdup).lower(),

                     stranded = args.stranded,
                     count_multi = str(args.multi).lower(),
                     ase = str(args.ase).lower(),
                     ril = str(args.ril).lower(),
                     cage = str(args.cage).lower(),
                     salmon = str(args.salmon).lower(),
                     stringtie = str(args.stringtie).lower(),

                     narrow_peak = str(args.narrow_peak).lower()
    )

    rundir = op.join(args.projdir, args.lib, 'nf', yid)
    workdir = op.join(args.workdir, args.lib, yid)
    rawdir = op.join(args.rawdir, args.lib, yid)
    mkdir(rundir, overwrite=True)
    mkdir(workdir, overwrite=True)
    mkdir(rawdir, overwrite=True)
    os.chdir(rundir)

    # metadir = op.join(barn, genome, '05_excel')
    # xls = "%s/%s.xlsx" % (metadir, yid)
    # assert op.isfile(xls)
    sh(f"nf.barn.py {yid} -o design.tsv")
    #check_fastq(design, args.paired, args.source)

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
    allowed_sources = ['sra','local','s3','mixed']
    ps.add_argument('lib', choices = libs, help = 'library type')
    ps.add_argument('yid', help = 'study/project id')
    ps.add_argument('--projdir', default=os.environ['proj'], help = 'project dir')
    ps.add_argument('--cfgdir', default="%s/configs/templates" % os.environ['nf'], help = 'nextflow template config dir')
    ps.add_argument('--workdir', default=os.environ['NXF_WORK'], help = 'nextflow work dir')
    ps.add_argument('--rawdir', default="%s/raw" % os.environ['NXF_CACHE'], help = 'nextflow raw output dir')
    ps.add_argument('--source', default='local', choices=allowed_sources, help='sequence source')
    ps.add_argument('--read_type', default='illumina', choices=['illumina','nanopore'], help='read type')
    ps.add_argument('--paired', default='SE', choices=['SE','PE','mixed'], help='single end, paired end or mixed')
    ps.add_argument('--interleaved', action='store_true', help='interleaved format?')
    ps.add_argument('--metadir', default=os.environ['ba'], help = 'meta table directory')
    ps.add_argument('--genome', default='Zmays_B73v5', help = 'reference genome')
    ps.add_argument('--keep', action='store_true', help='keep previous results?')
    ps.add_argument('--save_fastq', action='store_true', help='save fastq files?')
    ps.add_argument('--trimmer', default='trim_galore', choices=['no','trim_galore'], help='trimming software')
    ps.add_argument('--save_trimmed', action='store_true', help='save trimmed fastq files?')
    ps.add_argument('--aligner', default='auto', help='aligning software')
    ps.add_argument('--saveBAM', action='store_true', help='save bam files?')
    ps.add_argument('--preseq', action='store_true', help='run preseq?')
    ps.add_argument('--markdup', action='store_true', help='mark PCR duplicates?')

    g1 = ps.add_argument_group('rnaseq', 'RNA-Seq specific arguments')
    g1.add_argument('--stranded', default='no', choices=['no','forward','reverse'], help = 'read strandedness')
    g1.add_argument('--multi', action='store_true', help='count multi-mapping reads?')
    g1.add_argument('--ase', action='store_true', help='allele specific expression?')
    g1.add_argument('--ril', action='store_true', help='genotype (ril) samples?')
    g1.add_argument('--cage', action='store_true', help='run CAGE pipeline?')
    g1.add_argument('--salmon', action='store_true', help='run salmon?')
    g1.add_argument('--stringtie', action='store_true', help='run stringtie?')

    g2 = ps.add_argument_group('chipseq', 'chip-seq specific arguments')
    g2.add_argument('--narrow_peak', action='store_true', help = 'turn off broad peak calling mode in MACS2')

    args = ps.parse_args()
    nf_start(args)

