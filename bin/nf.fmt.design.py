#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import logging
import numpy as np
import pandas as pd

def read_samplelist(fs):
    assert op.isfile(fs), "samplelist not found: %s" % fs
    cvts = dict(SampleID=str,Tissue=str,Genotype=str,group=str,control=str,antibody=str)
    sl = pd.read_csv(fs, sep="\t", header=0, converters=cvts)
    return sl

def fmt_chipseq_1(args):
    yid = op.splitext(op.basename(args.fi))[0]
    sl = read_samplelist(args.fi)
    HEADER = ['group', 'replicate', 'fastq_1', 'fastq_2', 'antibody', 'control', 'spots']
    print(",".join(HEADER))
    for i in range(len(sl)):
        sid, paired = sl['SampleID'][i], sl['paired'][i]
        group, rep, ab, ctrl = sl['group'][i], sl['Replicate'][i], sl['antibody'][i], sl['control'][i]
        spots = sl['spots'][i]
        fq1, fq2 = '', ''
        if paired:
            f1 = "%s/%s/%s_1.fq.gz" % (args.seqdir, yid, sid)
            f1b = "%s/%s/%s_R1.fq.gz" % (args.seqdir, yid, sid)
            f2 = "%s/%s/%s_2.fq.gz" % (args.seqdir, yid, sid)
            f2b = "%s/%s/%s_R2.fq.gz" % (args.seqdir, yid, sid)
            assert op.isfile(f1) or op.isfile(f1b), "fastq not found: %s" % f1
            assert op.isfile(f2) or op.isfile(f2b), "fastq not found: %s" % f2
            fq1 = f1 if op.isfile(f1) else f1b
            fq2 = f2 if op.isfile(f2) else f2b
        else:
            f1 = "%s/%s/%s.fq.gz" % (args.seqdir, yid, sid)
            f1b = "%s/%s/%s_R0.fq.gz" % (args.seqdir, yid, sid)
            assert op.isfile(f1) or op.isfile(f1b), "fastq not found: %s" % f1
            fq1 = f1 if op.isfile(f1) else f1b
        print(",".join([group, str(rep), fq1, fq2, ab, ctrl, str(spots)]))
    yid = op.basename(args.fi)

def fmt_chipseq(args):
    fi, fo, fm = args.fi, args.fo, args.fm
    ERROR_STR = 'ERROR: Please check design file'
    yid = op.splitext(op.basename(args.fi))[0]
    sl = read_samplelist(args.fi)

    sampleMappingDict, antibodyDict = {}, {}
    for i in range(len(sl)):
        sid, paired = sl['SampleID'][i], sl['paired'][i]
        group, replicate, antibody, control = sl['group'][i], sl['Replicate'][i], sl['antibody'][i], sl['control'][i]
        ## CHECK GROUP ID DOESNT CONTAIN SPACES
        if group.find(' ') != -1:
            print("{}: Group id contains spaces!\nLine: '{}'".format(ERROR_STR,line.strip()))
            sys.exit(1)

        ## CREATE GROUP MAPPING DICT = {GROUP_ID: {REPLICATE_ID: [SAMPLE_ID(s)]}
        if not group in sampleMappingDict:
            sampleMappingDict[group] = {}
        if not replicate in sampleMappingDict[group]:
            sampleMappingDict[group][replicate] = []
        sampleMappingDict[group][replicate].append([i, sid])

        ## CHECK BOTH ANTIBODY AND CONTROL COLUMNS HAVE VALID VALUES
        if antibody:
            if antibody.find(' ') != -1:
                print("{}: Antibody id contains spaces!\nLine: '{}'".format(ERROR_STR,line.strip()))
                sys.exit(1)
            if not control:
                print("{}: both Antibody and Control must be specified!\nLine: '{}'".format(ERROR_STR,line.strip()))
                sys.exit(1)
        if control:
            if control.find(' ') != -1:
                print("{}: Control id contains spaces!\nLine: '{}'".format(ERROR_STR,line.strip()))
                sys.exit(1)
            if not antibody:
                print("{}: both Antibody and Control must be specified!\nLine: '{}'".format(ERROR_STR,line.strip()))
                sys.exit(1)

        ## CREATE ANTIBODY MAPPING CONTROL DICT
        if antibody and control:
            antibodyDict[group] = (antibody,control)

    ## CHECK IF DATA IS PAIRED-END OR SINGLE-END AND NOT A MIXTURE

    ## CHECK IF ANTIBODY AND CONTROL COLUMNS HAVE BEEN SPECIFIED AT LEAST ONCE
    if len(antibodyDict) == 0:
        print("{}: Antibody and Control must be specified at least once!".format(ERROR_STR))
        sys.exit(1)

    ## WRITE READ MAPPING FILE
    antibodyGroupDict = {}
    sl.assign(acc = [''] * len(sl))
    for group in sorted(sampleMappingDict.keys()):
        ## CHECK THAT REPLICATE IDS ARE IN FORMAT 1..<NUM_REPLICATES>
        uniq_rep_ids = set(sampleMappingDict[group].keys())
        if len(uniq_rep_ids) != max(uniq_rep_ids):
            print("{}: Replicate IDs must start with 1..<num_replicates>\nGroup: {}, Replicate IDs: {}".format(ERROR_STR,group,list(uniq_rep_ids)))
            sys.exit(1)

        ## RECONSTRUCT LINE FOR SAMPLE IN DESIGN
        for replicate in sorted(sampleMappingDict[group].keys()):
            for idx in range(len(sampleMappingDict[group][replicate])):
                i, sid = sampleMappingDict[group][replicate][idx]

                ## GET SAMPLE_ID,FASTQ_1,FASTQ_2 COLUMNS
                sample_id = "{}_R{}_T{}".format(group,replicate,idx+1)
                sl.at[i, 'SampleID'] = sample_id
                sl.at[i, 'acc'] = sid

                ## EXTRAPOLATE CONTROL COLUMN
                if group in antibodyDict:
                    antibody,control = antibodyDict[group]
                    if control in sampleMappingDict:
                        control_id = "{}_R1".format(control)
                        if replicate in sampleMappingDict[control]:
                            control_id = "{}_R{}".format(control,replicate)
                        if not antibody in antibodyGroupDict:
                            antibodyGroupDict[antibody] = {}
                        if not group in antibodyGroupDict[antibody]:
                            antibodyGroupDict[antibody][group] = []
                        antibodyList = [sample_id[:-3],control_id]
                        if not antibodyList in antibodyGroupDict[antibody][group]:
                            antibodyGroupDict[antibody][group].append(antibodyList)
                    else:
                        print("{}: Control id not a valid group\nControl id: {}, Valid Groups: {}".format(ERROR_STR,control,sorted(sampleMappingDict.keys())))
                        sys.exit(1)
    #fho = open(fo,'w')
    #fho.close()
    sl.to_csv(fo, index=False, sep="\t")

    ## WRITE SAMPLE TO CONTROL MAPPING FILE
    fhm = open(fm,'w')
    fhm.write('\t'.join(['sample_id','control_id','antibody','replicatesExist','multipleGroups']) + '\n')
    for antibody in sorted(antibodyGroupDict.keys()):
        repsExist = '0'
        if max([len(x) for x in antibodyGroupDict[antibody].values()]) > 1:
            repsExist = '1'
        multipleGroups = '0'
        if len(antibodyGroupDict[antibody].keys()) > 1:
            multipleGroups = '1'
        for group in sorted(antibodyGroupDict[antibody].keys()):
            for antibodyList in antibodyGroupDict[antibody][group]:
                fhm.write('\t'.join(antibodyList+[antibody,repsExist,multipleGroups]) + '\n')
    fhm.close()

def fmt_atacseq(args):
    yid = op.splitext(op.basename(args.fi))[0]
    sl = read_samplelist(args.fi)
    HEADER = ['group', 'replicate', 'fastq_1', 'fastq_2', 'spots']
    print(",".join(HEADER))
    for i in range(len(sl)):
        sid, paired = sl['SampleID'][i], sl['paired'][i]
        group, rep = sl['group'][i], sl['Replicate'][i]
        spots = sl['spots'][i]
        fq1, fq2 = '', ''
        if paired:
            f1 = "%s/%s/%s_1.fq.gz" % (args.seqdir, yid, sid)
            f1b = "%s/%s/%s_R1.fq.gz" % (args.seqdir, yid, sid)
            f2 = "%s/%s/%s_2.fq.gz" % (args.seqdir, yid, sid)
            f2b = "%s/%s/%s_R2.fq.gz" % (args.seqdir, yid, sid)
            assert op.isfile(f1) or op.isfile(f1b), "fastq not found: %s" % f1
            assert op.isfile(f2) or op.isfile(f2b), "fastq not found: %s" % f2
            fq1 = f1 if op.isfile(f1) else f1b
            fq2 = f2 if op.isfile(f2) else f2b
        else:
            f1 = "%s/%s/%s.fq.gz" % (args.seqdir, yid, sid)
            f1b = "%s/%s/%s_R0.fq.gz" % (args.seqdir, yid, sid)
            assert op.isfile(f1) or op.isfile(f1b), "fastq not found: %s" % f1
            fq1 = f1 if op.isfile(f1) else f1b
        print(",".join([group, str(rep), fq1, fq2, str(spots)]))


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'nextflow format meta/design table'
    )
    sp = parser.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("chipseq",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            help = "format chipseq design table")
    #sp1.add_argument('yid', help = 'study/project ID')
    sp1.add_argument('fi', help = 'design/meta file')
    sp1.add_argument('fo', help = 'output design/meta file')
    sp1.add_argument('fm', help = 'control mapping file')
    #sp1.add_argument('--metadir', default='/home/springer/zhoux379/projects/barn/data/15_read_list', help = 'meta dir')
    #sp1.add_argument('--seqdir', default='/scratch.global/zhoux379/barn/data/fastq', help = 'seq dir')
    sp1.set_defaults(func = fmt_chipseq)

    sp1 = sp.add_parser("atacseq",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            help = "format atacseq design table")
    sp1.add_argument('fi', help = 'design/meta file')
    sp1.add_argument('fo', help = 'output design/meta file')
    sp1.set_defaults(func = fmt_atacseq)

    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        parser.print_help()

