#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import logging
import re

from pyfaidx import Fasta, Sequence
#from scipy.stats import fisher_exact
import pybedtools

def locate(args):
    kmer, fd, fo = args.kmer, args.db, args.out
    seqs = Fasta(fd)

    fho = open(fo, 'w')
    #fho.write('kmer\tsid\tsrd\tstart\n')
    i = 1
    for seqid in seqs.keys():
        seq = seqs[seqid][0:].seq
        kseq = kmer
        kseq2 = Sequence(name='kmer',seq=kseq).reverse.complement.seq
        for m in re.finditer(kseq, seq):
            fho.write("%s\t%d\t%s\t+\t%d\t%d\n" % (kseq, i, seqid, m.start()+1, m.start()+len(kmer)))
            i += 1
        for m in re.finditer(kseq2, seq):
            fho.write("%s\t%d\t%s\t-\t%d\t%d\n" % (kseq, i, seqid, m.start()+1, m.start()+len(kmer)))
            i += 1
    fho.close()

def read_kmer(fk, nfea):
    kms = []
    fhk = open(fk,'r')
    for line in fhk:
        line = line.rstrip("\n")
        if not line: continue
        i, opt, bin, epi, pval, fid, fname, kmers = line.split()[:8]
        if i == 'i': continue
        i = int(i)
        if nfea == 'top30' and i > 30: break
        if nfea == 'top50' and i > 50: break
        if nfea == 'top100' and i > 100: break
        kseqs = kmers.split(',')
        kseqs2 = [Sequence(name='kmer',seq=kseq).reverse.complement.seq for kseq in kseqs]
        ptn = "|".join([ "("+k+")" for k in kseqs+kseqs2 ])
        kms.append([fid,ptn,set(kseqs)])
    fhk.close()
    return kms

off = 4050
range_dict = {
    'TSS:-500': [[1500,2000]],
    'TSS:+500': [[2000,2500]],
    'TSS:-/+500':[[1500,2500]],
    'TSS:-1k':[[1000,2000]],
    'TSS:+1k':[[2000,3000]],
    'TSS:-/+1k':[[1000,3000]],
    'TSS:-2k':[[0,2000]],
    'TSS:+2k':[[2000,4000]],
    'TSS:-/+2k':[[0,4000]],
    'TTS:-500':[[off+1500,off+2000]],
    'TTS:+500':[[off+2000,off+2500]],
    'TTS:-/+500':[[off+1500,off+2500]],
    'TTS:-1k':[[off+1000,off+2000]],
    'TTS:+1k':[[off+2000,off+3000]],
    'TTS:-/+1k':[[off+1000,off+3000]],
    'TTS:-2k':[[off+0,off+2000]],
    'TTS:+2k':[[off+2000,off+4000]],
    'TTS:-/+2k':[[off+0,off+4000]],
    'TSS:-/+2k&TTS:-500':[[0,4000],[off+1500,off+2000]],
    'TSS:-/+2k&TTS:-2k':[[0,4000],[off+0,off+2000]],
    'TSS:-/+2k&TTS:-/+500':[[0,4000],[off+1500,off+2500]],
    'TSS:-/+2k&TTS:-/+1k':[[0,4000],[off+1000,off+3000]],
    'TSS:-/+2k&TTS:-/+2k':[[0,4000],[off+0,off+4000]]
}

def get_gene_region(fg, db, bin, epi, umr):
    gr1, gr2 = [], []
    fhg = open(fg,'r')
    dict_status = {}
    for line in fhg:
        line = line.rstrip("\n")
        if not line: continue
        gid, status = line.split()[:2]
        if gid == 'gid': continue
        dict_status[gid] = status
        if gid not in db: continue
        size = len(db[gid])
        gr1.append([gid, 0, size])
        rgs = range_dict[bin]
        gr2 += [[gid, s, e] for s, e in rgs]

    fhg.close()
    gr1s = "\n".join([f"{gid} {s} {e}" for gid, s, e in gr1])
    gr2s = "\n".join([f"{gid} {s} {e}" for gid, s, e in gr2])
    a = pybedtools.BedTool(gr1s, from_string=True)
    b = pybedtools.BedTool(gr2s, from_string=True)
    ab = a.intersect(b, stream=True)
    if epi == "umr":
        u = pybedtools.BedTool(umr)
        ab = ab.intersect(u)

    grd = {}
    for fea in ab:
        if fea.chrom not in grd:
            grd[fea.chrom] = []
        grd[fea.chrom].append([fea.start, fea.stop])
    grs = [[gid, dict_status[gid], grd[gid] if gid in grd else ''] for gid, x0, x1 in gr1]
    return grs

def get_gene_regions(fg, db, bins, epi, umr):
    gr1, gr2 = [], {bin: [] for bin in bins}
    fhg = open(fg,'r')
    dict_status = {}
    for line in fhg:
        line = line.rstrip("\n")
        if not line: continue
        gid, status = line.split()[:2]
        if gid == 'gid': continue
        dict_status[gid] = status
        if gid not in db: continue
        size = len(db[gid])
        gr1.append([gid, 0, size])
        for bin in bins:
            rgs = range_dict[bin]
            gr2[bin] += [[gid, s, e] for s, e in rgs]
    fhg.close()

    gr1s = "\n".join([f"{gid} {s} {e}" for gid, s, e in gr1])
    a = pybedtools.BedTool(gr1s, from_string=True)

    grd = {gid: {} for gid, x0, x1 in gr1}
    for bin in bins:
        gr2s = "\n".join([f"{gid} {s} {e}" for gid, s, e in gr2[bin]])
        b = pybedtools.BedTool(gr2s, from_string=True)
        ab = a.intersect(b, stream=True)
        if epi == "umr":
            u = pybedtools.BedTool(umr)
            ab = ab.intersect(u)
        for fea in ab:
            gid, s, e = fea.chrom, fea.start, fea.stop
            if bin not in grd[gid]:
                grd[gid][bin] = []
            grd[gid][bin].append([s, e])

    grs = [[gid, dict_status[gid], grd[gid]] for gid, x0, x1 in gr1]
    return grs

def prepare_ml(args):
    fk, fg, fd, fo, fmt = args.kmer, args.gene, args.db, args.out, args.fmt
    bin_str,epi,umr,nfea,mod = args.bin, args.epi, args.umr,args.nfea, args.mod
    bins = bin_str.split(",")
    kms = read_kmer(fk, nfea)
    db = Fasta(fd)
    grs = get_gene_regions(fg, db, bins, epi, umr)

    fho = open(fo, 'w')
    if fmt == 'wide':
        kms_str = "\t".join([e[0]+"_"+str(i) for e in kms for i in range(len(bins))])
        fho.write(f"gid\tstatus\t{kms_str}\n")
    elif fmt == 'long':
        fho.write(f"gid\tstatus\tfid\tstart\tend\tsrd\n")
    for gid, status, range_dict in grs:
        kmcs = []
        for i in range(len(bins)):
            bin = bins[i]
            kmc = {fid:0 for fid, ptn, kseqs in kms}
            if bin in range_dict:
                ranges = range_dict[bin]
                for b, e in ranges:
                    seq = db[gid][b:e].seq
                    for fid, ptn, kseqs in kms:
                        for m in re.finditer(ptn, seq):
                            if fmt == 'long':
                                start, end = m.start()+1, m.end()
                                srd = "+" if m.group(0) in kseqs else "-"
                                fho.write(f"{gid}\t{status}\t{fid}\t{start+b-1}\t{end+b-1}\t{srd}\n")
                            kmc[fid] += 1
                if mod == 'zoops':
                    kmc = {fid: min(1, cnt) for fid, cnt in kmc.items()}
            kmcs += [kmc[fid] for fid, ptn, kseqs in kms]
        if fmt == 'wide':
            kmc_str = "\t".join([str(kmc) for kmc in kmcs])
            fho.write(f"{gid}\t{status}\t{kmc_str}\n")
    fho.close()

if __name__ == "__main__":
    import argparse
    ps = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = "kmer utilities"
    )
    sp = ps.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("locate",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            help = "find given kmers in sequence database and report locations")
    sp1.add_argument('kmer', help = 'kmer sequence')
    sp1.add_argument('db', help = 'database (fasta file) to search')
    sp1.add_argument('out', help = 'output file (*.tsv)')
    #sp1.add_argument('--envdir', default='~/git/nf/configs/environments', help = 'config folder')
    sp1.set_defaults(func = locate)

    sp1 = sp.add_parser("prepare_ml",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            help = "locate given kmers for given sequences in sequence db using various filters and prepare output for ML")
    sp1.add_argument('db', help = 'database (fasta file) to search')
    sp1.add_argument('gene', help = 'gene status table')
    sp1.add_argument('kmer', help = 'kmer table')
    sp1.add_argument('out', help = 'output file (*.tsv)')
    sp1.add_argument('--bin', default='+/-2k', help = 'window option')
    sp1.add_argument('--epi', default='raw', choices=['raw','umr'], help = 'epigenetic filter')
    sp1.add_argument('--nfea', default='top30', choices=['top30','top50','top100','all'], help = 'number features/motifs to use')
    sp1.add_argument('--mod', default='zoops', choices=['zoops','anr'], help = 'encoding option')
    sp1.add_argument('--umr', default=False, help = 'UMR bed file')
    sp1.add_argument('--fmt', default='wide', choices=['long','wide'], help = 'output format')
    sp1.set_defaults(func = prepare_ml)

    args = ps.parse_args()
    args = ps.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        ps.print_help()

