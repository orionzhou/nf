#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import pandas as pd
#import simplejson as json
import yaml

def get_gsize(fs):
    cl = pd.read_csv(fs, sep="\t", header=None, names=['chrom','size'])
    return sum(cl['size'])

def dump(args):
    cvts = dict(genome=str,species=str,source=str)
    gl = pd.read_csv(args.fo, sep="\t", header=0, converters=cvts,
                     true_values=['1','Y','Yes','T','True'],
                     false_values=['0','N','No','F','False'])

    jd = dict()
    for i in range(len(gl)):
        genome, species, source, run = \
            gl['genome'][i], gl['species'][i], gl['source'][i], gl['run'][i]
        if not run in ['C','T']: continue
        jd1 = dict()
        pre = "%s/%s" % (args.dirg, genome)
        jd1['fasta'] = "%s/10.fasta" % pre
        jd1['fasta_idx'] = "%s/10.fasta.fai" % pre
        jd1['genome_bed'] = "%s/15_intervals/01.chrom.bed" % pre
        jd1['genome_sizes'] = "%s/15_intervals/01.chrom.sizes" % pre
        jd1['macs_gsize'] = get_gsize(jd1['genome_sizes'])
        # annotation
        jd1['gff'] = "%s/50_annotation/10.gff" % pre
        jd1['gtf'] = "%s/50_annotation/10.gtf" % pre
        jd1['bed12'] = "%s/50_annotation/10.bed" % pre
        jd1['transcript_fasta'] = "%s/50_annotation/10.fna" % pre
        jd1['tss_bed'] = "%s/50_annotation/10.tss.bed" % pre
        if gl['star'][i]:
            jd1['star'] = "%s/21_dbs/star/" % pre
        if gl['hisat2'][i]:
            jd1['hisat2'] = "%s/21_dbs/hisat2/db" % pre
            if genome in ['Zmays_B73']:
                jd1['hisat2'] = "%s/21_dbs/hisat2/B73_vt01/db" % pre
        if gl['bwa'][i]:
            jd1['bwa'] = "%s/21_dbs/bwa/db" % pre
        if gl['bismark'][i]:
            jd1['bismark'] = "%s/21_dbs/bismark" % pre
        if gl['salmon'][i]:
            jd1['salmon'] = "%s/21_dbs/salmon/db" % pre
            jd1['tx2gene'] = "%s/21_dbs/salmon/tx2gene.csv" % pre
        if gl['rcfg'][i]:
            jd1['rcfg'] = "%s/55.rds" % pre
        win11 = "%s/15_intervals/20.win11.tsv" % pre
        win56 = "%s/15_intervals/20.win56.tsv" % pre
        win127 = "%s/15_intervals/20.win127.tsv" % pre
        if op.isfile(win11): jd1['win11'] = win11
        if op.isfile(win56): jd1['win56'] = win56
        if op.isfile(win127): jd1['win127'] = win127
        jd1['fc_group_features'] = 'gene_id'
        jd1['fc_group_features_type'] = 'gene_biotype'
        jd[genome] = jd1

    #j = dict(params = dict(genomes = jd))
    j = dict(genomes = jd)
    # with open(args.json, 'w') as outfile:
        # json.dump(j, outfile)
    with open(args.yaml, 'w') as outfile:
        yaml.dump(j, outfile)

def main(args):
    os.system("readGs.py --sheet genomes > %s" % args.fo)
    dump(args)

if __name__ == "__main__":
    import argparse
    ps = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = "prepare genome config file for nextflow pipelines"
    )

    ps.add_argument('--fo', default='genomes.tsv', help = 'output tsv file')
#    ps.add_argument('--json', default='genomes.json', help = 'output json file')
    ps.add_argument('--yaml', default='genomes.yml', help = 'output yaml file')
    ps.add_argument('--dirg', default='/home/springer/zhoux379/projects/genome/data', help='genome directory')

    args = ps.parse_args()
    main(args)

