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

def main(args):
    cvts = dict(genome=str,species=str,source=str)
    gl = pd.read_csv(args.fi, sep="\t", header=0, converters=cvts,
                     true_values=['1','Y','Yes','T','True'],
                     false_values=['0','N','No','F','False'])

    jd = dict()
    for i in range(len(gl)):
        genome, species, source, status = \
            gl['genome'][i], gl['species'][i], gl['source'][i], gl['status'][i]
        #if not status in ['C','T']: continue
        jd1 = dict()
        pre = "%s/data/%s" % (args.dirg, genome)
        jd1['fasta'] = "%s/10.fasta" % pre
        jd1['fasta_idx'] = "%s/10.fasta.fai" % pre
        jd1['genome_bed'] = "%s/15_intervals/01.chrom.bed" % pre
        jd1['genome_sizes'] = "%s/15_intervals/01.chrom.sizes" % pre
        if op.isfile(jd1['genome_sizes']):
            jd1['macs_gsize'] = get_gsize(jd1['genome_sizes'])
        # annotation
        jd1['gff'] = "%s/50_annotation/10.gff" % pre
        jd1['gtf'] = "%s/50_annotation/10.gtf" % pre
        jd1['bed'] = "%s/50_annotation/10.bed" % pre
        jd1['fna'] = "%s/50_annotation/10.nt.fasta" % pre
        jd1['faa'] = "%s/50_annotation/10.aa.fasta" % pre
        jd1['tss'] = "%s/50_annotation/10.tss.bed" % pre
        jd1['pgff'] = "%s/50_annotation/15.gff" % pre
        jd1['pgtf'] = "%s/50_annotation/15.gtf" % pre
        jd1['pbed'] = "%s/50_annotation/15.bed" % pre
        jd1['pfna'] = "%s/50_annotation/15.nt.fasta" % pre
        jd1['pfaa'] = "%s/50_annotation/15.aa.fasta" % pre
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
    with open(args.fo, 'w') as outfile:
        yaml.dump(j, outfile)
    # with open(args.json, 'w') as outfile:
        # json.dump(j, outfile)

if __name__ == "__main__":
    import argparse
    ps = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = "prepare genome config file for nextflow pipelines"
    )

    ps.add_argument('fi', help = 'input genome tsv')
    ps.add_argument('fo', help = 'output yaml file')
    ps.add_argument('--dirg', default=os.environ["genome"], help='genome directory')
#    ps.add_argument('--json', default='genomes.json', help = 'output json file')

    args = ps.parse_args()
    main(args)

