#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import pandas as pd
import logging
#import simplejson as json
import yaml

from jcvi.apps.base import sh, mkdir

def get_gsize(fs):
    cl = pd.read_csv(fs, sep="\t", header=None, names=['chrom','size'])
    return sum(cl['size'])

def tsv2yml(args):
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
        pre = "%s/%s" % (args.dirg, genome)
        jd1['alias'] = gl['alias'][i]
        jd1['prefix'] = gl['prefix'][i]
        jd1['fasta'] = "%s/10.fasta" % pre
        jd1['fasta_idx'] = "%s/10.fasta.fai" % pre
        jd1['genome_bed'] = "%s/15_intervals/01.chrom.bed" % pre
        jd1['genome_sizes'] = "%s/15_intervals/01.chrom.sizes" % pre
        jd1['gap_bed'] = "%s/15_intervals/11.gap.bed" % pre
        if op.isfile(jd1['genome_sizes']):
            jd1['macs_gsize'] = get_gsize(jd1['genome_sizes'])
        # annotation
        jd1['gff'] = "%s/50_annotation/10.gff" % pre
        jd1['gff_db'] = "%s/50_annotation/10.gff.db" % pre
        jd1['gtf'] = "%s/50_annotation/10.gtf" % pre
        jd1['bed'] = "%s/50_annotation/10.bed" % pre
        jd1['fna'] = "%s/50_annotation/10.nt.fasta" % pre
        jd1['faa'] = "%s/50_annotation/10.aa.fasta" % pre
        jd1['tss'] = "%s/50_annotation/10.tss.bed" % pre
        jd1['pgff'] = "%s/50_annotation/15.gff" % pre
        jd1['pgff_db'] = "%s/50_annotation/15.gff.db" % pre
        jd1['pgtf'] = "%s/50_annotation/15.gtf" % pre
        jd1['pbed'] = "%s/50_annotation/15.bed" % pre
        jd1['pfna'] = "%s/50_annotation/15.nt.fasta" % pre
        jd1['pfaa'] = "%s/50_annotation/15.aa.fasta" % pre
        if gl['blat'][i]:
            jd1['blat'] = "%s/21_dbs/blat/db.2bit" % pre
        if gl['gatk'][i]:
            jd1['gatk'] = f"{pre}/21_dbs/gatk/"
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
        if gl['blast'][i]:
            jd1['blastp'] = f"{pre}/21_dbs/blastp"
            jd1['blastn'] = f"{pre}/21_dbs/blastn"
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

def download(args):
    cvts = dict(genome=str,species=str,source=str)
    gl = pd.read_csv(args.cfg, sep="\t", header=0, converters=cvts,
                     true_values=['1','Y','Yes','T','True'],
                     false_values=['0','N','No','F','False'])

    url_pre = "http://ftp.ebi.ac.uk/ensemblgenomes/pub"
    for i in range(len(gl)):
        if pd.isna(gl['status'][i]) or not gl['status'][i]:
            logging.warning(f"{gl['genome'][i]}: skipped")
            continue

        genome,species,source,version,assembly,url_fas,url_gff = \
            gl['genome'][i], gl['species'][i], gl['source'][i], gl['version'][i], \
            gl['assembly'][i], gl['url_fas'][i], gl['url_gff'][i]

        dirw = f"{args.dirg}/{genome}/raw"
        if not op.isdir(dirw):
            mkdir(dirw)
        os.chdir(dirw)

        if source == 'ensembl':
            version = int(version)
            species = species.replace(" ", "_")
            assembly = assembly.replace(" ", "_")
            url_fas = f"{url_pre}/release-{version}/plants/fasta/{species.lower()}/dna/{species}.{assembly}.dna.toplevel.fa.gz"
            url_gff = f"{url_pre}/release-{version}/plants/gff3/{species.lower()}/{species}.{assembly}.{version}.gff3.gz"
        fn1, fn2 = op.basename(url_fas), op.basename(url_gff)

        comp1, fn1c = False, fn1
        comp2, fn2c = False, fn2
        if not fn1.endswith(".gz"):
            comp1 = True
            fn1c = f"{fn1}.gz"
        if not fn2.endswith(".gz"):
            comp2 = True
            fn2c = f"{fn2}.gz"

        if op.isfile(fn1c) and os.stat(fn1c).st_size > 0 and op.isfile(fn2c) and os.stat(fn2c).st_size > 0:
            logging.warning(f"{genome}: already done")
            continue
        else:
            logging.debug(f"{genome}: working")

        if source == 'local':
            if not op.isfile(url_fas):
                url_fas = op.join(args.dirg, url_fas)
            if not op.isfile(url_gff):
                url_gff = op.join(args.dirg, url_gff)
            if not op.isfile(url_fas) or not op.isfile(url_gff):
                logging.error(f"no fasta/gff found: {url_fas} {url_gff}")
            sh(f"cp {url_fas} {fn1}")
            sh(f"cp {url_gff} {fn2}")
        else:
            sh(f"axel -n {args.thread} {url_fas} -o {fn1}")
            sh(f"axel -n {args.thread} {url_gff} -o {fn2}")

        if comp1:
            sh(f"gzip {fn1}")
        if comp2:
            sh(f"gzip {fn2}")
        sh(f"ln -sf {fn1c} raw.fasta.gz")
        sh(f"ln -sf {fn2c} raw.gff.gz")


if __name__ == "__main__":
    import argparse
    ps = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = "help utilities to run the nextflow genome pipeline"
    )
    sp = ps.add_subparsers(title = 'available commands', dest = 'command')

    sp1 = sp.add_parser("conv",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            help = "prepare genome config file for nextflow pipelines")
    sp1.add_argument('fi', help = 'input genome tsv')
    sp1.add_argument('fo', help = 'output yaml file')
    sp1.add_argument('--dirg', default="%s/zhoup-genome" % os.environ["s3"], help='genome directory')
#    sp1.add_argument('--json', default='genomes.json', help = 'output json file')
    sp1.set_defaults(func = tsv2yml)

    sp1 = sp.add_parser("download",
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            help = "download raw genome (fasta+gff) files")
    sp1.add_argument('--cfg', default=f"{os.environ['genome']}/nf/genomes.tsv", help='genome config table')
    sp1.add_argument('--dirg', default="%s/zhoup-genome" % os.environ["s3"], help='genome directory')
    sp1.add_argument('--thread', type=int, default=4, help='downloading threads')
    sp1.set_defaults(func = download)

    args = ps.parse_args()
    if args.command:
        args.func(args)
    else:
        print('Error: need to specify a sub command\n')
        ps.print_help()

