#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import os.path as op
import pandas as pd
import gspread
from oauth2client.service_account import ServiceAccountCredentials

def gs2tsv(args):
    scope = ['https://spreadsheets.google.com/feeds',
             'https://www.googleapis.com/auth/drive']
    cred = ServiceAccountCredentials.from_json_keyfile_name(args.cred, scope)
    gc = gspread.authorize(cred)

    #bk = gc.open_by_key('1SacBnsUW4fzqGYl0k5FVV5AFq2TJvUlBOIWp1NzLh88')
    #print(gc.list_spreadsheet_files())
    bk = gc.open(args.book)
    st = bk.worksheet(args.sheet)
    df = pd.DataFrame(st.get_all_records())
    if args.workflow != '':
        df = df[df.workflow.eq(args.workflow) & df.Run.eq('C')]
        if args.workflow == 'rnaseq':
            df = df.filter(items=['yid','author','year','source','accession','study','genotype','tissue','n','ASE','stress','RIL','Run'])
        else:
            df = df.filter(items=['yid','author','year','source','accession','study','genotype','tissue','n','stress','Run'])
    print(df.to_csv(index=False, sep="\t"), end='')

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'read google spreadsheet'
    )
    p.add_argument('--book', default='maize_studies', help = 'name of google spreadsheet')
    p.add_argument('--sheet', default='all', help = 'sheet/tab name')
    p.add_argument('--cred', default='$HOME/.config/google_account_token.json', help = 'google credential')
    p.add_argument('--workflow', default='', help='filter by workflow')

    args = p.parse_args()
    if not op.isfile(args.cred):
        p.print_help()
    else:
        gs2tsv(args)
