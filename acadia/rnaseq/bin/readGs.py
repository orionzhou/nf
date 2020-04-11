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
    df2 = df[df.libtype.eq(args.libtype) & df.Run.eq('T')]
    df3 = df2.filter(items=['yid','libtype','author','year','ASE','stress','RIL','Run'])
    print(df3.to_csv(index=False, sep="\t"), end='')

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'read google spreadsheet'
    )
    p.add_argument('--libtype', default='rnaseq', help='library type')
    p.add_argument('--book', default='maize_studies', help = 'name of google spreadsheet')
    p.add_argument('--sheet', default='all', help = 'sheet/tab name')
    p.add_argument('--cred', default='/home/springer/zhoux379/.config/google_account_token.json', help = 'google credential')

    args = p.parse_args()
    if not op.isfile(args.cred):
        p.print_help()
    else:
        gs2tsv(args)
