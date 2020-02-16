#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import os.path as op
import pandas as pd
import gspread
from oauth2client.service_account import ServiceAccountCredentials

def read_gs(book, sheet, f_cred):
    scope = ['https://spreadsheets.google.com/feeds',
             'https://www.googleapis.com/auth/drive']
    cred = ServiceAccountCredentials.from_json_keyfile_name(f_cred, scope)
    gc = gspread.authorize(cred)

    #bk = gc.open_by_key('1SacBnsUW4fzqGYl0k5FVV5AFq2TJvUlBOIWp1NzLh88')
    print(gc.list_spreadsheet_files())
    bk = gc.open(book)
    st = bk.worksheet(sheet)
    df = pd.DataFrame(st.get_all_records())
    df.to_csv(index=False, sep="\t")

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = 'read google spreadsheet'
    )
    p.add_argument('book', help = 'name of google spreadsheet')
    p.add_argument('sheet', help = 'worksheet/tab name')
    p.add_argument('--cred', default='/home/springer/zhoux379/.config/google_service_account.json', help = 'google credential')

    args = p.parse_args()
    if not op.isfile(args.cred):
        p.print_help()
    else:
        read_gs(args.book, args.sheet, args.cred)
