#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
import lenticrispr_v2 as P
import batch_page
OUT=HERE.parent / 'lenticrispr-v2-screening.html'
if __name__=='__main__':
 OUT.write_text(batch_page.render(P),encoding='utf-8')
 print('wrote',OUT)
