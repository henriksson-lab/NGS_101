#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
import moffat_readout as P
# Preserve the existing self-test entry points while this migration awaits review.
from moffat_readout import pcr1_product, D501_F, D701_R
import batch_page
OUT=HERE.parent / 'moffat-two-step-crispr-screen-readout.html'
if __name__=='__main__':
 OUT.write_text(batch_page.render(P),encoding='utf-8')
 print('wrote',OUT)
