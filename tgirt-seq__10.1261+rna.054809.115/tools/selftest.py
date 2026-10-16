#!/usr/bin/env python3
"""Original S1 transcription checks; molecular invariants live in constructors."""
from pathlib import Path
import re
import sys
H = Path(__file__).resolve().parent
sys.path[:0] = [str(H), str(H.parents[1] / 'lib')]
from storage import data_root
import tgirt_seq as M
from checks import Check, Source, have

check = Check()
bundle = data_root() / 'sources' / H.parent.name
source = Source(bundle / 'PMC4691826_supp_supp_054809.115_FigS1.eps.txt',
                'Original Qin2016 Fig.S1 EPS via pmccloud:PMC4691826; ref/MANIFEST.md describes text extraction.')
if have(check, source, label='original Supplementary Figure S1 text'):
    text = re.sub(r'\s+', '', source.path.read_text())
    # S1 typesets the C3-spacer name /3'SpC3/ in the RNA row. Methods calls
    # the same block 3SpC3. Normalize this one named typography difference.
    text = text.replace("/3'SpC3/", '/3SpC3/').replace('/3′SpC3/', '/3SpC3/').replace('/3’SpC3/', '/3SpC3/')
    text = text.replace('PCRprimer', '')  # left-column label splits the wrapped P7 row
    check('complete original RNA starter and C3-spacer designation', M.R2_RNA_WRITTEN in text)
    check('complete original R2R primer including mixed overhang', M.R2R_DNA in text)
    check('complete ordered R1R including 5-phosphate and 3-block', M.R1R_WRITTEN in text)
    check('complete original unindexed P5 primer', M.P5_PRIMER in text)
    reverse = M.P7_PARTS[0].top + 'BARCODE' + M.P7_PARTS[-1].top
    check('complete original indexed P7 primer with source placeholder', reverse in text)
check.report()
