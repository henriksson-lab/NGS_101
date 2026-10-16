#!/usr/bin/env python3
"""Check the corrected original SI, rather than generic TruSeq defaults."""
from pathlib import Path
import re
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
from storage import data_root
import duplex_seq as model
from checks import Check, Source, have

check=Check()
data=data_root()
source=Source(data/'sources'/HERE.parent.name/'source_pdf_837ee6137323e8b633b0.pdf.txt',
    'Original Schmitt2012 article + Supporting Information corrected March27,2013; registered public USPTO exhibit PDF route in catalogue/extra_sources.tsv')
if have(check, source, label='Corrected original Supporting Information'):
    text=source.path.read_text()
    si=text[text.index('Supporting Information Corrected March 27, 2013'):]
    flat=re.sub(r'\s+', '', si)
    check('adapter primer strand matches original corrected SI', model.PRIMER_STRAND in flat)
    check('phosphorylated template strand matches original corrected SI',
          '/5phos/'+model.TEMPLATE_STRAND in flat)
    check('both PCR primers preserve source-specific arms and six-position index',
          model.PCR_FORWARD in flat and model.PCR_REVERSE.replace('N'*6,'X'*6) in flat)
    tag_length=int(re.search(r'a (\d+)-\s*nucleotide random sequence', si)[1])
    spacer_length=int(re.search(r'a (\d+)-nucleotide ﬁxed', si)[1])
    source_spacer=re.search(r'ﬁxed sequence ([ACGT]+) present after exactly', si)[1]
    check('copied ACTG and added A reproduce the verbatim source spacer', model.READ_SPACER, source_spacer)
    for primer in (model.sp.TRUSEQ['R1'],model.sp.TRUSEQ['R2']):
        hit=model.sp.locate(model.FINAL_LIBRARY, primer, n=tag_length+spacer_length)
        check(primer.role+' begins with the source-length tag and fixed spacer',
              len(hit.reads[:tag_length+spacer_length]) == tag_length+spacer_length
              and hit.reads[:tag_length].upper() == 'N'*tag_length
              and hit.reads[tag_length:tag_length+spacer_length] == source_spacer)
check.report()
