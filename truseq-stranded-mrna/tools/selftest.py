#!/usr/bin/env python3
"""Source transcription and local workflow boundaries."""
from pathlib import Path
import os
import re
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / 'lib')]
import chem as m
from checks import Check, Source, have
c = Check()
bundle = Path(os.environ.get('CHEM_DATA', HERE.parents[1] / '_data')) / 'sources' / HERE.parent.name
workflow = Source(bundle / 'illumina_truseq-stranded-mrna-workflow-reference-1000000040498-00.pdf.txt', 'Registered v00 workflow original')
adapters = Source(bundle / 'illumina_1000000002694_22_ilmn-adapter-sequences.pdf.txt', 'Registered v22 adapter original')
if have(c, adapters, label='actual source CD templates'):
    text = next(block for block in re.findall(r'TruSeq DNA and RNA CD Indexes(.*?)TruSeq Single Indexes', adapters.path.read_text(), re.S) if 'D701' in block)
    left = m.il.P5 + '[i5]' + m.il.TRUSEQ_READ1
    right = m.il.INDEX1_PRIMER + '[i7]' + m.il.P7_RC
    c('complete i5 adapter equals original CD template', left in text)
    c('complete i7 adapter equals original CD template', right in text)
    for name in ('D701', 'D501'):
        member = re.search(name + r'\s+([ACGT]+)', text)[1]
        index = next(s for s in m.FINAL_LIBRARY if s.name.startswith('i7 CD' if name == 'D701' else 'i5 CD'))
        c(name + ' source family width retained without selecting member', index.length_bp, len(member))
if have(c, workflow, label='original workflow choices'):
    text = workflow.path.read_text()
    c('PCR cycles transcribed from original', m.PCR_CYCLES, int(re.search(r'u (\d+) cycles of:', text)[1]))
    c('SMM explicitly produces blunt-ended cDNA', 'The result is blunt-ended cDNA.' in text)
    c('source explicitly removes RNA during SMM', 'This process removes the RNA template' in text)
    c('source ligation cleanup has two rounds', 'Repeat steps 2 through 17 with the new plate using the Round 2 volumes.' in text)
# These regressions catch the prior shared recipe's extra indexing PCR and omitted recoveries.
g = m.protocol_graph()
order = list(g.steps)
c('both ligation cleanups precede selective PCR', order.index('ligate') < order.index('clean_ligation_1') < order.index('clean_ligation_2') < order.index('pcr') < order.index('clean_pcr'))
c('historical adapter assignment is inferred by default', all(s.inferred for s in m.cd_library() if s.name in ('P5', 'Read 1 adapter arm', 'i7 adapter arm', 'P7 reverse complement')))
c.report()
