#!/usr/bin/env python3
"""Original-source checks for i-BLESS declarations and restored reaction order."""
import os
import re
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
from checks import Check, Source, have
import chem
c = Check()
bundle = Path(os.environ.get('CHEM_DATA', HERE.parents[1] / '_data')) / 'sources' / HERE.parent.name
paper = Source(bundle / 'article.xml.txt', 'Original i-BLESS Methods')
supp = Source(bundle / 'supplement.pdf.txt', 'Original i-BLESS Supplementary Table 4')
if have(c, paper, supp, label='Original i-BLESS sources'):
    text = paper.path.read_text()
    table = supp.path.read_text().split('Supplementary Table 4.', 1)[1]
    letters = re.sub(r'[^A-Z]', '', table.replace('[Phos]', '').replace('[BtndT]', 'T'))
    c('all four whole published oligos are molecular graph inputs',
      all(seq in letters for seq in (chem.BLESS_PROXIMAL, chem.BLESS_DISTAL, chem.BLESS_PCR_PROXIMAL, chem.BLESS_PCR_DISTAL)) and
      all(chem.protocol_graph().materials[name].state is not None for name in ('proximal_adapter', 'distal_adapter', 'insert_forward', 'insert_reverse')))
    c('proximal biotin and both terminal phosphate modifications retained',
      table.count('[Phos]') == 2 and '[BtndT]' in table and '/iBiodT/' in chem.oligos()[0])
    steps = list(chem.protocol_graph().steps)
    c('source two sonications remain separated by DNA recovery',
      'once again fragmented by sonication' in text and
      steps.index('initial_sonication') < steps.index('recover_dna') < steps.index('fragment') < steps.index('capture'))
    c('source insert PCR and purification precede XhoI',
      'Purified PCR products were subsequently treated with XhoI' in text and
      steps.index('linearize') < steps.index('insert_pcr') < steps.index('purify_insert') < steps.index('trim_tags'))
    c('source paired 75-base acquisition stays distinct from index configuration',
      '2 × 75' in text.replace('\u2009', ' ') and chem.READ_LENGTHS == {'Read 1': 75, 'Read 2': 75})
manual = Source(bundle / 'ThruPLEX-DNA-seq-Kit-Manual-QAM-108-003-updated-9_16_16.pdf.txt', 'Original September 2016 ThruPLEX guide')
if have(c, manual, label='Contemporaneous kit reference'):
    text = manual.path.read_text()
    c('guide stem-loop mechanism replaces inherited dA/Y-adapter endpoint',
      'stem-loop adaptors with blocked 5’ ends' in text and 'leaving a nick at the 3’ end' in text and
      all(segment.is_role_token() for segment in chem.FINAL_LIBRARY))
    c('guide reagent indexes do not silently define study index cycles',
      'only the first 6 bases of the 8nt sequence will be sequenced' in text and
      'dual indexes are 8nt long' in text and not any('Index' in key for key in chem.READ_LENGTHS))
c.report()
