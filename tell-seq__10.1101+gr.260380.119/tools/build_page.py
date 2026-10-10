#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import tell_seq as T
import seqprimers as sp
from chemdraw import Scene, annotation_rows, panel, workflow_panel
from page import head, info
OUT=H.parent/'tell-seq.html'

def render():
    return '\n'.join([head('TELL-seq linked-read chemistry'),'<div class="wrap">',
      '<h1>TELL-seq &mdash; single-tube bead-linked reads</h1>',
      '<p class="research-notes"><a href="01_tell-seq.html">Research notes</a></p>',
      info('Defining source: <a href="https://doi.org/10.1101/gr.260380.119">Chen et al., <i>Genome Research</i> (2020)</a>.'),
      info('Stable MuA strand-transfer complexes and hybrid capture keep fragments from one long DNA molecule associated with one clonally barcoded bead in an unpartitioned tube.'),
      '<h2>Reaction workflow</h2>',workflow_panel(T.workflow(),cls='long'),
      '<h2>Final Illumina library</h2>',panel(T.rows(T.FINAL,'TELL-seq'),cls='long',caption='Index 1 is the 18-base linked-molecule barcode; Index 2 is the 8-base sample index. Adaptor ordering is shown from the published architecture; proprietary bases remain role tokens.'),
      '<h2>Sequencing-primer binding</h2>',
      info('TELL-seq requires vendor custom Read 1, Read 2 and Index 1 primers, plus a custom Index 2 primer on instruments that do not use grafted P5 for i5. Their sequences are not disclosed in the paper or public sequencing guide, so exact annealing cannot be computed without fabricating bases. The guide specifies 2 × 146 reads, an 18-cycle Index 1 read and an 8-cycle Index 2 read.'),
      sp.unavailable_diagram(T.FINAL,'custom primer sequence is proprietary',roles=('Read 1','Read 2','Index 1 (18 cycles)','Index 2 (8 cycles)')),'</div>'])
if __name__=='__main__': OUT.write_text(render(),encoding='utf-8'); print(f'wrote {OUT}')
