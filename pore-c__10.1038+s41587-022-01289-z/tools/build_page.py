#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/'lib')]
import pore_c as P
from chemdraw import panel,workflow_panel
from page import head,info
OUT=H.parent/'pore-c.html'
def render():
 return '\n'.join([head('Pore-C library chemistry'),'<div class="wrap">','<h1>Pore-C &mdash; intact multiway proximity concatemers</h1>','<p class="research-notes"><a href="01_pore-c.html">Research notes</a></p>',info('Defining source: <a href="https://doi.org/10.1038/s41587-022-01289-z">Deshpande et al., <i>Nature Biotechnology</i> (2022)</a>; molecular procedure cross-checked against Oxford Nanopore\'s restriction-enzyme Pore-C protocol.'),info('Restriction fragments that were close in fixed nuclei ligate into one long chimeric polymer. Nanopore sequencing keeps three or more contact partners on the same physical read.'),'<h2>Reaction workflow</h2>',workflow_panel(P.workflow(),cls='long'),'<h2>Final nanopore library</h2>',panel(P.final_rows(),cls='long',caption='A three-fragment contact is drawn. Real molecules contain a variable number and order of restriction fragments; ** marks checked T:A adapter ligations.'),'<h2>Nanopore entry</h2>',info('There is no sequencing primer. A motor-bearing ONT ligation adapter feeds the intact proximity concatemer through the pore; software later splits its alignment at restriction-fragment boundaries while retaining one read identity.'),'</div>'])
if __name__=='__main__': OUT.write_text(render(),encoding='utf-8'); print(f'wrote {OUT}')
