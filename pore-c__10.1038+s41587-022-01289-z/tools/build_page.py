#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/'lib')]
import pore_c as P
import seqprimers as sp
from chemdraw import panel,protocol_graph_panels,workflow_panel
from page import head,info
OUT=H.parent/'pore-c.html'
def render():
 graph=P.protocol_graph()
 return '\n'.join([head('Pore-C library chemistry'),'<div class="wrap" data-protocol-complete="true">','<h1>Pore-C &mdash; intact multiway proximity concatemers</h1>','<p class="research-notes"><a href="01_pore-c.html">Research notes</a></p>',info('Defining source: <a href="https://doi.org/10.1038/s41587-022-01289-z">Deshpande et al., <i>Nature Biotechnology</i> (2022)</a>. The extraction steps follow Oxford Nanopore\'s later restriction-enzyme Pore-C procedure. The optional size-selection and dA/dT ligation panels depict its later native-library branch; the original paper\'s kit revision remains unverified.'),info('Restriction fragments that were close in fixed nuclei ligate into one long chimeric polymer. Nanopore sequencing keeps three or more contact partners on the same physical read.'),info('SOURCE DISCREPANCY — The ONT Pore-C info sheet prints document version RPC_S1015_v1_revL_12Nov2019 in its overview, yet its library-preparation section recommends SQK-LSK114 and calls SQK-LSK109/110 older alternatives. The model uses the later kit guidance only for a labelled example; the defining paper\'s kit and adapter remain unresolved.'),'<h2>Reaction workflow</h2>',protocol_graph_panels(graph,cls='long'),'<h2>Molecular details</h2>',workflow_panel(P.workflow(),cls='long'),'<h2>Final nanopore library</h2>',panel(P.final_rows(),cls='long',caption='A three-fragment contact is drawn. Real molecules contain a variable number and order of restriction fragments; ** marks dA/dT ligation geometry established for a later ONT kit, not yet for the 2022 paper.'),'<h2>Nanopore entry</h2>',info('There is no sequencing primer. A motor-bearing ONT ligation adapter feeds the intact proximity concatemer through the pore; software later splits its alignment at restriction-fragment boundaries while retaining one read identity.'),sp.read_spec_table([sp.read_spec('Sequencing output',None,'intact multiway proximity concatemer; alignments split at restriction-fragment boundaries','include')],'Nanopore sequencing is primer-free; read length and fragment count are molecule-dependent.'),'</div>'])
if __name__=='__main__': OUT.write_text(render(),encoding='utf-8'); print(f'wrote {OUT}')
