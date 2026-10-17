#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import ncats as N
import seqprimers as sp
from chemdraw import panel, protocol_graph_panels, workflow_panel
from page import head, info
from protocol_graph import AnnotationSlot
from protocol_views import render_annotations
OUT=H.parent/'ncats.html'

def render():
    graph=N.protocol_graph()
    return '\n'.join([head('nCATS targeted nanopore chemistry'),'<div class="wrap" data-protocol-complete="true">',
      '<h1>nCATS &mdash; Cas9-guided nanopore adapter ligation</h1>',
      '<p class="research-notes"><a href="01_ncats.html">Research notes</a></p>',
      info('Defining source: <a href="https://doi.org/10.1038/s41587-020-0407-5">Gilpatrick et al., <i>Nature Biotechnology</i> (2020)</a>.'),
      info('Dephosphorylate every pre-existing DNA end, then let Cas9 expose fresh phosphorylated ends at chosen loci. Only those new ends efficiently receive nanopore adapters.'),
      render_annotations(graph,AnnotationSlot.BEFORE_WORKFLOW),'<h2>Reaction workflow</h2>',protocol_graph_panels(graph,cls='long'),
      '<h2>Molecular details</h2>',workflow_panel(N.workflow(),cls='long'),
      '<h2>Final targeted library</h2>',panel(N.final_rows(),cls='long',caption='The selected native target remains unamplified, preserving long-range sequence and base modifications. ** marks the exact checked T:A adapter-ligation boundaries; adapter roles are proprietary.'),
      '<h2>Nanopore entry</h2>',info('There is no sequencing primer. An LSK109 motor-loaded adapter engages a pore and controls passage of the intact native target strand. Reads begin at a Cas9-created adapter-bearing end and extend into the region of interest.'),
      sp.read_spec_table([sp.read_spec('Sequencing output',None,
          'native target strand from a Cas9-created adapter-bearing end','include')],
          'Nanopore sequencing is primer-free; read length is molecule-dependent.'),'</div>'])
if __name__=='__main__': OUT.write_text(render(),encoding='utf-8'); print(f'wrote {OUT}')
