#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import playr as M
import seqprimers as sp
from chemdraw import workflow_panel
from page import head,info,table
OUT=H.parent/"playr.html"
def render():
  oligos=[("Probe 1 reporter tail",f"5′-{M.PROBE1_TAIL}-3′"),("Probe 2 reporter tail",f"5′-{M.PROBE2_TAIL}-3′"),("Insert 1",f"5′-p-{M.INSERT}-3′"),("Universal backbone",f"5′-p-{M.BACKBONE}-3′"),("Detection sequence",f"5′-{M.DETECTION}-3′")]
  return "\n".join([head(M.TITLE),'<div class="wrap" data-assay-complete="true">',f'<h1>{M.TITLE}</h1>','<p class="research-notes"><a href="01_playr.html">Research notes</a></p>',info(M.SOURCE),
    info('System 1 from the published 27-system table is shown. Gene-specific target arms vary with the chosen transcript and are therefore symbolic.'),
    '<h2>Probe-pair topology</h2>',M.TOPOLOGY.panel('Probe 2 binds upstream and Probe 1 downstream. Their 5′ target arms anneal antiparallel to RNA; the 3′ tails turn upward. The inner insert and outer backbone each bridge both probes, and the two dots are the T4-ligase junctions.'),
    '<h2>Nucleotide-resolved fold</h2>',M.TOPOLOGY.folded_panel('The same ligation-ready complex with every modelled nucleotide shown. N denotes the unpublished gene-specific target sequence. The target RNA is constrained straight; paired spans determine the helices, while unpaired spans are routed as flexible joints. The renderer rejects crossings and incompatible base pairs.'),
    '<h2>Reaction workflow</h2>',workflow_panel(M.workflow(),cls="long"),
    '<h2>Published system-1 reporter oligos</h2>',table(("Oligo","Sequence as modelled"),oligos),
    sp.read_spec_table([sp.read_spec("Sequencing output",None,"No sequencing library or nucleotide reads","not_applicable")],"PLAYR is quantified by flow or mass cytometry."),'</div>'])
if __name__=="__main__": OUT.write_text(render(),encoding="utf-8"); print(f"wrote {OUT}")
