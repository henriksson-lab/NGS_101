#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import qbc as M
import seqprimers as sp
from chemdraw import duplex_rows, panel, protocol_graph_panels
from page import head,info,table
OUT=H.parent/"qbc.html"
def render():
  rows=[("1: SC1","7-nt cell subcode"),("2: SCN","7-nt cell subcode"),("3: SC2","7-nt cell subcode"),("4: SC3","10-nt UMI + 7-nt cell subcode")]
  return "\n".join([head("Quantum Barcoding chemistry"),'<div class="wrap" data-protocol-complete="true">','<h1>Quantum Barcoding (QBC)</h1>',
    '<p class="research-notes"><a href="01_qbc.html">Research notes</a></p>',info('Defining source: <a href="https://doi.org/10.1038/s42003-020-0896-2">Roy et al., <i>Communications Biology</i> (2020)</a>.'),
    '<h2>Splint-directed split-pool indexing</h2>',panel(M.splint_scene().rows(),cls="long"),
    '<h3>Nucleotide-resolved Splint6 assembly</h3>',M.SPLINT6_TOPOLOGY.folded_panel(
      'Splint6 is one DNA scaffold interrupted by four non-nucleotide C3 spacers. Exact '
      'published docking arms place SC1, SCN, SC2 and SC3 antiparallel; N marks only the '
      'published variable UMI and cell-subcode positions. Dots mark the four ligation '
      'junctions that produce the continuous reporter shown above the splint.'),
    table(("Round","Molecular content"),rows),
    '<h2>Reaction workflow</h2>',protocol_graph_panels(M.protocol_graph(),cls="long",annotations=True),
    '<h2>Antibody reporter library</h2>',panel(duplex_rows(M.PROTEIN_LIBRARY),cls="long"),sp.section(M.PROTEIN_LIBRARY,M.PRIMERS),
    '<h2>Targeted-RNA library</h2>',panel(duplex_rows(M.RNA_LIBRARY),cls="long"),sp.section(M.RNA_LIBRARY,M.PRIMERS),
    info('Index structure: the Methods name NEBNext E7600S (targeted-RNA PCR) and E7335 (protein-only PCR) primer sets. The vendor index layouts were not inspected in this pass, so both libraries are drawn with placeholder dual 8-nt indices; the acquired index cycles are not stated in the article.'),
    info('The source reports 151–258-base single- or paired-end runs; paired reads were overlapped with FLASH. Read 1 traverses UMI, four cell subcodes, then antibody code or target cDNA; paired Read 2 approaches from the opposite adapter.'),'</div>'])
if __name__=="__main__": OUT.write_text(render(),encoding="utf-8"); print(f"wrote {OUT}")
