#!/usr/bin/env python3
"""Build the diagram-centric Seq-Well chemistry page."""
from __future__ import annotations
import html,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE)); sys.path.insert(0,str(HERE.parents[1]/"lib"))
import seq_well as S
import seqprimers as sp
from chemdraw import Segment,oligo,protocol_graph_panels
from page import head,info,table
OUT=HERE.parent/"seq-well.html"


def sg(name,top,tag=None,**kw): return Segment(name=name,top=top,tag=tag,**kw)


def oligos()->str:
    rows=[oligo("Barcoded bead oligo",S.bead_primer(),five="5'-Bead-Linker-"),
          oligo("Template-switching oligo",S.tso(),three="-3' (terminal rGrGrG)"),
          oligo("SMART PCR primer",[sg("SMART handle",S.SMART_PCR,"tso")]),
          oligo("P5-SMART selective PCR primer",S.p5_hybrid(),three="-3' (two terminal phosphorothioates)"),
          oligo("Custom Read 1 primer",[sg("spacer + handle + AC",S.CUSTOM_R1,"r1")])]
    return "<h2>Key oligos</h2>\n<seq>\n"+"\n".join(rows)+"\n</seq>"


def construction()->str:
    return '<h2>Library construction</h2>\n'+protocol_graph_panels(S.protocol_graph(),cls="long")


def sequencing()->str:
    lib=S.final_library()
    rows=[]
    for p,h in S.primer_landings(): rows.append((html.escape(p.role),html.escape(p.name),html.escape(", ".join(h.covers)),f"<code>{html.escape(h.reads)}</code>&hellip;"))
    return f'''<h2>Sequencing</h2>
{info('Read 1: 20 cycles with the custom primer; Index 1: 8; Read 2: 50. Read 1 reports cell barcode + UMI, and Read 2 reports transcript sequence.')}
{sp.diagram(lib,S.SEQ_PRIMERS)}
{table(("Read","Primer","Primer site","First bases"),rows)}
<h3>Read 1 layout</h3>{table(("Cycles","Content"),S.read1_layout(),scroll=False)}
{sp.read_content_table(lib,S.SEQ_PRIMERS,
    {"Read 1":20,"Index 1 (i7)":8,"Read 2":50},
    acquired_roles=tuple(p.role for p in S.SEQ_PRIMERS))}
'''


def main()->None:
    page="\n".join([head("Seq-Well library chemistry"),'<div class="wrap" data-protocol-complete="true">','<h1>Seq-Well &mdash; barcoded-bead RNA capture in sealed microwells</h1>',info('Defining source: <a href="https://doi.org/10.1038/nmeth.4179">Gierahn et al., <i>Nature Methods</i> (2017)</a>.'),oligos(),construction(),sequencing(),'</div>'])
    OUT.write_text(page,encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__=="__main__": main()
