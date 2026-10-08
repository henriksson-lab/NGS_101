#!/usr/bin/env python3
"""Build the diagram-centric Seq-Well chemistry page."""
from __future__ import annotations
import html,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE)); sys.path.insert(0,str(HERE.parents[1]/"lib"))
import seq_well as S
import seqprimers as sp
from chemdraw import Scene,Segment,annotation_rows,oligo,panel,strand_row
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
    first=S.first_strand(); wta=S.wta_product(); tagged=S.selected_tagmented_strand(); lib=S.final_library()
    return f'''<h2>Library construction</h2>
<h3>(1) Confine one cell and one bead in a membrane-sealed microwell</h3>
{panel(S.capture_scene().rows(),cls="small",caption="After lysis, poly(A) RNA hybridizes to that well's bead. Reverse transcription is performed after beads are recovered and pooled.")}
{panel([strand_row(first),*annotation_rows(first)],cls="small",caption="One bead-bound first strand after RT adds the untemplated CCC switch site.")}

<h3>(2) Template-switch, remove unused bead primers, and amplify full-length cDNA</h3>
{panel(S.template_switch_scene().rows(),cls="small",caption="The TSO's rGrGrG pairs with CCC and contributes a second SMART handle.")}
{panel([*Scene.duplex(list(wta),label="WTA").rows(),*annotation_rows(wta)],cls="small",caption="Template-switched cDNA before SMART-handle PCR.")}

<h3>(3) Tagment amplified cDNA; select the bead end with P5-SMART + N7xx PCR</h3>
{panel([strand_row(tagged),*annotation_rows(tagged)],cls="small",caption="Only the bead end carries handle+AC for the selective P5 primer; the retained fragment ends at the nearest s7 insertion.")}
{panel([*Scene.duplex(list(lib),label="library").rows(),*annotation_rows(lib)],cls="small",caption="Final Seq-Well library, P5 to P7'. Samples are indexed on i7 only.")}
'''


def sequencing()->str:
    lib=S.final_library()
    rows=[]
    for p,h in S.primer_landings(): rows.append((html.escape(p.role),html.escape(p.name),html.escape(", ".join(h.covers)),f"<code>{html.escape(h.reads)}</code>&hellip;"))
    return f'''<h2>Sequencing</h2>
{info('Read 1: 20 cycles with the custom primer; Index 1: 8; Read 2: 50. Read 1 reports cell barcode + UMI, and Read 2 reports transcript sequence.')}
{sp.diagram(lib,S.SEQ_PRIMERS)}
{table(("Read","Primer","Primer site","First bases"),rows)}
<h3>Read 1 layout</h3>{table(("Cycles","Content"),S.read1_layout(),scroll=False)}
'''


def main()->None:
    page="\n".join([head("Seq-Well library chemistry"),'<div class="wrap">','<h1>Seq-Well &mdash; barcoded-bead RNA capture in sealed microwells</h1>',info('Defining source: <a href="https://doi.org/10.1038/nmeth.4179">Gierahn et al., <i>Nature Methods</i> (2017)</a>.'),oligos(),construction(),sequencing(),'</div>'])
    OUT.write_text(page,encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__=="__main__": main()
