#!/usr/bin/env python3
"""Build the diagram-centric Drop-seq chemistry page."""
from __future__ import annotations
import html,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE)); sys.path.insert(0,str(HERE.parents[1]/"lib"))
import drop_seq as D
from chemdraw import Scene,Segment,annotation_rows,oligo,panel,strand_row
from page import caveat,head,info,table
OUT=HERE.parent/"drop-seq.html"


def sg(name,top,tag=None,**kw): return Segment(name=name,top=top,tag=tag,**kw)


def oligos()->str:
    rows=[oligo("Barcoded bead oligo, batch A",D.bead_primer("A"),five="5'-Bead-Linker-"),
          oligo("Barcoded bead oligo, batch B",D.bead_primer("B"),five="5'-Bead-Linker-"),
          oligo("Template-switching oligo",D.tso(),three="-3' (terminal rGrGrG)"),
          oligo("P5-SMART selective PCR primer, batch B",D.p5_hybrid("B")),
          oligo("Custom Read 1 primer, batch B",[D.inferred("spacer + handle + AC",D.custom_r1("B"),"r1")])]
    return "<h2>Key oligos</h2>\n"+info("INFERRED &mdash; exact bases are reconstructed from the secondary schematic and the later Seq-Well oligo table; Drop-seq Table S6 was unavailable.")+"<seq>\n"+"\n".join(rows)+"\n</seq>"


def construction()->str:
    first=D.first_strand(); selected=D.selected_tagmented_strand(); lib=D.final_library()
    return f'''<h2>Library construction</h2>
<h3>(1) Co-encapsulate one cell and one bead in a nanoliter droplet</h3>
{panel(D.capture_scene().rows(),cls="small",caption="INFERRED — reconstructed batch-B bead sequence. Poly(A) RNA hybridizes inside the droplet; RT occurs after pooled beads are recovered.")}
{panel([strand_row(first),*annotation_rows(first)],cls="small",caption="INFERRED — one bead-bound first strand after RT adds the CCC switch site.")}

<h3>(2) Template-switch, remove unused bead primers, and amplify cDNA</h3>
{panel(D.template_switch_scene().rows(),cls="small",caption="INFERRED — reconstructed TSO sequence; its rGrGrG pairs with CCC and supplies the second SMART handle.")}

<h3>(3) Tagment amplified cDNA; select the bead end with P5-SMART + N7xx PCR</h3>
{panel([strand_row(selected),*annotation_rows(selected)],cls="small",caption="INFERRED — reconstructed bead-end fragment. Only this end has the batch-matched handle/constant sequence for the selective P5 primer.")}
{panel([*Scene.duplex(list(lib),label="library").rows(),*annotation_rows(lib)],cls="small",caption="INFERRED — reconstructed batch-B final library, P5 to P7'. Samples are indexed on i7 only.")}
'''


def sequencing()->str:
    rows=[]
    for p,h in D.primer_landings(): rows.append((html.escape(p.role),html.escape(p.name),html.escape(", ".join(h.covers)),f"<code>{html.escape(h.reads)}</code>&hellip;"))
    return f'''<h2>Sequencing</h2>
{info('Read 1 reports the 12-base cell barcode + 8-base UMI; paired Read 2 is typically 50 bases of transcript sequence. The exact custom primer bases shown are inferred.')}
{table(("Read","Primer","Primer site","First bases"),rows)}
<h3>Read 1 layout</h3>{table(("Cycles","Content"),D.read1_layout(),scroll=False)}
'''


def main()->None:
    page="\n".join([head("Drop-seq library chemistry"),'<div class="wrap">','<h1>Drop-seq &mdash; barcoded-bead RNA capture in nanoliter droplets</h1>',info('Defining source: <a href="https://doi.org/10.1016/j.cell.2015.05.002">Macosko et al., <i>Cell</i> (2015)</a>.'),caveat('<b>Inferred oligo structure.</b> The defining paper’s oligo table was unavailable. Dotted regions reconstruct the published workflow using the secondary schematic and Seq-Well’s later primary-source table.'),oligos(),construction(),sequencing(),'</div>'])
    OUT.write_text(page,encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__=="__main__": main()
