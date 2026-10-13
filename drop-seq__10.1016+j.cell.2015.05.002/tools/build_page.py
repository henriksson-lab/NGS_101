#!/usr/bin/env python3
"""Build the diagram-centric Drop-seq chemistry page."""
from __future__ import annotations
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE)); sys.path.insert(0,str(HERE.parents[1]/"lib"))
import drop_seq as D
import seqprimers as sp
from chemdraw import Segment,oligo,protocol_graph_panels
from page import caveat,head,info
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
    return '<h2>Library construction</h2>\n'+protocol_graph_panels(D.protocol_graph(),cls="long")


def sequencing()->str:
    lib=D.final_library()
    return f'''{sp.section(lib,D.seq_primers(),heading="Sequencing primers and reads",
            intro="Read 1 reports the 12-base cell barcode + 8-base UMI; paired Read 2 is typically 50 bases of transcript sequence. The exact custom primer bases shown are inferred.",
            required_roles=tuple(p.role for p in D.seq_primers()),
            read_lengths=D.READ_LENGTHS,show_read_content=True)}
'''


def main()->None:
    page="\n".join([head("Drop-seq library chemistry"),'<div class="wrap" data-protocol-complete="true">','<h1>Drop-seq &mdash; barcoded-bead RNA capture in nanoliter droplets</h1>',info('Defining source: <a href="https://doi.org/10.1016/j.cell.2015.05.002">Macosko et al., <i>Cell</i> (2015)</a>.'),caveat('<b>Inferred oligo structure.</b> The defining paper’s oligo table was unavailable. Dotted regions reconstruct the published workflow using the secondary schematic and Seq-Well’s later primary-source table.'),oligos(),construction(),sequencing(),'</div>'])
    OUT.write_text(page,encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__=="__main__": main()
