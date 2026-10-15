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
    rows=[oligo("Barcoded bead oligo",D.bead_primer(),five="5'-Bead-Linker-"),
          oligo("Template-switching oligo",D.tso(),three="-3' (terminal rGrGrG)"),
          oligo("P5-SMART selective PCR primer",D.p5_hybrid(),three="-*A*C-3'"),
          oligo("Custom Read 1 primer",[D.inferred("spacer + handle + AC",D.custom_r1(),"r1")])]
    return "<h2>Key oligos</h2>\n"+info("INFERRED &mdash; the bases are the ones Seq-Well\u2019s Supplementary Table 1 prints for the same ChemGenes MACOSKO-2011-10 bead; Drop-seq\u2019s own Table S6 is unavailable, so their transfer to Drop-seq is an inference and is drawn dotted.")+"<seq>\n"+"\n".join(rows)+"\n</seq>"


def construction()->str:
    return '<h2>Library construction</h2>\n'+protocol_graph_panels(D.protocol_graph(),cls="long")


def sequencing()->str:
    lib=D.final_library()
    return f'''{sp.section(lib,D.seq_primers(),heading="Sequencing primers and reads",
            intro="Read 1 reports the 12-base cell barcode + 8-base UMI; paired Read 2 is typically 50 bases of transcript sequence. The custom Read 1 primer bases come from Seq-Well\u2019s oligo table for the same bead and are drawn as inferred for Drop-seq.",
            required_roles=tuple(p.role for p in D.seq_primers()),
            read_lengths=D.READ_LENGTHS,show_read_content=True)}
'''


def main()->None:
    page="\n".join([head("Drop-seq library chemistry"),'<div class="wrap" data-protocol-complete="true">','<h1>Drop-seq &mdash; barcoded-bead RNA capture in nanoliter droplets</h1>',info('Defining source: <a href="https://doi.org/10.1016/j.cell.2015.05.002">Macosko et al., <i>Cell</i> (2015)</a>.'),caveat('<b>Inferred oligo structure.</b> The defining paper prints no oligo: its Table S6 is behind the PMC download gate. Dotted regions carry the bases that the later Seq-Well paper (doi:10.1038/nmeth.4179) prints in its own Supplementary Table 1 for the same ChemGenes MACOSKO-2011-10 bead; that transfer is an inference, not a Drop-seq statement.'),oligos(),construction(),sequencing(),'</div>'])
    OUT.write_text(page,encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__=="__main__": main()
