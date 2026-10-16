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
    return "<h2>Key oligos</h2>\n"+info("INFERRED for the defining experiment &mdash; Macosko and Goldman\u2019s inventor-lab Drop-seq protocol v3.1 (28 Dec 2015), required-primer table p. 3, prints these oligos and maps their names to Macosko et al.; Seq-Well Supplementary Table 1 corroborates the same bead family. Original Cell Table S6 remains inaccessible, so exact 2015 experiment identity is unverified and the bases are drawn dotted.")+"<seq>\n"+"\n".join(rows)+"\n</seq>"


def construction()->str:
    return '<h2>Library construction</h2>\n'+protocol_graph_panels(D.protocol_graph(),cls="long")


def sequencing()->str:
    lib=D.final_library()
    return f'''{sp.section(lib,D.seq_primers(),heading="Sequencing primers and reads",
            intro="Read 1 reports the 12-base cell barcode + 8-base UMI; paired Read 2 is typically 50 bases of transcript sequence. The inventor-lab v3.1 protocol prints the custom Read 1 primer; its transfer to the original 2015 experiment remains inferred until Table S6 is checked.",
            required_roles=tuple(p.role for p in D.seq_primers()),
            read_lengths=D.READ_LENGTHS,show_read_content=True)}
'''


def main()->None:
    page="\n".join([head("Drop-seq library chemistry"),'<div class="wrap" data-protocol-complete="true">','<h1>Drop-seq &mdash; barcoded-bead RNA capture in nanoliter droplets</h1>',info('Defining source: <a href="https://doi.org/10.1016/j.cell.2015.05.002">Macosko et al., <i>Cell</i> (2015)</a>.'),caveat('<b>Inferred oligo identity for the defining experiment.</b> The available Cell main article gives bead architecture but refers ordering oligos to its inaccessible Extended Procedures/Table S6. <a href="https://mccarrolllab.org/wp-content/uploads/2015/05/Online-Dropseq-Protocol-v.-3.1-Dec-2015.pdf">Macosko and Goldman’s inventor-lab v3.1 protocol</a> (28 Dec 2015, p. 3 required-primer table) prints bead, TSO, PCR and custom Read 1 strands, mapped by name to Macosko et al.; its opening notice says the optimized v3.1 workflow does not precisely match the Cell experiments. <a href="https://doi.org/10.1038/nmeth.4179">Seq-Well</a> later prints matching oligos for the same mRNA capture beads. The original Table S6 and its relationship to the later lab revision remain unverified, so source-transferred bases are dotted.'),oligos(),construction(),sequencing(),'</div>'])
    OUT.write_text(page,encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__=="__main__": main()
