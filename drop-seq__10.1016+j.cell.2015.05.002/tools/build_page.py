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
    rows=[oligo("Barcoded Bead SeqB (drawn)",D.bead_primer("SeqB"),five="5'-Bead-Linker-"),
          oligo("Barcoded Bead SeqA",D.bead_primer("SeqA"),five="5'-Bead-Linker-"),
          oligo("Template_Switch_Oligo",D.tso(),three="-3' (terminal rGrGrG)"),
          oligo("TSO_PCR",D.tso_pcr()),
          oligo("P5-TSO_Hybrid",D.p5_hybrid(),three="-3' (*A*C phosphorothioates)"),
          oligo("Nextera_N701",D.n701()),
          oligo("Read1CustomSeqB",[sg("spacer + handle + AC",D.custom_r1("SeqB"),"r1")]),
          oligo("Read1CustomSeqA",[sg("spacer + handle + ACGT",D.custom_r1("SeqA"),"r1")])]
    return "<h2>Key oligos</h2>\n"+info("Verbatim from Macosko et al. 2015 Supplemental Information, Table S6. Nextera_N702 and N703 differ from N701 only in the i7 window and were used for the three-condition species-contamination library.")+"<seq>\n"+"\n".join(rows)+"\n</seq>"


def discrepancies()->str:
    return ("<h2>Published discrepancies</h2>\n"+info(
        "<b>Read 1 versus read 2 for SeqB.</b> The Supplemental sequencing paragraph (Preparation of Drop-Seq cDNA Library for Sequencing) says Read1CustSeqA or Read1CustSeqB is used &ldquo;for priming of read 1&rdquo;, but the closing note on bead surface primers says &ldquo;to prime read 2 for Drop-Seq libraries produced using Barcoded Bead SeqB beads, Read1CustSeqB was used&rdquo;. Read1CustomSeqB ends on the SeqB bead constant and lands immediately before the cell barcode, so this page uses it as the Read 1 primer; the &ldquo;read 2&rdquo; wording is treated as a slip, which the authors have not confirmed.<br>"
        "<b>Oligo names.</b> The Supplemental Experimental Procedures call the oligos Template_Switch_PCR primer, P5_TSO_Hybrid and Read1CustSeqA/B, while Table S6 heads them TSO_PCR, P5-TSO_Hybrid and Read1CustomSeqA/B. Each name pair has one printed sequence; this page uses the Table S6 names."))


def construction()->str:
    return '<h2>Library construction</h2>\n'+protocol_graph_panels(D.protocol_graph(),cls="long")


def sequencing()->str:
    lib=D.final_library()
    alt=D.final_library(variant="SeqA")
    return f'''{sp.section(lib,D.seq_primers(),heading="Sequencing primers and reads",
            intro="NextSeq 500 with Read1CustomSeqB (SeqB beads) or Read1CustomSeqA (SeqA beads) as the custom Read 1 primer. Read 1 was 20 bp (bases 1-12 cell barcode, 13-20 UMI); paired Read 2 was 50 bp for the human-mouse experiment and 60 bp for the retina experiment (Supplemental Experimental Procedures). The Supplemental Information does not state the index-read length; the 8-base Index 1 shown here follows the inventor-lab Drop-seq v3.1 specification and is needed only when libraries are multiplexed.",
            required_roles=tuple(p.role for p in D.seq_primers()),
            read_lengths=D.READ_LENGTHS,show_read_content=True)}
{sp.section(alt,D.seq_primers("SeqA"),heading="SeqA libraries (human-mouse and retina replicates 1-3)",
            intro="The SeqA bead carries ACGT after the SMART handle; Read1CustomSeqA ends on that ACGT, so Read 1 again starts on the first cell-barcode base. P5-TSO_Hybrid ends at the shared AC and amplifies both bead batches.",
            required_roles=tuple(p.role for p in D.seq_primers("SeqA")),
            read_lengths=D.READ_LENGTHS,show_read_content=True)}
'''


def main()->None:
    page="\n".join([head("Drop-seq library chemistry"),'<div class="wrap" data-protocol-complete="true">','<h1>Drop-seq &mdash; barcoded-bead RNA capture in nanoliter droplets</h1>',info('Defining source: <a href="https://doi.org/10.1016/j.cell.2015.05.002">Macosko et al., <i>Cell</i> (2015)</a>.'),info('Oligos and conditions: the paper\u2019s Supplemental Information (Table S6 and Supplemental Experimental Procedures). Two bead batches were used; the construction diagram draws Barcoded Bead SeqB.'),caveat('<b>Inferred bases.</b> Only the untemplated CCC added by the reverse transcriptase (paired by the TSO rGrGrG) is inferred from the standard template-switching mechanism; its length is not printed. Steps carrying it are marked INFERRED.'),oligos(),discrepancies(),construction(),sequencing(),'</div>'])
    OUT.write_text(page,encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__=="__main__": main()
