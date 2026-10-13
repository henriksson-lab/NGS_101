#!/usr/bin/env python3
"""Build the diagram-centric sci-RNA-seq3 chemistry page."""
from __future__ import annotations
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE)); sys.path.insert(0,str(HERE.parents[1]/"lib"))
import sci_rna_seq3 as S
import seqprimers as sp
from chemdraw import oligo,protocol_graph_panels
from page import head,info
OUT=HERE.parent/"sci-rna-seq3.html"


def oligos()->str:
    rows=[oligo("5'-phosphorylated RT primer",S.rt_primer(),mods="/5Phos/"),
          oligo("Barcoded dU hairpin adaptor",S.hairpin_oligo()),
          oligo("Indexed P5 PCR primer",S.p5_primer()),oligo("Indexed P7 PCR primer",S.p7_primer())]
    return "<h2>Key oligos</h2>\n<seq>\n"+"\n".join(rows)+"\n</seq>"


def construction()->str:
    return (S.HAIRPIN_TOPOLOGY.folded_panel(
        "The published adaptor folded base by base. Its ten-base ligation barcode and "
        "reverse complement form the stem; dU plus the Read 1 handle form the loop. "
        "The six-base 5′ overhang remains outside the stem and splints ligation to the "
        "RT-primer end.")
        + '<h2>Library construction</h2>\n'
        + protocol_graph_panels(S.protocol_graph(),cls="long"))


def sequencing()->str:
    lib=S.final_library()
    return (sp.section(
        lib,S.SEQ_PRIMERS,heading="Sequencing",
        intro=("Read 1 contains the ligation barcode, UMI and RT barcode; Read 2 "
               "contains transcript sequence. The two index reads identify the PCR well."),
        read_lengths=S.READ_LENGTHS,read_uses=S.READ_USES,show_read_content=True)
        +info('The published adaptor set mixes 9- and 10-base ligation barcodes; a 9-base member shifts the following fields one cycle earlier.'))


def main()->None:
    page="\n".join([head("sci-RNA-seq3 library chemistry"),'<div class="wrap" data-protocol-complete="true">','<h1>sci-RNA-seq3 &mdash; three-level combinatorial indexing of nuclei</h1>',info('Defining source: <a href="https://doi.org/10.1038/s41586-019-0969-x">Cao et al., <i>Nature</i> (2019)</a>.'),oligos(),construction(),sequencing(),'</div>'])
    OUT.write_text(page,encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__=="__main__": main()
