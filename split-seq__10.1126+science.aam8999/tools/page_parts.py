"""Individual sparse pages for SPLiT-seq and microSPLiT."""
from __future__ import annotations
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import splitseq as S
import seqprimers as sp
from chemdraw import (Construct, Scene, Segment, annotation_rows, oligo, panel,
                      protocol_graph_panels, strand_row)
from page import head, info, table

PAPERS={"SPLiT-seq":("10.1126/science.aam8999","Rosenberg et al.","2018"),
        "microSPLiT":("10.1126/science.aba5257","Kuchina et al.","2021")}
def seg(n,t,g=None,**kw): return Segment(n,t,g,**kw)

def render_page(protocol):
    doi,author,year=PAPERS[protocol]; micro=protocol=="microSPLiT"
    notes=("01_split-seq.html" if not micro else
           "../split-seq__10.1126+science.aam8999/01_split-seq.html")
    caveat=("The published oligo supplement was unavailable. Dotted protocol-specific "
            "regions show architecture only; upstream sequences are not promoted to fact."
            if not micro else
            "Barcoding oligos come from the authors' preprint table. The published downstream "
            "library supplement was unavailable, so outer library arms stop at that boundary.")
    oligos=[oligo("round-1 oligo-dT primer",S.round1(protocol))]
    if micro:
        oligos += [oligo("round-2 barcode oligo",S.round2(),mods="/5Phos/"),
                   oligo("round-3 barcode oligo",S.round3(),mods="/5Biosg/"),
                   oligo("template-switching oligo",[seg("",S.TSO_BODY,"tso"),
                                                       seg("","rGrG+G",placeholder=True)])]
    b1=S.barcoded_first_strand(protocol,1); b2=S.barcoded_first_strand(protocol,2)
    b3=S.barcoded_first_strand(protocol,3); final=S.final_boundary(protocol)
    r1=Construct([seg("transcript", "XXXXXXXX...",placeholder=True)])
    r2=Construct([seg("UMI", "N"*10,"umi",placeholder=True),
                  seg("barcode 3", "N"*8,"cbc",placeholder=True),
                  seg("link/barcode 2/barcode 1", "[LINK / BC2 / BC1]", "cbc",
                      placeholder=True,inferred=not micro)])
    rows=[("Read 1","transcript sequence"),
          ("Read 2","UMI, barcode 3, barcode 2, barcode 1"),
          ("Index","sublibrary barcode (round 4)")]
    graph=S.protocol_graph(protocol)
    return "\n".join([head(f"{protocol} library chemistry"),'<div class="wrap" data-protocol-complete="true">',
      f'<h1>{protocol} &mdash; split-pool ligation barcoding in fixed cells</h1>',
      f'<p class="research-notes"><a href="{notes}">Research notes</a></p>',
      info(f'Defining source: <a href="https://doi.org/{doi}">{author}, <i>Science</i> ({year})</a>.'),
      f'<div class="caveat"><b>Sequence boundary.</b> {caveat}</div>',
      f'<h2>Key oligos</h2><seq>{"".join(oligos)}</seq>',
      '<h2>Reaction workflow</h2>',protocol_graph_panels(graph,cls="long"),
      sp.unavailable_diagram(final,"the downstream library-primer sequences are unavailable"),
      '<h2>Read layout</h2>',panel([strand_row(r1,prefix="Read 1  ",suffix=""),strand_row(r2,prefix="Read 2  ",suffix="")],cls="small",caption=("Published microSPLiT run: Read 1 74 nt, Read 2 86 nt, index 6 nt." if micro else "Published SPLiT-seq architecture: transcript in Read 1; combinatorial identity in Read 2.")),
      table(("Read","Reports"),rows),
      sp.read_spec_table([
          sp.read_spec("Read 1", 74 if micro else None, "transcript sequence", "include"),
          sp.read_spec("Read 2", 86 if micro else None,
                       "UMI, barcode 3, barcode 2, barcode 1", "extract"),
          sp.read_spec("Index 1 (i7)", 6 if micro else None,
                       "sublibrary barcode (round 4)", "demultiplex"),
      ], "Read roles and contents are source-declared because exact downstream primer sequences are unavailable."),
      '</div>'])
