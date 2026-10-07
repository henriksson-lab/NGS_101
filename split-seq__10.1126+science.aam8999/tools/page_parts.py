"""Individual sparse pages for SPLiT-seq and microSPLiT."""
from __future__ import annotations
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import splitseq as S
from chemdraw import Construct, Scene, Segment, annotation_rows, oligo, panel, strand_row
from page import head, info, table

PAPERS={"SPLiT-seq":("10.1126/science.aam8999","Rosenberg et al.","2018"),
        "microSPLiT":("10.1126/science.aba5257","Kuchina et al.","2021")}
def seg(n,t,g=None,**kw): return Segment(n,t,g,**kw)

def render_page(protocol):
    doi,author,year=PAPERS[protocol]; micro=protocol=="microSPLiT"
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
    return "\n".join([head(f"{protocol} library chemistry"),'<div class="wrap">',
      f'<h1>{protocol} &mdash; split-pool ligation barcoding in fixed cells</h1>',
      info(f'Defining source: <a href="https://doi.org/{doi}">{author}, <i>Science</i> ({year})</a>.'),
      f'<div class="caveat"><b>Sequence boundary.</b> {caveat}</div>',
      f'<h2>Key oligos</h2><seq>{"".join(oligos)}</seq>',
      '<h2>Barcode construction</h2>','<h3>(1) Barcoded in-cell reverse transcription</h3>',
      panel(S.rt_scene(protocol).rows(),cls="long",caption="Round 1 installs the first cell barcode."),
      panel([strand_row(b1),*annotation_rows(b1)],cls="long"),
      '<h3>(2) First split-pool ligation</h3>',
      (panel(S.splint_scene(2).rows(),cls="long",caption="The linker splints the round-2 oligo to the phosphorylated cDNA end.") if micro else ""),
      panel([strand_row(b2),*annotation_rows(b2)],cls="long",caption=("Round 2 adds barcode 2." if micro else "INFERRED — round-2 architecture; exact oligos unavailable.")),
      '<h3>(3) Second split-pool ligation</h3>',
      (panel(S.splint_scene(3).rows(),cls="long",caption="Round 3 adds the UMI and barcode 3; its biotin enables later capture.") if micro else ""),
      panel([strand_row(b3),*annotation_rows(b3)],cls="long",caption=("Three cell barcodes and the UMI now share one first strand." if micro else "INFERRED — round-3 architecture; exact oligos unavailable.")),
      '<h3>(4) Lyse, capture biotinylated cDNA, template-switch and amplify</h3>',
      panel([strand_row(Construct([seg("round-3 end", "[R3 / UMI / BC3 / BC2 / BC1]",placeholder=True),
                                  seg("cDNA", "XXXXXXXX...XXXXXXXX",placeholder=True),
                                  seg("template-switch handle", "[TSO HANDLE]", "tso",placeholder=True)]))],cls="long"),
      '<h3>(5) Build the indexed sequencing library</h3>',
      panel([strand_row(final),*annotation_rows(final)],cls="long",caption="INFERRED — supported inner barcode order; unpublished or unresolved outer library arms are dotted."),
      '<h2>Read layout</h2>',panel([strand_row(r1,prefix="Read 1  ",suffix=""),strand_row(r2,prefix="Read 2  ",suffix="")],cls="small",caption=("Published microSPLiT run: Read 1 74 nt, Read 2 86 nt, index 6 nt." if micro else "Published SPLiT-seq architecture: transcript in Read 1; combinatorial identity in Read 2.")),
      table(("Read","Reports"),rows),'</div>'])
