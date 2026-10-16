"""Individual sparse pages for SPLiT-seq and microSPLiT."""
from __future__ import annotations
import sys
import html
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

def render_page(protocol, *, graph_factory=None):
    doi,author,year=PAPERS[protocol]; micro=protocol=="microSPLiT"
    if micro:
        # Keep the published adapter transcription in the microSPLiT module.
        from microsplit_graph import BC_0243_CORE, BC_0243_TOP, BC_0244_BOTTOM
    notes=("01_split-seq.html" if not micro else
           "../split-seq__10.1126+science.aam8999/01_split-seq.html")
    caveat=("Supplementary Table S12 supplies the barcoding oligos, splints, blocking "
            "strands and PCR primers. The downstream library and primer-binding geometry "
            "has not yet been reconstructed in this model from that table."
            if not micro else
            "A mirrored copy of the authors' 2021 Supplementary Table S3 supplies the "
            "barcoding oligos, adapter duplex and terminal indexed PCR pair. The later 2024 "
            "protocol supplies BC_0062/BC_0108 cDNA-amplification primers. The original "
            "publisher Table S3 attachment remains inaccessible for byte comparison. "
            "TSBC07 is the representative sublibrary index below.")
    oligos=[oligo("round-1 oligo-dT primer",S.round1(protocol),mods="/5Phos/")]
    if not micro:
        oligos += [oligo("round-2 barcode oligo",S.round2(),mods="/5Phos/"),
                   oligo("round-3 barcode oligo",S.round3(),mods="/5Biosg/"),
                   oligo("BC_0215 round-2 splint",[seg("",S.SPLIT_LINK2)]),
                   oligo("BC_0216 round-2 blocker",[seg("",S.revcomp(S.SPLIT_LINK2))]),
                   oligo("BC_0060 round-3 splint",[seg("",S.SPLIT_LINK3)]),
                   oligo("BC_0066 round-3 blocker",[seg("",S.revcomp(S.SPLIT_LINK3))]),
                   oligo("BC_0127 template-switch oligo",[seg("",S.TSO_BODY,"tso"),
                                                              seg("","rGrG+G",placeholder=True)]),
                   oligo("BC_0062 cDNA PCR",[seg("",S.MICRO_PREAMP_FWD)]),
                   oligo("BC_0108 cDNA PCR",[seg("",S.MICRO_PREAMP_REV)]),
                   oligo("BC_0118 library PCR",[seg("",S.SPLIT_LIBRARY_FWD)]),
                   oligo("BC_0076 / TSBC07 library PCR",[seg("",S.SPLIT_LIBRARY_REV)])]
    if micro:
        oligos += [oligo("round-2 barcode oligo",S.round2(),mods="/5Phos/"),
                   oligo("round-3 barcode oligo",S.round3(),mods="/5Biosg/"),
                   oligo("template-switching oligo",[seg("",S.TSO_BODY,"tso"),
                                                       seg("","rGrG+G",placeholder=True)]),
                   oligo("BC_0062 cDNA PCR",[seg("",S.MICRO_PREAMP_FWD)]),
                   oligo("BC_0108 cDNA PCR",[seg("",S.MICRO_PREAMP_REV)]),
                   oligo("BC_0243 adapter top — terminal C*T phosphorothioate",
                         [seg("",BC_0243_TOP)]),
                   oligo("BC_0244 adapter bottom, written 5′ to 3′",
                         [seg("",BC_0244_BOTTOM)]),
                   oligo("BC_0027 library PCR",[seg("",S.MICRO_LIBRARY_FWD)]),
                   oligo("BC_0076 / TSBC07 library PCR",[seg("",S.MICRO_LIBRARY_REV)])]
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
    graph=graph_factory() if graph_factory else S.graph_for(protocol)
    return "\n".join([head(f"{protocol} library chemistry"),'<div class="wrap" data-protocol-complete="true">',
      f'<h1>{protocol} &mdash; split-pool ligation barcoding in fixed cells</h1>',
      f'<p class="research-notes"><a href="{notes}">Research notes</a></p>',
      info(f'Defining source: <a href="https://doi.org/{doi}">{author}, <i>Science</i> ({year})</a>.'),
      f'<div class="caveat"><b>Sequence boundary.</b> {caveat}</div>',
      (f'<div class="caveat"><b>Source-version difference.</b> The mirrored 2021 '
       f'<i>Science</i> Table S3, General Oligonucleotides row BC_0243, prints '
       f'<code>{html.escape(BC_0243_CORE + "*T")}</code> and labels the asterisk '
       f'“phosphorothioated base.” The 2019 bioRxiv v2 Auxiliary Table line 14 and '
       f'2024 <i>Nature Protocols</i> Table 3 line 121 print the same bases as '
       f'<code>{html.escape(BC_0243_TOP)}</code> without the mark. This model uses the '
       f'2021 Table S3 terminal C*T modification. Whether the difference reflects a '
       f'version change or omitted notation is unresolved; the official 2021 publisher '
       f'attachment has not been compared byte for byte.</div>' if micro else ''),
      f'<h2>Key oligos</h2><seq>{"".join(oligos)}</seq>',
      '<h2>Reaction workflow</h2>',protocol_graph_panels(graph,cls="long"),
      (sp.diagram(final,S.MICRO_SEQ_PRIMERS) if micro else
       sp.unavailable_diagram(final,"the downstream primer-binding geometry has not yet been reconstructed")),
      '<h2>Read layout</h2>',panel([strand_row(r1,prefix="Read 1  ",suffix=""),strand_row(r2,prefix="Read 2  ",suffix="")],cls="small",caption=("Published microSPLiT run: Read 1 74 nt, Read 2 86 nt, index 6 nt." if micro else "Published SPLiT-seq architecture: transcript in Read 1; combinatorial identity in Read 2.")),
      table(("Read","Reports"),rows),
      sp.read_spec_table([
          sp.read_spec("Read 1", 74 if micro else None, "transcript sequence", "include"),
          sp.read_spec("Read 2", 86 if micro else None,
                       "UMI, barcode 3, barcode 2, barcode 1", "extract"),
          sp.read_spec("Index 1 (i7)", 6 if micro else None,
                       "sublibrary barcode (round 4)", "demultiplex"),
      ], ("Read contents and primer landing sites are derived from Supplementary Table S3."
          if micro else "Read roles and contents are source-declared; Table S12 prints the PCR oligos, but this model has not reconstructed the exact primer-binding geometry.")),
      '</div>'])
