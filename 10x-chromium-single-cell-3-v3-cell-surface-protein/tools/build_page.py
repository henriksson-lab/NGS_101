#!/usr/bin/env python3
from __future__ import annotations
import html,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import tenx_csp as C
import seqprimers as sp
from chemdraw import Scene,annotation_rows,oligo,panel,protocol_graph_panels
from page import head,info,table
OUT=HERE.parent/"10x-3prime-v3-cell-surface-protein.html"
def sequencing(graph):
    plan = graph.sequencing["library"]
    lib = plan.library
    rows=[]
    for p in plan.primers:
        h=sp.locate(lib,p)
        if h is None: raise ValueError(f"{p.name} has no site")
        rows.append((html.escape(p.role),html.escape(p.name),html.escape(", ".join(h.covers)),f"<code>{html.escape(h.reads)}</code>&hellip;"))
    return '<h2>Sequencing</h2>'+info(f'Nextera Read 1 reports the 16-base cell barcode and 12-base UMI. TruSeq Read 2 reaches the 15-base antibody feature barcode after 10 diversity bases; i7 identifies the sample. CG000185 Rev B recommends 28/8/0/91 cycles (Read 1/i7/i5/Read 2) and states a minimum Read 2 length of {C.READ2_MINIMUM} bp for Cell Surface Protein libraries; at 91 cycles Read 2 runs on through Capture Sequence 1 into the UMI and cell barcode (computed).')+sp.diagram(lib,plan.primers)+table(("Read","Primer","Primer site","First bases"),rows)+table(("Cycles","Content"),C.read_layout(),scroll=False)+sp.read_content_table(lib,plan.primers,plan.read_lengths,acquired_roles=plan.acquired_roles)
def render():
    product=C.amplified_feature(); lib=C.final_library(); graph=C.protocol_graph()
    return "\n".join([head("10x 3' v3 Cell Surface Protein chemistry"),'<div class="wrap" data-protocol-complete="true">','<h1>10x Chromium Single Cell 3&#39; v3 &mdash; Cell Surface Protein</h1>',info(graph.authority),
      '<h2>Reaction workflow</h2>'+protocol_graph_panels(graph,cls="long"),
      '<h2>Key oligos</h2><seq>'+oligo("Capture Sequence 1 gel-bead primer",C.bead_primer(),mods="bead")+oligo("antibody Feature Barcode oligo",C.antibody_oligo(),mods="5' antibody conjugate")+'</seq>',
      '<h2>Library construction</h2><h3>(1) Capture antibody tags inside each GEM</h3>'+panel(C.capture_scene().rows(),cls="small",caption="Capture Sequence 1 pairs the antibody tag with a bead primer carrying the cell barcode and UMI; extension copies the feature barcode and Read 2 handle."),
      '<h3>(2) Co-amplify the short feature cDNA, then size-separate it from gene-expression cDNA</h3>'+panel([*Scene.duplex(list(product),label="feature cDNA").rows(),*annotation_rows(product)],cls="small",caption="Feature cDNA Primers 2 (PN-2000097) amplify between Partial Read 1N and Partial Read 2; a 0.6X SPRIselect leaves this short product in the supernatant, which is recovered with 2.0X SPRIselect."),
      '<h3>(3) Add P5, P7 and the sample index by Feature SI PCR</h3>'+panel([*Scene.duplex(list(lib),label="library").rows(),*annotation_rows(lib)],cls="small",caption="Final Cell Surface Protein library from CG000185 Rev B appendix Step 4.1 (Feature SI Primers 2 PN-2000099 + Chromium i7 Sample Index)."),sequencing(graph),'</div>'])
def main(): OUT.write_text(render(),encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")
if __name__=="__main__": main()
