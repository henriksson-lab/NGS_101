#!/usr/bin/env python3
from __future__ import annotations
import html,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import tenx_crispr as C
import seqprimers as sp
from chemdraw import Scene,annotation_rows,oligo,panel,protocol_graph_panels
from page import head,info,table
from protocol_graph import AnnotationSlot
from protocol_views import assemble_page
OUT=HERE.parent/"10x-3prime-v3-crispr-screening.html"
def sequencing(graph):
    plan = graph.sequencing["library"]
    lib = plan.library
    rows=[]
    for p in plan.primers:
        h=sp.locate(lib,p)
        if h is None: raise ValueError(f"{p.name} has no site")
        rows.append((html.escape(p.role),html.escape(p.name),html.escape(", ".join(h.covers)),f"<code>{html.escape(h.reads)}</code>&hellip;"))
    return '<h2>Sequencing</h2>'+info(f'Nextera Read 1 reports the 16-base cell barcode and 12-base UMI. TruSeq Read 2 crosses 30 TSO-derived bases, then reads the protospacer; i7 identifies the sample. CG000184 Rev A recommends 28/8/0/91 cycles (Read 1/i7/i5/Read 2) and states a minimum Read 2 length of {C.READ2_MINIMUM} bp for CRISPR Screening libraries.')+sp.diagram(lib,plan.primers)+table(("Read","Primer","Primer site","First bases"),rows)+table(("Cycles","Content"),C.read_layout(),scroll=False)+sp.read_content_table(lib,plan.primers,plan.read_lengths,acquired_roles=plan.acquired_roles)
def render():
    product=C.amplified_cdna(); lib=C.final_library(); graph=C.protocol_graph()
    fwd,rev=C.feature_cdna_primers(); si_fwd,si_rev=C.feature_si_primers()
    return assemble_page(graph,[head("10x 3′ gRNA capture with Feature Barcoding"),'<div class="wrap" data-protocol-complete="true">','<h1>10x 3&#8242; Feature Barcode &mdash; direct gRNA capture</h1>','<p class="research-notes"><a href="01_10x-3p-crispr-screening.html">Research notes</a></p>',info('Commercial protocol: <a href="https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/CG000184_ChromiumSingleCellSingleCell3v3_FeatureBarcodingtechnology_CRISPR_RevA.pdf">10x Genomics CG000184 Rev A</a>; guide designs: <a href="https://cdn.10xgenomics.com/image/upload/v1660261286/support-documents/CG000197_GuideRNA_SpecificationsCompatible_withFeatureBarcodingtechnology_forCRISPRScreening_Rev-A.pdf">CG000197 Rev A</a>. This page shows the 3&#8242;-end CS1 guide configuration.'),AnnotationSlot.INTRO,
      '<h2>Reaction workflow</h2>'+protocol_graph_panels(graph,cls="long"),
      '<h2>Capture oligos</h2><seq>'+oligo("Capture Sequence 1 gel-bead primer",C.bead_primer(),mods="bead")+oligo("template-switch oligo",C.tso(),mods="3' rGrGrG")+'</seq>',
      '<h2>Library construction</h2><h3>(1) Capture the engineered sgRNA inside each GEM</h3>'+panel(C.capture_scene().rows(),cls="small",caption="The guide carries CS1 reverse complement before its Pol III terminator. It pairs with the CS1 bead primer, which already carries the cell barcode and UMI."),
      '<h3>(2) Copy the sgRNA and template-switch</h3>'+panel(C.switch_scene().rows(),cls="small",caption="Reverse transcription copies the sgRNA feature, adds CCC and switches onto the TSO.")+panel([*Scene.duplex(list(product),label="feature cDNA").rows(),*annotation_rows(product)],cls="small",caption="Feature cDNA Primers 1 (PN-2000096) amplify between Partial Read 1N and the partial TSO; 0.6X SPRIselect leaves this short product in the supernatant (recovered at 1.2X, then 1.0X Guide RNA cDNA Cleanup)."),
      '<h3>(3) Feature cDNA Primers 1 used in cDNA amplification</h3><seq>'+oligo("Feature cDNA forward",fwd)+oligo("Feature cDNA reverse",rev)+'</seq>',
      '<h3>(4) Feature PCR adds P5 and the Read 2 arm (Feature SI Primers 1)</h3><seq>'+oligo("Feature SI forward",si_fwd)+oligo("Feature SI reverse",si_rev)+'</seq>'+info('Feature SI reverse is TruSeq Read 2 followed by the partial TSO. It converts the template-switched end into a standard Read 2 end. A 1.0X SPRIselect cleanup follows.'),
      '<h3>(5) Sample Index PCR adds P7 and the i7 index (Feature SI Primers 2 + Chromium i7)</h3>'+panel([*Scene.duplex(list(lib),label="library").rows(),*annotation_rows(lib)],cls="small",caption="Final Feature Barcode library after 0.7X/1.0X double-sided SPRIselect. Read 1 identifies the cell and molecule; Read 2 identifies the guide."),sequencing(graph),'</div>'])
def main(): OUT.write_text(render(),encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")
if __name__=="__main__": main()
