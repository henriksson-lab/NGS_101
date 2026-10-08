#!/usr/bin/env python3
from __future__ import annotations
import html,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import tenx_crispr as C
import seqprimers as sp
from chemdraw import Scene,annotation_rows,oligo,panel
from page import head,info,table
OUT=HERE.parent/"10x-3prime-v3-crispr-screening.html"
def sequencing(lib):
    rows=[]
    for p in C.SEQ_PRIMERS:
        h=sp.locate(lib,p)
        if h is None: raise ValueError(f"{p.name} has no site")
        rows.append((html.escape(p.role),html.escape(p.name),html.escape(", ".join(h.covers)),f"<code>{html.escape(h.reads)}</code>&hellip;"))
    return '<h2>Sequencing</h2>'+info('Nextera Read 1 reports the 16-base cell barcode and 12-base UMI. TruSeq Read 2 crosses 30 TSO-derived bases, then reads the protospacer; i7 identifies the sample.')+sp.diagram(lib,C.SEQ_PRIMERS)+table(("Read","Primer","Primer site","First bases"),rows)+table(("Cycles","Content"),C.read_layout(),scroll=False)
def render():
    product=C.feature_cdna(); lib=C.final_library()
    fwd,rev=C.feature_cdna_primers(); si_fwd,si_rev=C.feature_si_primers()
    return "\n".join([head("10x 3′ gRNA capture with Feature Barcoding"),'<div class="wrap">','<h1>10x 3&#8242; Feature Barcode &mdash; direct gRNA capture</h1>','<p class="research-notes"><a href="01_10x-3p-crispr-screening.html">Research notes</a></p>',info('Commercial protocol: <a href="https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/CG000184_ChromiumSingleCellSingleCell3v3_FeatureBarcodingtechnology_CRISPR_RevA.pdf">10x Genomics CG000184 Rev A</a>; guide designs: <a href="https://cdn.10xgenomics.com/image/upload/v1660261286/support-documents/CG000197_GuideRNA_SpecificationsCompatible_withFeatureBarcodingtechnology_forCRISPRScreening_Rev-A.pdf">CG000197 Rev A</a>. This page shows the 3&#8242;-end CS1 guide configuration.'),
      '<h2>Capture oligos</h2><seq>'+oligo("Capture Sequence 1 gel-bead primer",C.bead_primer(),mods="bead")+oligo("template-switch oligo",C.tso(),mods="3' rGrGrG")+'</seq>',
      '<h2>Library construction</h2><h3>(1) Capture the engineered sgRNA inside each GEM</h3>'+panel(C.capture_scene().rows(),cls="small",caption="The guide carries CS1 reverse complement before its Pol III terminator. It pairs with the CS1 bead primer, which already carries the cell barcode and UMI."),
      '<h3>(2) Copy the sgRNA and template-switch</h3>'+panel(C.switch_scene().rows(),cls="small",caption="Reverse transcription copies the sgRNA feature, adds CCC and switches onto the TSO.")+panel([*Scene.duplex(list(product),label="feature cDNA").rows(),*annotation_rows(product)],cls="small",caption="Feature cDNA Primers 1 amplify between Partial Read 1N and the partial TSO."),
      '<h3>(3) Enrich the Feature Barcode cDNA</h3><seq>'+oligo("Feature cDNA forward",fwd)+oligo("Feature cDNA reverse",rev)+'</seq>',
      '<h3>(4) Add P5 and the Read 2 arm</h3><seq>'+oligo("Feature SI forward",si_fwd)+oligo("Feature SI reverse",si_rev)+'</seq>'+info('Feature SI reverse is TruSeq Read 2 followed by the partial TSO. It converts the template-switched end into a standard Read 2 end.'),
      '<h3>(5) Add P7 and the i7 sample index</h3>'+panel([*Scene.duplex(list(lib),label="library").rows(),*annotation_rows(lib)],cls="small",caption="Final Feature Barcode library. Read 1 identifies the cell and molecule; Read 2 identifies the guide."),sequencing(lib),'</div>'])
def main(): OUT.write_text(render(),encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")
if __name__=="__main__": main()
