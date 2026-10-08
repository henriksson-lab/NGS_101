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
    return '<h2>Sequencing</h2>'+info('Nextera Read 1 reports the 16-base cell barcode and 12-base UMI. TruSeq Read 2 crosses the TSO-derived segment before reaching the sgRNA feature; i7 identifies the sample.')+sp.diagram(lib,C.SEQ_PRIMERS)+table(("Read","Primer","Primer site","First bases"),rows)+table(("Cycles","Content"),C.read_layout(),scroll=False)
def render():
    product=C.feature_cdna(); lib=C.final_library()
    return "\n".join([head("10x 3' v3 CRISPR Screening chemistry"),'<div class="wrap">','<h1>10x Chromium Single Cell 3&#39; v3 &mdash; CRISPR Screening</h1>',info('Primary chemistry source: <a href="https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/CG000184_ChromiumSingleCellSingleCell3v3_FeatureBarcodingtechnology_CRISPR_RevA.pdf">10x Genomics User Guide CG000184 Rev A</a>.'),
      '<div class="caveat"><b>Final-product drawing.</b> CG000184 prints the TSO-derived segment in its Feature PCR product, then omits it from the following Sample Index PCR sequence line. This schematic carries the existing segment through the indexing PCR.</div>',
      '<h2>Key oligos</h2><seq>'+oligo("Capture Sequence 1 gel-bead primer",C.bead_primer(),mods="bead")+oligo("template-switch oligo",C.tso(),mods="3' rGrGrG")+'</seq>',
      '<h2>Library construction</h2><h3>(1) Capture an engineered sgRNA inside each GEM</h3>'+panel(C.capture_scene().rows(),cls="small",caption="Capture Sequence 1 pairs the sgRNA with a bead primer carrying the cell barcode and UMI."),
      '<h3>(2) Copy the sgRNA and template-switch</h3>'+panel(C.switch_scene().rows(),cls="small",caption="Reverse transcription copies the sgRNA feature, adds CCC and switches onto the TSO.")+panel([*Scene.duplex(list(product),label="feature cDNA").rows(),*annotation_rows(product)],cls="small",caption="Feature cDNA Primers 1 amplify between Partial Read 1N and the partial TSO."),
      '<h3>(3) Add P5 and the Read 2 arm by Feature PCR</h3>'+info('The source-listed Feature SI reverse primer is TruSeq Read 2 followed by the partial TSO, so the TSO end becomes a standard Read 2 end.'),
      '<h3>(4) Add P7 and the i7 sample index</h3>'+panel([*Scene.duplex(list(lib),label="library").rows(),*annotation_rows(lib)],cls="small",caption="Final CRISPR Screening library, with the sgRNA feature between Capture Sequence 1 and the TSO-derived segment."),sequencing(lib),'</div>'])
def main(): OUT.write_text(render(),encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")
if __name__=="__main__": main()
