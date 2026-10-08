#!/usr/bin/env python3
"""Build the diagram-centric sci-RNA-seq chemistry page."""
from __future__ import annotations
import html
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parents[1] / "lib"))
import sci_rna_seq as S
import seqprimers as sp
from chemdraw import Scene, annotation_rows, oligo, panel, strand_row
from page import head, info, table

OUT = HERE.parent / "sci-rna-seq.html"


def oligos() -> str:
    rows = [oligo("Barcoded anchored oligo-dT RT primer", S.rt_primer()),
            oligo("Indexed P5 PCR primer", S.p5_primer()),
            oligo("Indexed P7 PCR primer", S.p7_primer())]
    return "<h2>Key oligos</h2>\n<seq>\n" + "\n".join(rows) + "\n</seq>"


def construction() -> str:
    first = S.first_strand(); tagged = S.selected_tagmented_strand(); lib = S.final_library()
    return f'''<h2>Library construction</h2>
<h3>(1) Barcode first-strand cDNA by in-situ reverse transcription</h3>
{panel(S.rt_scene().rows(), cls="small", caption="Each RT well contributes a 10-base cell barcode; the adjacent 8-base UMI labels the captured transcript molecule.")}
{panel([strand_row(first), *annotation_rows(first)], cls="small", caption="One first-strand cDNA after reverse transcription.")}

<h3>(2) Pool cells, sort into PCR wells, and synthesize the second strand</h3>
{panel(S.second_strand_scene().rows(), cls="small", caption="The fixed cell remains the reaction vessel while the RT barcode stays attached to its cDNA.")}

<h3>(3) Tagment; retain the RT-primer end nearest an s7 insertion</h3>
{panel([strand_row(tagged), *annotation_rows(tagged)], cls="small", caption="Only a handle-to-s7 fragment has both PCR primer sites. Other Tn5 products and carrier DNA lack the RT handle.")}

<h3>(4) Add the PCR-well i5/i7 pair</h3>
{panel([*Scene.duplex(list(lib), label="library").rows(), *annotation_rows(lib)], cls="small", caption="Final two-level sci-RNA-seq library, P5 to P7'. Cell identity is RT barcode + PCR-well i5/i7 pair.")}
'''


def sequencing() -> str:
    lib=S.final_library()
    rows=[]
    for p,h in S.primer_landings():
        rows.append((html.escape(p.role),html.escape(p.name),html.escape(", ".join(h.covers)),f"<code>{html.escape(h.reads)}</code>&hellip;"))
    return f'''<h2>Sequencing</h2>
{info('Read 1: 18 cycles; Index 1: 10; Index 2: 10; Read 2: 52. Read 1 contains UMI + RT barcode; Read 2 contains transcript sequence.')}
{sp.diagram(lib,S.SEQ_PRIMERS)}
{table(("Read","Primer","Primer site","First bases"),rows)}
<h3>Read 1 barcode layout</h3>{table(("Cycles","Content"),S.read1_layout(),scroll=False)}
'''


def main() -> None:
    page="\n".join([head("sci-RNA-seq library chemistry"),'<div class="wrap">','<h1>sci-RNA-seq &mdash; two-level combinatorial indexing of fixed cells</h1>',info('Defining source: <a href="https://doi.org/10.1126/science.aam8940">Cao et al., <i>Science</i> (2017)</a>.'),oligos(),construction(),sequencing(),'</div>'])
    OUT.write_text(page,encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__": main()
