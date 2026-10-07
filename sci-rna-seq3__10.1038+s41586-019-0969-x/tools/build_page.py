#!/usr/bin/env python3
"""Build the diagram-centric sci-RNA-seq3 chemistry page."""
from __future__ import annotations
import html,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE)); sys.path.insert(0,str(HERE.parents[1]/"lib"))
import sci_rna_seq3 as S
from chemdraw import Scene,annotation_rows,oligo,panel,strand_row
from page import head,info,table
OUT=HERE.parent/"sci-rna-seq3.html"


def oligos()->str:
    rows=[oligo("5'-phosphorylated RT primer",S.rt_primer(),mods="/5Phos/"),
          oligo("Barcoded dU hairpin adaptor",S.hairpin_oligo()),
          oligo("Indexed P5 PCR primer",S.p5_primer()),oligo("Indexed P7 PCR primer",S.p7_primer())]
    return "<h2>Key oligos</h2>\n<seq>\n"+"\n".join(rows)+"\n</seq>"


def construction()->str:
    first=S.first_strand(); lig=S.ligated_product(); opened=S.user_cleaved_tagmented(); lib=S.final_library()
    return f'''<h2>Library construction</h2>
<h3>(1) Add the RT-well barcode and UMI in situ</h3>
{panel(S.rt_scene().rows(),cls="small",caption="The 5'-phosphorylated primer begins with the six-base ligation site, then UMI, RT barcode and oligo-dT.")}
{panel([strand_row(first),*annotation_rows(first)],cls="small",caption="One barcoded first-strand cDNA.")}

<h3>(2) Ligate a second well-specific barcode through the hairpin splint</h3>
{panel(S.ligation_scene().rows(),cls="small",caption="The hairpin's GCTCTG overhang pairs with the RT primer's CAGAGC. Its barcode arms are reverse complements by construction.")}
{panel([strand_row(lig),*annotation_rows(lig)],cls="small",caption="Continuous ligated strand shown linearly; in solution the barcode arms form the hairpin stem around dU + the Read 1 handle.")}

<h3>(3) Make the second strand, tagment with N7-only Tn5, then open the hairpin with USER</h3>
{panel([strand_row(opened),*annotation_rows(opened)],cls="small",caption="USER removes dU and releases the short outer hairpin arm. The retained handle-to-s7 fragment carries both in-nucleus barcodes.")}

<h3>(4) Add the PCR-well i5/i7 pair</h3>
{panel([*Scene.duplex(list(lib),label="library").rows(),*annotation_rows(lib)],cls="small",caption="Final sci-RNA-seq3 library, P5 to P7'. Cell identity is ligation barcode + RT barcode + PCR-well i5/i7 pair.")}
'''


def sequencing()->str:
    rows=[]
    for p,h in S.primer_landings(): rows.append((html.escape(p.role),html.escape(p.name),html.escape(", ".join(h.covers)),f"<code>{html.escape(h.reads)}</code>&hellip;"))
    return f'''<h2>Sequencing</h2>
{info('Read 1: 34 cycles; Index 1: 10; Index 2: 10; Read 2: 52. Read 1 contains ligation barcode + CAGAGC + UMI + RT barcode; Read 2 contains transcript sequence.')}
{table(("Read","Primer","Primer site","First bases"),rows)}
<h3>Read 1 layout (representative 10-base ligation barcode)</h3>{table(("Cycles","Content"),S.read1_layout(),scroll=False)}
{info('The published adaptor set mixes 9- and 10-base ligation barcodes; a 9-base member shifts the following fields one cycle earlier.')}
'''


def main()->None:
    page="\n".join([head("sci-RNA-seq3 library chemistry"),'<div class="wrap">','<h1>sci-RNA-seq3 &mdash; three-level combinatorial indexing of nuclei</h1>',info('Defining source: <a href="https://doi.org/10.1038/s41586-019-0969-x">Cao et al., <i>Nature</i> (2019)</a>.'),oligos(),construction(),sequencing(),'</div>'])
    OUT.write_text(page,encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__=="__main__": main()
