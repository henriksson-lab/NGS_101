#!/usr/bin/env python3
"""Build the diagram-centric Seq-Well S3 chemistry page."""
from __future__ import annotations

import html
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import seqprimers as sp
import seqwell_s3 as S
from chemdraw import oligo, protocol_graph_panels
from page import head, info, table
from protocol_graph import AnnotationSlot
from protocol_views import assemble_page

OUT = HERE.parent / "seq-well-s3.html"


def oligos() -> str:
    rows = [
        oligo("Seq-Well bead oligo", list(S.bead_oligo()), five="bead-|-5'-"),
        oligo("Template-switching oligo", S.tso_oligo()),
        oligo("S3 randomer", S.randomer()),
        oligo("SMART PCR primer", [S.seg("SMART handle", S.rt.SMART_HANDLE, "tso")]),
        oligo("P5-SMART hybrid", list(S.p5_primer()),
              three="-3' (terminal *A*C bonds)"),
        oligo("Custom Read 1 primer", [S.seg("spacer", S.P5_SPACER),
                                        S.seg("SMART handle", S.rt.SMART_HANDLE, "tso"),
                                        S.seg("AC", S.BEAD_SUFFIX)]),
        oligo("N7xx library-PCR primer", list(S.n7xx_primer())),
    ]
    return "<h2>Key oligos</h2><seq>" + "".join(rows) + "</seq>"


def construction(graph) -> str:
    return '<h2>Library construction</h2>\n'+protocol_graph_panels(graph,cls="long")


def sequencing() -> str:
    lib = S.final_library()
    rows = []
    for primer in S.SEQ_PRIMERS:
        hit = sp.locate(lib, primer)
        rows.append((html.escape(primer.role), html.escape(primer.name),
                     html.escape(", ".join(hit.covers)),
                     f"<code>{html.escape(hit.reads)}</code>&hellip;"))
    return f'''<h2>Sequencing</h2>
{info('Read 1 reports the 12-nt cell barcode and 8-nt UMI; Index 1 reports the N700 sample index; Read 2 reports cDNA from the Tn5 end.')}
{sp.diagram(lib, S.SEQ_PRIMERS)}
{table(("Read", "Primer", "Primer site", "First bases"), rows)}
{sp.read_content_table(lib,S.SEQ_PRIMERS,
    {"Read 1":20,"Index 1 (i7)":8,"Read 2":50},
    acquired_roles=tuple(p.role for p in S.SEQ_PRIMERS))}'''


def main() -> None:
    graph = S.protocol_graph()
    page = assemble_page(graph, [head("Seq-Well S3 library chemistry"), '<div class="wrap" data-protocol-complete="true">',
        '<h1>Seq-Well S3 &mdash; randomly primed second-strand capture</h1>',
        info("""Defining source: <a href="https://doi.org/10.1016/j.immuni.2020.09.015">Hughes et al., <i>Immunity</i> (2020)</a>. The <a href="https://shaleklab.com/wp-content/uploads/2024/02/Hughes.et_.al_.Master.Protocol.pdf">author Seq-Well S3 Master Protocol</a> (June 2019, pp. 18–24 and 39) prints the bead and custom Read 1 oligos and the Exonuclease I → NaOH → second-strand sequence. INFERRED — the Nextera right-arm bases and standard Index 1/Read 2 binding sites use official vendor family context; the defining experiment's kit revision remains unassigned. The displayed run is the author guide's MiSeq/NextSeq 20/8/50 option. The paper names an N700 index family but not a selected index member; the final library therefore shows eight variable index bases."""),
        AnnotationSlot.INTRO, oligos(), construction(graph), AnnotationSlot.AFTER_WORKFLOW,
        sequencing(), '</div>'])
    OUT.write_text(page, encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
