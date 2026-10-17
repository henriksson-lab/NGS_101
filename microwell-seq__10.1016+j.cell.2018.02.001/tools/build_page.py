#!/usr/bin/env python3
"""Build the diagram-centric Microwell-seq chemistry page."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import microwell_seq as M
import seqprimers as sp
from chemdraw import annotation_rows, oligo, panel, strand_row, protocol_graph_panels
from page import head, info

OUT = HERE.parent / "microwell-seq.html"


def oligos() -> str:
    rows = [
        oligo("Final seqA split-pool bead oligo", list(M.bead_oligo()), five="bead-|-5'-"),
        *(oligo(f"Split-pool {o.name.removeprefix('Table S1 ')} (Table S1)", list(o))
          for o in M.barcode_oligo_sets()),
        oligo("Template-switching oligo", M.tso_oligo()),
        oligo("TSO-PCR primer", [M.seg("SMART handle", M.rt.SMART_HANDLE, "tso")]),
        oligo("P5 bead-end-selective primer", [M.seg("P5", M.il.P5, "p5"),
                                                M.seg("spacer", M.P5_SPACER),
                                                M.seg("SMART handle", M.rt.SMART_HANDLE, "tso"),
                                                M.seg("AC", "AC")], three="-*A*C-3'"),
        oligo("N701 library-PCR primer", list(M.p7_primer(M.I7_OLIGO_INDEX))),
    ]
    return "<h2>Key oligos</h2><seq>" + "".join(rows) + "</seq>"


def construction() -> str:
    bead = M.bead_oligo()
    lib = M.final_library()
    return f'''<h2>Library construction</h2>
<h3>(1) Build the bead barcode by three split-pool extension rounds</h3>
{panel([strand_row(bead), *annotation_rows(bead)], cls="small",
       caption="The seqA bead carries three 6-nt barcode blocks, two fixed 15-nt linkers, a 6-nt UMI and T30.")}
<h3>(2) Capture microwell RNA, reverse-transcribe and template-switch</h3>
{panel(M.rt_scene().rows(), cls="small",
       caption="Collected beads are pooled for reverse transcription; the TSO installs the second SMART handle.")}
<h3>(3) Single-primer WTA</h3>
{panel([strand_row(bead), *annotation_rows(bead)], cls="small",
       caption="TSO-PCR amplifies SMART-handle-flanked cDNA before library construction.")}
<h3>(4) INFERRED &mdash; single-insertion-sequence tagmentation and bead-end-selective PCR</h3>
{panel([*M.final_scene().rows(), *annotation_rows(lib)], cls="small",
       caption="INFERRED — STAR Methods: a customized TruePrep V2 (Vazyme) transposase carrying two identical insertion sequences, with the kit N5xx primers replaced by the Table S2 P5 primer so only bead-end fragments amplify; other fragments self-loop. The insertion-sequence bases are not printed: dotted s7/ME and P7-side regions follow the N7xx primer.")}'''


def sequencing() -> str:
    lib = M.final_library()
    return sp.section(lib, M.SEQ_PRIMERS, heading="Sequencing",
        intro='Read 1 needs at least 54 cycles to cover the three cell-barcode blocks and UMI; Index 1 reports the 8-nt sample index; Read 2 reports cDNA. These cycle numbers are informative minima computed from the Table S2 oligos, not acquired run lengths: the 2018 STAR Methods say only that samples were sequenced \u201con an Illumina Hiseq system\u201d (GEO GSE108097 samples: HiSeq 2500), and print no cycle counts. The Guo lab\u2019s undated web protocol instead names NextSeq 500 with the custom Read 1 primer.',
        required_roles=M.RUN_ROLES,
        read_lengths={"Read 1":54,"Index 1 (i7)":8},
        read_uses={"Read 1":"extract"}, show_read_content=True)


def main() -> None:
    graph = M.protocol_graph()
    page = "\n".join([head("Microwell-seq library chemistry"), '<div class="wrap" data-protocol-complete="true">',
        '<h1>Microwell-seq &mdash; split-pool-barcoded bead capture</h1>',
        info('Defining source: <a href="https://doi.org/10.1016/j.cell.2018.02.001">Han et al., <i>Cell</i> (2018)</a>.'),
        '<div class="caveat"><b>Sources and library-adaptor boundary.</b> Reagents, conditions and step order come from the <i>Cell</i> STAR Methods (with the 2018 correction, which changes only the bead supplier and the lysis-buffer Tris concentration); oligos from the Key Resources Table and Supplementary Tables S1/S2. The STAR Methods use a customized TruePrep V2 (Vazyme) transposase \u201cthat carries two identical insertion sequences\u201d but do not print them, so the s7/ME end is inferred from the N7xx primer and drawn dotted. The Guo lab\'s undated web protocol at bis.zju.edu.cn/MCA describes a different implementation (Nextera XT, Exo I 45 min, 14 WTA cycles, 13\u201315 library cycles, NextSeq 500); the model follows the published Methods (Exo I 60 min, 4 + 10\u201314 WTA cycles, 14 library cycles, HiSeq). Table S2 also prints a seqB bead and Read 1 primer ending \u2026GAGTAC (no extra GT); the page draws seqA.</div>',
        oligos(), '<h2>Reaction workflow</h2>', protocol_graph_panels(graph, cls="long"),
        sequencing(), '</div>'])
    OUT.write_text(page, encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
