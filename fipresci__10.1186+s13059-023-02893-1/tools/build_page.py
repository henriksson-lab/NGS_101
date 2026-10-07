#!/usr/bin/env python3
"""Build the diagram-centric FIPRESCI chemistry page."""
from __future__ import annotations

import html
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import fipresci as F
import seqprimers as sp
from chemdraw import annotation_rows, oligo, panel, strand_row
from page import head, info, table

OUT = HERE.parent / "fipresci.html"


def oligos() -> str:
    rows = [
        oligo("Round-1 indexed Tn5 top", F.indexed_tn5_top()),
        oligo("Tn5 ME bottom", [F.seg("ME reverse complement", F.nx.ME_RC, "me")],
              mods="/5Phos/", three="-/3ddC/-3'"),
        oligo("10x 5-prime gel-bead TSO", F.bead_tso()),
        oligo("S5R-P5 enrichment primer", [F.seg("P5", F.il.P5, "p5"),
                                            F.seg("Read-1 tail", F.S5R_TAIL, "r1")]),
        oligo("S-P7 sample-index primer", [F.seg("P7", F.il.P7, "p7"),
                                            F.seg("i7", F.SAMPLE_INDEX_OLIGO, "cbc"),
                                            F.seg("TruSeq Read 2", F.il.TRUSEQ_READ2, "r2")],
              three="-s-T-3'"),
    ]
    return "<h2>Key oligos</h2><seq>" + "".join(rows) + "</seq>"


def construction() -> str:
    lib = F.final_library()
    return f'''<h2>Library construction</h2>
<h3>(1) Reverse-transcribe RNA inside fixed cells</h3>
{panel([strand_row([F.seg("first-strand cDNA", "X" * 28, placeholder=True),
                    F.seg("untemplated C", "CCC")])], cls="small",
       caption="No TSO is present yet; reverse transcriptase leaves the 3' C tail used in droplets.")}
<h3>(2) Tagment RNA/cDNA hybrids in 96 wells</h3>
{panel(F.hybrid_tagment_scene().rows(), cls="small",
       caption="A homodimeric Tn5 adaptor installs the 6-nt round-1 barcode and TruSeq Read-2 handle.")}
<h3>(3) Overload a 10x 5-prime channel and template-switch onto the bead TSO</h3>
{panel([strand_row(F.bead_tso()), *annotation_rows(F.bead_tso())], cls="small",
       caption="The 16-nt droplet barcode and 10-nt UMI form the second cell index on the transcript's 5' fragment.")}
<h3>(4) Enrich the TSO-bearing end, capture biotin products, then library-PCR</h3>
{panel([*F.final_scene().rows(), *annotation_rows(lib)], cls="small",
       caption="Only the 5'-terminal transcript fragment carries both the indexed Tn5 end and the droplet TSO end.")}'''


def sequencing() -> str:
    lib = F.final_library()
    rows = []
    for primer in F.SEQ_PRIMERS:
        hit = sp.locate(lib, primer)
        rows.append((html.escape(primer.role), html.escape(primer.name),
                     html.escape(", ".join(hit.covers)),
                     f"<code>{html.escape(hit.reads)}</code>&hellip;"))
    return f'''<h2>Sequencing</h2>
{info('Paired-end 150: Read 1 reports the 16-nt droplet barcode, 10-nt UMI and transcript 5-prime end; Index 1 reports the sample index; Read 2 starts with the 6-nt round-1 barcode.')}
{table(("Read", "Primer", "Primer site", "First bases"), rows)}'''


def main() -> None:
    page = "\n".join([head("FIPRESCI library chemistry"), '<div class="wrap">',
        '<h1>FIPRESCI &mdash; plate-indexed hybrid tagmentation plus 10x 5-prime barcoding</h1>',
        info('Defining source: <a href="https://doi.org/10.1186/s13059-023-02893-1">Luo et al., <i>Genome Biology</i> (2023)</a>.'),
        oligos(), construction(), sequencing(), '</div>'])
    OUT.write_text(page, encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
