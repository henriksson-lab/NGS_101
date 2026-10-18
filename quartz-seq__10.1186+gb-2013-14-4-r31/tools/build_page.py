#!/usr/bin/env python3
"""Build the diagram-centric Quartz-Seq chemistry page."""
from __future__ import annotations

import html
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import quartz_seq as Q
import seqprimers as sp
from chemdraw import Construct, Segment, oligo, protocol_graph_panels
from page import head, info, table


OUT = HERE.parent / "quartz-seq.html"


def seg(name, top, tag=None, **kw):
    return Segment(name=name, top=top, tag=tag, **kw)


def oligos() -> str:
    rows = [
        oligo("RT primer", Q.rt_primer()),
        oligo("Tagging primer", Q.tagging_primer()),
        oligo("Suppression primer", [seg("G", "G"), seg("M handle", Q.M, "tso")],
              mods="/5Amine/"),
        oligo("TRSU adapter strand", [seg("P5 + Read 1", Q.TRSU, "p5+r1")],
              three="-3' (terminal phosphorothioate)"),
        oligo("TRSI-2 adapter strand", [seg("Index 1 primer", Q.TRSI2[:33], "r2"),
                                         seg("i7", Q.TRSI2_INDEX, "cbc"),
                                         seg("P7 reverse complement", Q.TRSI2[-24:], "p7")],
              mods="/5Phos/", three="-3' (terminal phosphorothioate)"),
    ]
    return "<h2>Key oligos</h2>\n<seq>\n" + "\n".join(rows) + "\n</seq>"



def sequencing() -> str:
    lib = Q.final_library()
    rows = []
    for primer, hit in Q.primer_landings():
        rows.append((html.escape(primer.role), html.escape(primer.name),
                     html.escape(", ".join(hit.covers)),
                     f"<code>{html.escape(hit.reads)}</code>&hellip;"))
    return f'''<h2>Sequencing</h2>
{info('Paired-end 50-base reads with a source 6-base index sequence; acquired index-cycle count remains unassigned in this capped pass. Cell identity is the library tube/index; there is no molecular cell barcode or UMI.')}
{sp.diagram(lib, Q.SEQ_PRIMERS)}
{table(("Read", "Primer", "Primer site", "First bases"), rows)}
{sp.read_content_table(lib, Q.SEQ_PRIMERS,
    {"Read 1":50,"Index 1 (i7)":None,"Read 2":50},
    acquired_roles=tuple(p.role for p in Q.SEQ_PRIMERS))}
'''


def main() -> None:
    graph = Q.protocol_graph()
    construction = "<h2>Library construction</h2>\n" + protocol_graph_panels(graph)
    page = "\n".join([head("Quartz-Seq library chemistry"),
        '<div class="wrap" data-protocol-complete="true">',
        '<h1>Quartz-Seq &mdash; poly(A)-tagging whole-transcript amplification</h1>',
        info('Defining source: <a href="https://doi.org/10.1186/gb-2013-14-4-r31">Sasagawa et al., <i>Genome Biology</i> (2013)</a>.'),
        info('RNase H and TdT act in the same source tailing mixture. Native second-strand and adapter-ligated precursors enter the two PCR stages; residual original polymers remain explicitly bounded.'), oligos(), construction, sequencing(), '</div>'])
    OUT.write_text(page, encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
