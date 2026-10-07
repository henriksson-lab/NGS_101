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
from chemdraw import Scene, Segment, annotation_rows, oligo, panel, strand_row
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


def construction() -> str:
    tailed = Q.tailed_first_strand()
    wta = Q.wta_amplicon()
    lib = Q.final_library()
    return f'''<h2>Library construction</h2>
<h3>(1) Reverse-transcribe from the minimal oligo-dT primer</h3>
{panel(Q.mrna_rt_scene().rows(), cls="small",
       caption="The RT primer contributes the M amplification handle; exonuclease I removes unused primer after reverse transcription.")}

<h3>(2) Add a short poly(A) tail to first-strand cDNA</h3>
{panel([strand_row(tailed), *annotation_rows(tailed)], cls="small",
       caption="Terminal transferase creates the second priming site.")}

<h3>(3) Prime the new tail, synthesize the second strand, and suppress-amplify</h3>
{panel(Q.tagging_scene().rows(), cls="small",
       caption="The tagging primer contributes M at the opposite end.")}
{panel([*Scene.duplex(list(wta), label="WTA").rows(), *annotation_rows(wta)], cls="small",
       caption="One amplified full-length cDNA. A single G+M primer amplifies both ends; short M/M products self-anneal and are suppressed.")}

<h3>(4) Shear cDNA and ligate the single-index TruSeq adapter</h3>
{panel([*Scene.duplex(list(lib), label="library").rows(), *annotation_rows(lib)], cls="small",
       caption="Final Quartz-Seq library, P5 to P7'. The insert is a sheared fragment of amplified cDNA.")}
'''


def sequencing() -> str:
    rows = []
    for primer, hit in Q.primer_landings():
        rows.append((html.escape(primer.role), html.escape(primer.name),
                     html.escape(", ".join(hit.covers)),
                     f"<code>{html.escape(hit.reads)}</code>&hellip;"))
    return f'''<h2>Sequencing</h2>
{info('Paired-end 50-base reads with a 6-base Index 1. Cell identity is the library tube/index; there is no molecular cell barcode or UMI.')}
{table(("Read", "Primer", "Primer site", "First bases"), rows)}
'''


def main() -> None:
    page = "\n".join([head("Quartz-Seq library chemistry"), '<div class="wrap">',
        '<h1>Quartz-Seq &mdash; poly(A)-tagging whole-transcript amplification</h1>',
        info('Defining source: <a href="https://doi.org/10.1186/gb-2013-14-4-r31">Sasagawa et al., <i>Genome Biology</i> (2013)</a>.'),
        oligos(), construction(), sequencing(), '</div>'])
    OUT.write_text(page, encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
