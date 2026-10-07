#!/usr/bin/env python3
"""Build the diagram-centric snATAC-seq page."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import snatac as S
from chemdraw import Scene, Segment, annotation_rows, complement_segments, oligo, panel
import nextera as nx
from page import head, info, table
import seqprimers as sp

OUT = HERE.parent / "snatac-seq.html"


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def preamble() -> str:
    return f'''<div class="wrap">
<h1>snATAC-seq &mdash; two-level combinatorial indexing of accessible chromatin</h1>
{info('Defining source: <a href="https://doi.org/10.1038/s41593-018-0079-3">Preissl et al., 2018</a>.')}
{info('The first barcode pair is carried by the two Tn5 adapters; a second i5/i7 pair is added in PCR. One nucleus is identified by all four 8-nt barcodes. There is no UMI.')}
'''


def oligos() -> str:
    rows = [
        oligo("A / p5 transposon (8 variants)", S.p5_transposon()),
        oligo("B / p7 transposon (12 variants)", S.p7_transposon()),
        oligo("pMENTS", [seg("ME'", S.PMENTS, "me")], mods="/5Phos/"),
        oligo("Custom Read 1 primer", [seg("site", S.READ1_SITE, "r1"), seg("ME", nx.ME, "me")]),
        oligo("Custom Read 2 primer", [seg("site", S.READ2_SITE, "r1"), seg("ME", nx.ME, "me")]),
        oligo("Custom Index 1 primer", [seg("ME' + Read 2 site'", S.INDEX1_PRIMER, "r1")]),
        oligo("i5 PCR primer", S.p5_index_primer()),
        oligo("i7 PCR primer", S.p7_index_primer()),
    ]
    return "<h2>Key oligos</h2>\n<seq>\n" + "\n".join(rows) + "\n</seq>"


def loaded_transposon(parts: list[Segment], name: str) -> Scene:
    sc = Scene()
    sc.strand(name, parts, label=name)
    sc.anneal("pMENTS", [seg("ME'", S.PMENTS, "me")], to=name,
              pair=("ME'", "ME"), label="pMENTS", mod5="p")
    sc.mark(name, "p5 barcode" if name == "A / p5" else "p7 barcode",
            "round-1 barcode")
    return sc


def tagmented_scene() -> Scene:
    """Representative amplifiable A...B fragment before 9-nt gap fill."""
    g = nx.TAGMENTATION_GAP
    left = S.p5_transposon()
    right = S.p7_transposon()
    top = [*left, seg("top gap-facing", "X" * g, placeholder=True),
           seg("genomic core", "XXXXXXXX...XXXXXXXX", placeholder=True)]
    bottom = [*right, seg("bottom gap-facing", "x" * g, placeholder=True),
              *complement_segments([top[-1]])]
    sc = Scene()
    sc.strand("A end", top, label="A end")
    sc.anneal("B end", bottom, to="A end", pair=("genomic core'", "genomic core"), label="B end")
    sc.anneal("A pMENTS", [seg("ME", S.PMENTS, "me")], to="A end",
              pair=("ME", "ME"), label="")
    sc.anneal("B pMENTS", [seg("ME", S.PMENTS, "me")], to="B end",
              pair=("ME", "ME"), label="", above=True)
    sc.mark("A end", "top gap-facing", "9-nt gap opposite this end")
    sc.mark("B end", "bottom gap-facing", "9-nt gap opposite this end")
    return sc


def pcr_scene() -> Scene:
    tagged = S.tagged_fragment()
    sc = Scene.duplex(list(tagged), label="")
    sc.anneal("i5 primer", S.p5_index_primer(), to="bottom", pair=("s5", "s5'"))
    sc.arrow("i5 primer", "")
    sc.anneal("i7 primer", S.p7_index_primer(), to="top", pair=("s7", "s7'"), above=True)
    sc.arrow("i7 primer", "")
    return sc


def final() -> str:
    lib = S.final_library()
    duplex = Scene.duplex(lib.segments)
    duplex.strands["top"].label = duplex.strands["bottom"].label = ""
    rows = [*duplex.rows(), *annotation_rows(lib)]
    return panel(rows, caption="The four variable 8-nt blocks identify the nucleus: p5 + p7 record the tagmentation well; i5 + i7 record the PCR well.")


def steps() -> str:
    a = loaded_transposon(S.p5_transposon(), "A / p5")
    b = loaded_transposon(S.p7_transposon(), "B / p7")
    return f'''<h2>Library generation</h2>
<h3>(1) Assemble two indexed transposomes per tagmentation well</h3>
{panel(a.rows(), cls="small", caption="One of eight A / p5 adapters annealed to the shared phosphorylated pMENTS strand.")}
{panel(b.rows(), cls="small", caption="One of twelve B / p7 adapters annealed to pMENTS. Pairing one A and one B gives 96 tagmentation-well combinations.")}

<h3>(2) Tagment accessible genomic DNA</h3>
{panel(tagmented_scene().rows(), caption="Representative A...B product before gap fill. A...A and B...B molecules lack one PCR-primer site and do not amplify exponentially.")}

<h3>(3) Fill the 9-nt gaps and add the second barcode pair by PCR</h3>
{panel(pcr_scene().rows(), caption="The i5 primer lands on s5; the i7 primer lands on s7. Their P5/P7 and 8-nt index tails are copied into the product.")}

<h3>(4) Final library</h3>
{final()}
'''


def read_scene(primer, label: str) -> str:
    lib = S.final_library()
    hit = sp.locate(lib, primer)
    if hit is None:
        raise ValueError(f"{primer.name} has no site")
    sc = Scene.duplex(list(lib), label="")
    if primer.role == "Read 1":
        ps = [seg("Read 1 site", S.READ1_SITE, "r1"), seg("ME", nx.ME, "me")]
        sc.anneal("primer", ps, to="bottom", pair=("Read 1 site", "Read 1 site'"))
    elif primer.role == "Read 2":
        ps = [seg("Read 2 site", S.READ2_SITE, "r1"), seg("ME", nx.ME, "me")]
        sc.anneal("primer", ps, to="top", pair=("Read 2 site", "Read 2 site'"), above=True)
    elif primer.role == "Index 1 (i7)":
        ps = [seg("Index 1", S.INDEX1_PRIMER, "r1")]
        sc.anneal("primer", ps, to="bottom", pair=("Index 1", "ME''"))
    else:
        ps = [seg("P5", primer.seq, "p5")]
        sc.anneal("primer", ps, to="bottom", pair=("P5", "P5'"))
    sc.arrow("primer", label)
    return panel(sc.rows())


def sequencing() -> str:
    lib = S.final_library()
    # Locate every primer again at the consumer boundary; malformed read layouts fail the build.
    landings = {p.role: sp.locate(lib, p) for p in S.SEQ_PRIMERS}
    if any(v is None for v in landings.values()):
        raise ValueError("incomplete snATAC read layout")
    rows = [
        ("Read 1", "50", "genomic insert from the A end"),
        ("Index 1", "43", "p7 barcode (8) &middot; B linker (27) &middot; i7 (8)"),
        ("Index 2", "37", "i5 (8) &middot; A linker (21) &middot; p5 barcode (8)"),
        ("Read 2", "50", "genomic insert from the B end"),
    ]
    return f'''<h2>Read layout</h2>
{table(("Read", "Cycles", "Content"), rows)}
<h3>Read 1 &mdash; insert from the A end</h3>
{read_scene(S.SEQ_PRIMERS[0], "Read 1: genomic insert")}
<h3>Index 1 &mdash; p7 barcode, linker, i7</h3>
{read_scene(S.SEQ_PRIMERS[1], "43 cycles")}
<h3>Index 2 &mdash; i5, linker, p5 barcode</h3>
{read_scene(S.SEQ_PRIMERS[2], "37 cycles")}
<h3>Read 2 &mdash; insert from the B end</h3>
{read_scene(S.SEQ_PRIMERS[3], "Read 2: genomic insert")}
</div>'''


def main() -> None:
    OUT.write_text("\n".join([head("snATAC-seq library chemistry"), preamble(), oligos(),
                              steps(), sequencing()]), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
