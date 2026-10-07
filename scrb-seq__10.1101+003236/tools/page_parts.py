#!/usr/bin/env python3
"""Shared drawing implementation for individual SCRB-seq protocol pages."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import scrb as S
from chemdraw import Scene, Segment, annotation_rows, complement_segments, oligo, panel, strand_row
import illumina as il
import nextera as nx
from page import head, info, table
import seqprimers as sp

PAPERS = {
    S.SCRB: ("10.1101/003236", "SCRB-seq"),
    S.MCSCRB: ("10.1038/s41467-018-05347-6", "mcSCRB-seq"),
}


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def preamble(protocol: str) -> str:
    doi, title = PAPERS[protocol]
    distinguishing = ("The TSO has a 5' iso-base block; cells are barcoded in RT and pooled before amplification."
                      if protocol == S.SCRB else
                      "Molecular crowding (7.5% PEG 8000) enhances RT; the TSO is unblocked. Cells are barcoded in RT and pooled before amplification.")
    return f'''<div class="wrap">
<h1>{title} &mdash; plate-based 3' tag single-cell RNA sequencing</h1>
{info(f'Defining source: <a href="https://doi.org/{doi}">doi:{doi}</a>.')}
{info(distinguishing + ' Each molecule carries a 6-nt well barcode and 10-nt UMI.')}
'''


def oligos(protocol: str) -> str:
    rows = [
        oligo("Barcoded oligo-dT (E3V6NEXT)", S.oligo_dt_segments(), mods="5'-biotin"),
        oligo("Template-switching oligo (E5V6NEXT)", S.tso_segments(protocol)),
        oligo("PreAmp primer (SINGV6)", [seg("PCR handle", S.PCR_HANDLE, "r1")], mods="5'-biotin"),
        oligo("3' enrichment primer (P5NEXTPT5)",
              [seg("P5", il.P5, "p5"), seg("Read 1 remainder", S.FULL_HANDLE[4:], "r1")],
              mods="four 3'-terminal phosphorothioate bonds"),
        oligo("N7xx", [seg("P7", il.P7, "p7"), seg("i7", "I" * 8, "cbc", placeholder=True), seg("s7", nx.S7, "s7")]),
    ]
    return "<h2>Key oligos</h2>\n<seq>\n" + "\n".join(rows) + "\n</seq>"


def rt_scenes(protocol: str) -> tuple[Scene, Scene]:
    mrna = [seg("body", "XXXXXXXX...XXXXXXXX", placeholder=True), seg("poly(A)", "A" * S.DT_LEN)]
    primer = [seg("R1 handle", S.FULL_HANDLE, "r1"),
              seg("cell barcode", "B" * S.BARCODE_LEN, "cbc", placeholder=True),
              seg("UMI", "N" * S.UMI_LEN, "umi", placeholder=True),
              seg("dT", "T" * S.DT_LEN)]
    first = [*primer, *complement_segments([mrna[0]]), seg("CCC", "CCC", "tso")]
    tso = [seg("PCR handle", S.PCR_HANDLE, "r1"), seg("GGG", "GGG", "tso")]
    p = Scene(); p.strand("mRNA", mrna); p.anneal("barcoded oligo-dT", primer, to="mRNA", pair=("dT", "poly(A)"), mod5="biotin"); p.arrow("barcoded oligo-dT", "reverse transcriptase")
    sw = Scene(); sw.strand("mRNA", mrna); sw.anneal("first-strand cDNA", first, to="mRNA", pair=("dT", "poly(A)"), mod5="biotin"); sw.anneal("TSO", tso, to="first-strand cDNA", pair=("GGG", "CCC"), above=True); sw.arrow("first-strand cDNA", "template switch")
    return p, sw


def preamp_scene() -> Scene:
    top = [seg("PCR handle, TSO end", S.PCR_HANDLE, "r1"),
           seg("cDNA", "XXXXXXXX...XXXXXXXX", placeholder=True),
           seg("PCR handle, oligo-dT end", S.PCR_HANDLE, "r1")]
    sc = Scene.duplex(top)
    sc.anneal("PreAmp", [seg("PCR handle", S.PCR_HANDLE, "r1")], to="bottom",
              pair=("PCR handle", "PCR handle, TSO end'"))
    sc.arrow("PreAmp", "")
    return sc


def enrichment_scene() -> Scene:
    top = [seg("R1 overlap", S.FULL_HANDLE[:4], "r1"),
           seg("R1 anneal", S.FULL_HANDLE[4:], "r1"),
           seg("cell barcode", "B" * S.BARCODE_LEN, "cbc", placeholder=True),
           seg("UMI", "N" * S.UMI_LEN, "umi", placeholder=True),
           seg("poly(T)", "T" * S.DT_LEN), seg("VN", "VN", placeholder=True),
           seg("cDNA", "XXXXXXXX...XXXXXXXX", placeholder=True),
           seg("ME'", nx.ME_RC, "me"), seg("s7'", nx.S7_RC, "s7")]
    sc = Scene.duplex(top)
    p5 = [seg("P5 tail", il.P5, "p5"), seg("R1 anneal", S.FULL_HANDLE[4:], "r1")]
    sc.anneal("P5NEXTPT5", p5, to="bottom", pair=("R1 anneal", "R1 anneal'"),
              unpaired=("P5 tail",))
    sc.arrow("P5NEXTPT5", "")
    n7 = [seg("P7 tail", il.P7, "p7"), seg("i7", "I" * 8, "cbc", placeholder=True),
          seg("s7", nx.S7, "s7")]
    sc.anneal("N7xx", n7, to="top", pair=("s7", "s7'"), above=True,
              unpaired=("P7 tail", "i7"))
    sc.arrow("N7xx", "")
    return sc


def final_panel(protocol: str) -> str:
    lib = S.final_library(protocol)
    return panel([strand_row(lib, "top"), strand_row(lib, "bottom"), *annotation_rows(lib)],
                 caption="Only the oligo-dT / s7 tagmentation product receives both P5 and P7 and amplifies exponentially.")


def steps(protocol: str) -> str:
    prime, switch = rt_scenes(protocol)
    rt_caption = ("The blocked TSO contributes its 22-nt PCR handle; the 5' iso-bases are not copied."
                  if protocol == S.SCRB else
                  "The unblocked TSO contributes the same 22-nt PCR handle; PEG crowding is present during this reaction.")
    return f'''<h2>Library generation</h2>
<h3>(1) Barcoded reverse transcription</h3>
{panel(prime.rows(), cls="small", caption="The oligo-dT installs the well barcode and UMI before the poly(T) tract.")}
<h3>(2) Template switching</h3>
{panel(switch.rows(), cls="small", caption=rt_caption)}
<h3>(3) Pool wells, remove unused primer and amplify full-length cDNA</h3>
{panel(preamp_scene().rows(), cls="small", caption="SINGV6 amplifies from the shared 22-nt handle at both ends. The well barcode makes pooling irreversible but traceable.")}
<h3>(4) Tagment and selectively amplify the 3' cDNA end</h3>
{panel(enrichment_scene().rows(), caption="P5NEXTPT5 requires the full oligo-dT-side TruSeq Read-1 handle; N7xx lands on s7. TSO-end and s5 products lack one usable primer site.")}
<h3>(5) Final library</h3>
{final_panel(protocol)}
'''


def sequencing(protocol: str) -> str:
    lib = S.final_library(protocol)
    hits = {p.role: sp.locate(lib, p) for p in S.READ_PRIMERS}
    if any(h is None for h in hits.values()):
        raise ValueError("incomplete SCRB read layout")
    if protocol == S.SCRB:
        rows = [("Read 1", "17", "well barcode 1&ndash;6 &middot; UMI 7&ndash;16 &middot; first poly(T) base"),
                ("Index 1", "8", "plate i7"), ("Read 2", "34", "cDNA, sense to the mRNA")]
    else:
        rows = [("Read 1", "16", "well barcode 1&ndash;6 &middot; UMI 7&ndash;16"),
                ("Index 1", "8", "plate i7"), ("Read 2", "50", "cDNA, sense to the mRNA")]
    return '<h2>Read layout</h2>' + table(("Read", "Cycles", "Content"), rows) + '</div>'


def render_page(protocol: str, out: Path) -> None:
    out.write_text("\n".join([head(f"{protocol} library chemistry"), preamble(protocol),
                              oligos(protocol), steps(protocol), sequencing(protocol)]),
                   encoding="utf-8")
    print(f"wrote {out}  ({out.stat().st_size:,} bytes)")

