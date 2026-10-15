#!/usr/bin/env python3
"""Shared drawing implementation for individual SCRB-seq protocol pages."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import scrb as S
from chemdraw import (Construct, Scene, Segment, complement_segments, oligo,
                      protocol_graph_panels)
import illumina as il
import nextera as nx
from page import head, info
import seqprimers as sp
from scrb_graph import PAPERS, graph_for as protocol_graph



def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def preamble(protocol: str, *, complete: bool = False) -> str:
    doi, title = PAPERS[protocol]
    notes = ("01_scrb-seq.html" if protocol == S.SCRB else
             "../scrb-seq__10.1101+003236/01_scrb-seq.html")
    distinguishing = ("The TSO has a 5' iso-base block; cells are barcoded in RT and pooled before amplification."
                      if protocol == S.SCRB else
                      "Molecular crowding (7.5% PEG 8000) enhances RT; the TSO is unblocked. Cells are barcoded in RT and pooled before amplification.")
    complete_attr = ' data-protocol-complete="true"' if complete else ""
    return f'''<div class="wrap"{complete_attr}>
<h1>{title} &mdash; plate-based 3' tag single-cell RNA sequencing</h1>
<p class="research-notes"><a href="{notes}">Research notes</a></p>
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




def preamp_scene() -> Scene:
    top = list(S.preamp_product())
    sc = Scene.duplex(top)
    sc.anneal("PreAmp", [seg("PCR handle", S.PCR_HANDLE, "r1")], to="bottom",
              pair=("PCR handle", "oligo-dT-end PCR handle'"))
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




def sequencing(protocol: str) -> str:
    lib = S.final_library(protocol)
    if protocol == S.SCRB:
        lengths = {"Read 1": 17, "Index 1 (i7)": 8, "Read 2": 34}
        intro = ("Read 1 continues one base into poly(T) after its identifiers; Read 2 "
                 "reports 34 bases of cDNA in the mRNA sense.")
    else:
        lengths = {"Read 1": 16, "Index 1 (i7)": 8, "Read 2": 50}
        intro = "Read 1 ends with the UMI; Read 2 reports 50 bases of cDNA in the mRNA sense."
    return (sp.section(lib, S.READ_PRIMERS, heading="Read layout", intro=intro,
                       required_roles=tuple(p.role for p in S.READ_PRIMERS),
                       read_lengths=lengths) + '</div>')


def render_page(protocol: str, out: Path) -> None:
    graph = protocol_graph(protocol)
    workflow = "<h2>Library generation</h2>\n" + protocol_graph_panels(graph)
    out.write_text("\n".join([head(f"{protocol} library chemistry"),
                              preamble(protocol, complete=True), oligos(protocol),
                              workflow, sequencing(protocol)]),
                   encoding="utf-8")
    print(f"wrote {out}  ({out.stat().st_size:,} bytes)")
