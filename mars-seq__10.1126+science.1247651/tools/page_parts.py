"""Sparse individual page renderers for MARS-seq and MARS-seq2.0."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import marsseq as M
import seqprimers as sp
from chemdraw import (Construct, Segment, oligo, panel, protocol_graph_panels, revcomp,
                      strand_row)
from page import head, info, table


def seg(name, top, tag=None, **kw):
    return Segment(name=name, top=top, tag=tag, **kw)


from mars_graph import PAPERS, graph_for as protocol_graph


def render_page(protocol: str, *, graph=None) -> str:
    graph = graph if graph is not None else protocol_graph(protocol)
    doi, journal, year = PAPERS[protocol]
    author = "Jaitin et al." if protocol == "MARS-seq" else "Keren-Shaul et al."
    exact = protocol == "MARS-seq2.0"
    caveat = ('The <a href="https://www.researchgate.net/publication/279842956_Jaitin-SM">'
              'coauthor-uploaded original Supplementary Methods</a> gives the reaction '
              "order: pool barcoded first-strand cDNA, remove free primers with Exo I, "
              "then synthesize the second strand and perform IVT. The separate PMC "
              "figure PDF remains unavailable. The library below uses the printed "
              "Supplementary Table S7 oligos."
              if not exact else
              "The supplementary tables publish RT1 and pool-barcoded ligation adapters. "
              "The authors' detailed protocols.io procedure publishes RT2 and both "
              "enrichment primers (reagent setup, p. 29); the library flanks are assembled "
              "from those exact oligos. Separate instrument sequencing-primer oligos are "
              "not specified.")
    oligos = [oligo("RT1", M.rt1_segments(protocol))]
    if exact:
        oligos.append(oligo("representative pool ligation adapter",
                            list(M.mars2_ligation_adapter()), mods="/5Phos/ … /3SpC3/"))
        oligos.extend([
            oligo("RT2", [seg("second RT primer", M.MARS2_RT2, "r1")]),
            oligo("P5-Read1 enrichment primer",
                  [seg("P5 + Read 1", M.MARS2_P5_PCR, "r1")]),
            oligo("P7-Read2 enrichment primer",
                  [seg("P7 + Read 2", M.MARS2_P7_PCR, "r2")]),
        ])
    else:
        oligos.extend([
            oligo("representative pool ligation adapter",
                  [seg("pool barcode", M.POOL_ORDERED, "cbc"),
                   seg("random diversity", "N" * M.MARS1_POOL_RANDOM_NT),
                   seg("ligation/RT2 handle", M.LIG_CONSTANT, "r1")],
                  mods="/5Phos/ … /3SpC3/"),
            oligo("P5_Rd1 PCR", [seg("P5 + Read 1", M.MARS1_P5_PCR, "r1")]),
            oligo("P7_Rd2 PCR", [seg("P7 + Read 2", M.MARS1_P7_PCR, "r2")]),
        ])
    final = M.final_library(protocol)
    plan = graph.sequencing.get("library")
    lengths = dict(plan.read_lengths) if plan is not None else {}
    if exact:
        r1 = Construct([seg("random diversity", "N" * 5, placeholder=True),
                        seg("pool barcode", M.POOL_READ, "cbc"),
                        seg("cDNA", "XXXXXXXX...", placeholder=True)])
        # Author procedure Step 80 acquires rd2 = 15 cycles: well barcode + UMI only.
        r2 = Construct([seg("well barcode", "N" * 7, "cbc", placeholder=True),
                        seg("UMI", "N" * 8, "umi", placeholder=True)])
        # The Supplementary Information's read-design fields print the worked example
        # R1_design = 3I.4P.50M with R2_design = 7W.8R.5I.  Read 1's 3I matches the
        # MARS-seq1 N3 adapter, so 5I is what the 2.0 N5 adapter implies -- derived here,
        # never published, and labelled as such.
        read_caption = ("Published-source difference: MARS-seq2.0 Supplementary Information "
                        "(MOESM1, seq_batches.txt field description) prints 3I.4P.50M "
                        "and 7W.8R.5I, whereas Supplementary Table 2 (MOESM5, ligation "
                        "adapter rows) prints N5 before the 4-nt pool barcode. "
                        "DERIVED — this model uses the table's N5 adapter, so Read 1 "
                        "starts 5I.4P; whether 3I is a v1 carryover remains uncertain. "
                        "Acquired cycles follow the author procedure (protocols.io, Step 80: "
                        "rd1 = 75 bp, rd2 = 15 bp, no indices); the MOESM1 example's "
                        "57 + 20 cycles describe a different configuration, so the 5I "
                        "poly(T) cycles of 7W.8R.5I are not drawn.")
        read_rows = [("Read 1", "N5, 4-nt pool barcode, sense cDNA"),
                     ("Read 2", "7-nt well barcode, 8-nt UMI")]
    else:
        r1 = Construct([seg("random diversity", "N" * 3, placeholder=True),
                        seg("pool barcode", M.POOL_READ, "cbc"),
                        seg("cDNA", "XXXXXXXX...", placeholder=True)])
        r2 = Construct([seg("well barcode", "N" * 6, "cbc", placeholder=True),
                        seg("UMI", "N" * 4, "umi", placeholder=True),
                        seg("poly(T), ignored", "TTTTT")])
        read_caption = "Published MARS-seq read design: N3, 4-nt pool barcode, then cDNA; Read 2 reports cell barcode and UMI."
        read_rows = [("Read 1", "N3, 4-nt pool barcode, sense cDNA"),
                     ("Read 2", "6-nt well barcode, 4-nt UMI, poly(T)")]
    return "\n".join([
        head(f"{protocol} library chemistry"),
        '<div class="wrap" data-protocol-complete="true">',
        f'<h1>{protocol} &mdash; three-level barcoded RNA-seq with T7 amplification</h1>',
        info(f'Defining source: <a href="https://doi.org/{doi}">{author}, '
             f'<i>{journal}</i> ({year})</a>.'),
        (info('Author reagent setup: <a href="https://doi.org/10.17504/protocols.io.7hkhj4w">'
              'Keren-Shaul et al., protocols.io (2019), p. 29</a>.') if exact else ''),
        f'<div class="caveat"><b>Sequence boundary.</b> {caveat}</div>',
        (info('Original-source sequence conflict: the '
              '<a href="https://www.researchgate.net/publication/279842956_Jaitin-SM">'
              'coauthor-uploaded Supplementary Methods</a> (lines 370–375) call RT2 '
              'complementary to the ligated adapter, '
              'while Supplementary Table S7 prints RT2 '
              f'<code>{M.MARS1_RT2_PRINTED}</code>. The reverse complement of the '
              f'printed adapter constant <code>{M.LIG_CONSTANT}</code> is '
              f'<code>{revcomp(M.LIG_CONSTANT)}</code>. Read right to left, the printed '
              f'RT2 is <code>{M.MARS1_RT2_PRINTED_REVERSED}</code>'
              + (', exactly that complement: the table entry appears to be written '
                 '3′→5′.' if M.MARS1_RT2_REVERSAL_IS_ADAPTER_COMPLEMENT else '.')
              + ' This is an apparent published typographical discrepancy, not a '
              'verified correction. The model keeps the printed RT2 as the source '
              'value, marks the RT2 step INFERRED and does not choose a landing; '
              'whether the ordered primer was the reversed string remains uncertain. '
              'The later MARS-seq2.0 author protocol prints the complement 5′→3′ as '
              'its own RT2, which is a separate source and does not settle the 2014 '
              'order.') if not exact else ''),
        f'<h2>Key oligos</h2><seq>{"".join(oligos)}</seq>',
        '<h2>Library construction</h2>', protocol_graph_panels(graph),
        sp.diagram(final, M.SEQ_PRIMERS),
        '<h2>Read layout</h2>',
        panel([strand_row(r1, prefix="Read 1  ", suffix=""),
               strand_row(r2, prefix="Read 2  ", suffix="")], cls="small",
              caption=read_caption),
        table(("Read", "Reports"), read_rows),
        sp.read_spec_table([
            sp.read_spec("Read 1", lengths.get("Read 1"), read_rows[0][1],
                         "extract_then_include"),
            sp.read_spec("Read 2", lengths.get("Read 2"), read_rows[1][1], "extract"),
        ], ("Read contents are source-declared; cycles are the author procedure's Step 80 "
            "run (NextSeq 500 High Output v2, no index reads)." if exact else
            "Read contents are source-declared; the defining sources do not fix a complete run recipe.")),
        '</div>',
    ])
