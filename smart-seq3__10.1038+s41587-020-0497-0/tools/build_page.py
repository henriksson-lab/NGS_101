#!/usr/bin/env python3
"""Build the SMART-seq3 protocol schematic.

Molecular states (oligos, template-switched cDNA, Tn5 states, final libraries)
come from the shared SMART-seq drawing helpers. Reagents, conditions, cleanups
and step order are declared here from the SMART-seq3 author protocol V.3
(protocols.io DOI 10.17504/protocols.io.bcq4ivyw, 19-20 Feb 2020) and
Hagemann-Jensen et al. 2020 Supplementary Table 1; the published Methods are
behind a subscription preview and were not read.
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SHARED = HERE.parents[1] / "smart-seq__10.1038+nbt.2282" / "tools"
sys.path.insert(0, str(SHARED))

import page_parts as P  # noqa: E402
from chemdraw import Construct, Segment  # noqa: E402
from pcr import PCRGeometry, PCRPrimer  # noqa: E402
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph,  # noqa: E402
                            SequencingPlan, StepKind)

OUT = HERE.parent / "smart-seq3.html"
PROTOCOL = "SMART-seq3"
SOURCE = ("Defining source: doi:10.1038/s41587-020-0497-0 (Supplementary Table 1); "
          "author protocol V.3 doi:10.17504/protocols.io.bcq4ivyw")


def _cleanup(graph, identifier, label, materials, outputs, reason):
    """PEG/Sera-Mag bead cleanup: same molecules recovered, supernatant discarded."""
    beads = graph.input(identifier + "_beads",
                        "0.6:1 22% PEG / Sera-Mag Speed Bead suspension", MaterialKind.SELECTION)
    waste = graph.product(identifier + "_supernatant", "bead supernatant and ethanol washes",
                          MaterialKind.WASTE)
    graph.react(identifier, label, StepKind.CLEANUP,
                (*((m, InputRole.TEMPLATE) for m in materials), (beads, InputRole.SELECTION)),
                (*outputs, waste), changes_molecule=False)
    graph.discard(waste, reason)
    return outputs


def protocol_graph() -> ProtocolGraph:
    graph = ProtocolGraph(PROTOCOL, authority=SOURCE)
    mrna = P.input_mrna(graph)
    priming, first_scene, switch_scene = P.ss3_rt(PROTOCOL)
    template = P.ss3_preamp_template(PROTOCOL)

    # Author protocol V.3 steps 1-5 (PDF pp5-8): 3 uL lysis, 72 C 10 min; 1 uL RT mix.
    primer = graph.input("rt_primer", "Biotinylated Smartseq3_OligodT30VN (0.5 uM in RT)",
                         MaterialKind.OLIGO, state=priming)
    lysis = graph.input("lysis", "Lysis: 0.1% Triton X-100, 5% PEG 8000, 0.5 mM dNTP each, "
                        "RNase inhibitor; 72 C 10 min", MaterialKind.CONDITION)
    rt_buffer = graph.input("rt_buffer", "RT mix: 25 mM Tris-HCl pH 8.3, 30 mM NaCl, 2.5 mM MgCl2, "
                            "1 mM GTP, 8 mM DTT; 42 C 90 min, 10x(50 C 2 min, 42 C 2 min), "
                            "85 C 5 min", MaterialKind.CONDITION)
    enzyme = graph.input("reverse_transcriptase", "Maxima H Minus reverse transcriptase (2 U/uL)",
                         MaterialKind.CATALYST)
    first = graph.product("first_strand", "first-strand cDNA with untemplated CCC",
                          MaterialKind.MOLECULAR, state=first_scene)
    graph.react("reverse_transcription", "Lyse and reverse-transcribe from anchored oligo(dT)30VN",
                StepKind.REVERSE_TRANSCRIPTION,
                ((mrna, InputRole.TEMPLATE), (primer, InputRole.STRUCTURAL),
                 (lysis, InputRole.CONDITION), (rt_buffer, InputRole.CONDITION),
                 (enzyme, InputRole.CATALYST)), (first,))
    tso = graph.input("tso", "Biotinylated Smartseq3_N8_TSO (2 uM in RT mix)",
                      MaterialKind.OLIGO, state=switch_scene)
    switched = graph.product("switched", "template-switched full-length cDNA",
                             MaterialKind.MOLECULAR, state=template)
    graph.react("template_switch", "Template-switch at the transcript 5-prime end in the same RT",
                StepKind.TEMPLATE_SWITCHING,
                ((first, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL),
                 (enzyme, InputRole.CATALYST)), (switched,))

    # Step 6 (pp8-9): 6 uL KAPA HiFi mix, Fwd 0.5 uM, Rev 0.1 uM; 98/65/72 C, cycles input-dependent.
    # Bases as printed in Supplementary Table 1; the protocol adds two terminal
    # phosphorothioates to each (shown in the oligo table above the workflow).
    forward = graph.input("preamp_forward", "Fwd_PCR_primer (0.5 uM)", MaterialKind.OLIGO,
                          state=Construct([Segment("Fwd_PCR_primer", P.ss.SS3_FWD_PCR)]))
    reverse = graph.input("preamp_reverse", "Rev_PCR_primer (0.1 uM)", MaterialKind.OLIGO,
                          state=Construct([Segment("Rev_PCR_primer", P.ss.SS3_REV_PCR)]))
    polymerase = graph.input("preamp_polymerase", "KAPA HiFi HotStart polymerase; 98 C 3 min, "
                             "input-dependent cycles of 98 C 20 s / 65 C 30 s / 72 C 4 min, "
                             "72 C 5 min", MaterialKind.CATALYST)
    amplified = graph.product("amplified", "amplified full-length cDNA",
                              MaterialKind.MOLECULAR, state=template)
    graph.react("preamp", "Pre-amplify cDNA with the asymmetric primer pair", StepKind.PCR,
                ((switched, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (amplified,), geometry=PCRGeometry.exact(
                    template, PCRPrimer("Fwd PCR primer", P.ss.SS3_FWD_PCR),
                    PCRPrimer("Rev PCR primer", P.ss.SS3_REV_PCR),
                    statement="The published asymmetric forward and reverse primers face inward around the full-length template-switched cDNA."))
    purified = graph.product("purified_cdna", "bead-purified amplified cDNA",
                             MaterialKind.MOLECULAR, state=template)
    _cleanup(graph, "cdna_cleanup", "Purify cDNA with 0.6:1 PEG beads (protocol: preferable but optional)",
             (amplified,), (purified,),
             "supernatant removed on the magnet; primer and short-product removal is not quantified")

    # Step 12 (pp14-16): 100 pg cDNA, Nextera XT Amplicon Tagment Mix, 55 C 10 min, 0.2% SDS strip.
    transposome = graph.input("transposome", "Nextera XT Amplicon Tagment Mix (0.08 uL per 1 uL "
                              "of 100 pg/uL normalized cDNA; TDE1 stated to work as well)",
                              MaterialKind.ADAPTER)
    td = graph.input("td_buffer", "TD buffer: 10 mM Tris-HCl pH 7.5, 5 mM MgCl2, 5% DMF; 55 C 10 min",
                     MaterialKind.CONDITION)
    tn5 = graph.input("tn5", "Tn5 transposase activity in the Amplicon Tagment Mix",
                      MaterialKind.CATALYST)
    tagged = graph.product("tagmented", "heterologous s5/s7-tagmented cDNA fragments",
                           MaterialKind.MOLECULAR, state=P.tagmented("s5", "s7"))
    waste = graph.product("tagment_waste", "s5/s5, s7/s7 and other non-library Tn5 products",
                          MaterialKind.WASTE)
    graph.react("tagment", "Tagment normalized cDNA with Nextera s5/s7 adapters",
                StepKind.TAGMENTATION,
                ((purified, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (td, InputRole.CONDITION), (tn5, InputRole.CATALYST)), (tagged, waste))
    graph.discard(waste, "lacks the heterologous s5/s7 end combination required for amplification")
    sds = graph.input("sds", "0.5 uL 0.2% SDS, 5 min", MaterialKind.CONDITION)
    released = graph.product("tn5_released", "Tn5-released tagmented fragments",
                             MaterialKind.MOLECULAR, state=P.tagmented("s5", "s7"))
    graph.react("strip_tn5", "Strip Tn5 from the tagmented DNA with SDS", StepKind.DEPROTEINIZATION,
                ((tagged, InputRole.TEMPLATE), (sds, InputRole.CONDITION)), (released,),
                changes_molecule=False)

    # Step 12.7-12.10 (pp16-17): custom S50X/N70X 0.1 uM each, Phusion HF; 72 C 3 min, then 12 cycles.
    index_polymerase = graph.input("index_polymerase", "Phusion HF DNA polymerase (0.01 U/uL)",
                                   MaterialKind.CATALYST)
    filled = graph.product("gap_filled", "gap-filled heterologous Nextera fragments",
                           MaterialKind.MOLECULAR, state=P.gap_filled())
    graph.react("gap_fill", "Fill both nine-base Tn5 gaps (72 C 3 min before denaturation)",
                StepKind.EXTENSION,
                ((released, InputRole.TEMPLATE), (index_polymerase, InputRole.CATALYST)), (filled,))
    five = P.ss.smartseq3_library(PROTOCOL)
    internal = P.ss.smartseq3_library(PROTOCOL, five_prime=False)
    i5 = graph.input("index_forward", "Custom S50X Nextera index primer (0.1 uM; member unassigned)",
                     MaterialKind.OLIGO)
    i7 = graph.input("index_reverse", "Custom N70X Nextera index primer (0.1 uM; member unassigned)",
                     MaterialKind.OLIGO)
    library = graph.product("indexed_library", "SMART-seq3 5-prime library molecules",
                            MaterialKind.LIBRARY, state=five)
    other = graph.product("indexed_internal", "SMART-seq3 internal-fragment library molecules",
                          MaterialKind.LIBRARY, state=internal)
    graph.react("index_pcr", "Index PCR: 98 C 3 min, 12x(98 C 10 s, 55 C 30 s, 72 C 30 s), 72 C 5 min",
                StepKind.PCR,
                ((filled, InputRole.TEMPLATE), (i5, InputRole.STRUCTURAL),
                 (i7, InputRole.STRUCTURAL), (index_polymerase, InputRole.CATALYST)),
                (library, other), geometry=PCRGeometry.exact(
                    five, PCRPrimer("S50X", P.nx.n5xx_primer("N" * 8)),
                    PCRPrimer("N70X", P.nx.n7xx_primer("N" * 8)),
                    statement="Custom Nextera-architecture index primers face inward from the s5 and s7 entry points; the eight-N window is illustrative (the protocol recommends 10 bp or 8 bp)."))
    final = graph.product("library", "Pooled, bead-purified SMART-seq3 5-prime library molecules",
                          MaterialKind.LIBRARY, state=five)
    final_internal = graph.product("internal_library",
                                   "Pooled, bead-purified SMART-seq3 internal-fragment molecules",
                                   MaterialKind.LIBRARY, state=internal)
    _cleanup(graph, "library_cleanup", "Pool wells and purify with 0.6:1 PEG beads",
             (library, other), (final, final_internal),
             "supernatant removed on the magnet; optional 550 bp-2 kb gel selection not drawn")
    primers = tuple(P.ss.SS3_SEQ_PRIMERS)
    for material, construct in ((final, five), (final_internal, internal)):
        graph.sequence(material, SequencingPlan(library=construct, primers=primers,
                       acquired_roles=tuple(p.role for p in primers)))
    graph.declare_reviewed_scope(tuple(graph_ids(graph)), SOURCE)
    return graph.require_publishable()


def graph_ids(graph):
    """Every declared material and step; the whole graph is in this pass's scope."""
    return [*graph.materials, *graph.steps]


def main() -> None:
    OUT.write_text(P.render_page(PROTOCOL, graph_factory=protocol_graph), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
