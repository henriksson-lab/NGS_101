#!/usr/bin/env python3
"""Shared drawing implementation for the three individual STRT pages."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import strt as S
from adapter_ligation import AdapterLigationGeometry
from chemdraw import (Construct, Scene, Segment, complement_segments, oligo, panel,
                      protocol_graph_panels, revcomp)
from endprep import dA_tailed_scene, repair_and_dA_tail
import illumina as il
import nextera as nx
from page import head, info, table
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
import seqprimers as sp


PAPERS = {
    S.STRT: ("10.1101/gr.110882.110", "STRT-seq"),
    S.C1: ("10.1038/nmeth.2772", "STRT-C1"),
    S.TWO_I: ("10.1038/s41598-017-16546-4", "STRT-seq-2i"),
}


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def preamble(protocol: str, *, complete: bool = False) -> str:
    doi, title = PAPERS[protocol]
    notes = {
        S.STRT: "01_strt-seq.html",
        S.C1: "../strt-seq__10.1101+gr.110882.110/01_strt-seq.html",
        S.TWO_I: "../strt-seq__10.1101+gr.110882.110/01_strt-seq.html",
    }[protocol]
    caveat = ""
    if protocol == S.STRT:
        caveat = '<div class="caveat"><b>Inferred primer site.</b> The custom Read 1 primer assignment remains unverified in this capped pass. Its landing site is reconstructed from the published construct; the original supplement is deferred.</div>'
    elif protocol == S.C1:
        caveat = '<div class="caveat"><b>Partial methods.</b> The accessible Supplementary Note defines the molecule and says it is sequenced without library PCR; the paywalled Online Methods contain the unavailable selection details.</div>'
    elif protocol == S.TWO_I:
        caveat = '<div class="caveat"><b>Inferred junction.</b> The published P2 PCR primer overlaps the first three bases of the subarray barcode. The dotted barcode/P7 end shows the full-barcode outcome; the paper does not state whether proofreading instead trims the primer or overwrites those bases.</div>'
    complete_attr = ' data-protocol-complete="true"' if complete else ""
    return f'''<div class="wrap"{complete_attr}>
<h1>{title} &mdash; 5'-end single-cell RNA sequencing</h1>
<p class="research-notes"><a href="{notes}">Research notes</a></p>
{info(f'Defining source: <a href="https://doi.org/{doi}">doi:{doi}</a>.')}
{caveat}
'''


def oligos(protocol: str) -> str:
    if protocol == S.STRT:
        rows = [
            oligo("STRT-V3-T30", [seg("SMART handle", S.H, "tso"), seg("SalI junction", S.STRT_DT_LINK), seg("T30", "T" * S.STRT_DT_LEN), seg("VN", "VN", placeholder=True)], mods="biotin"),
            oligo("STRT-V2-n TSO", [seg("SMART handle", S.H, "tso"), seg("spacer", S.STRT_TSO_LINK, "r1"), seg("cell barcode", "B" * 6, "cbc", placeholder=True), seg("rGrGrG", "rGrGrG", placeholder=True)]),
            oligo("STRT-PCR", [seg("SMART handle", S.H, "tso")], mods="biotin"),
            oligo("P2 adapter, top", [seg("P7 fragment", il.P7[:22], "p7"), seg("stem + T", il.TRUSEQ_READ2[-12:], "r2")]),
            oligo("Library PCR 1", [seg("P5 fragment", il.P5[:25], "p5"), seg("SMART handle", S.H, "tso")]),
            oligo("Library PCR 2", [seg("P7 fragment", il.P7[:22], "p7")]),
        ]
    elif protocol == S.C1:
        rows = [
            oligo("C1-P1-T31", [seg("C1 handle", S.C1_HANDLE, "p5"), seg("PvuI junction", S.C1_DT_LINK), seg("T31", "T" * S.C1_DT_LEN)], mods="Bio"),
            oligo("C1-P1-RNA-TSO", [seg("C1 handle", S.C1_HANDLE, "p5"), seg("UMI", "N" * 5, "umi", placeholder=True), seg("rGrGrG", "rGrGrG", placeholder=True)], mods="Bio; all-ribo"),
            oligo("C1-P1-PCR-2", [seg("extra G", "G"), seg("C1 handle", S.C1_HANDLE, "p5")], mods="Bio"),
            oligo("C1-TN5-x", S.tn5_indexed_top()),
            oligo("C1-TN5-U", [seg("paired tail'", S.TN5_U, "me")], mods="5Phos"),
        ]
    else:
        rows = [
            oligo("STRT-P1-T31", [seg("C1 handle", S.C1_HANDLE, "p5"), seg("CG", S.C1_DT_LINK), seg("T31", "T" * S.C1_DT_LEN)], mods="Bio"),
            oligo("P1B-UMI-RNA-TSO", [seg("P1B", S.P1B, "r1"), seg("UMI", "N" * 6, "umi", placeholder=True), seg("rGrGrG", "rGrGrG", placeholder=True)], mods="Bio; all-ribo"),
            oligo("DI-P1A-idx-P1B", [seg("P5", il.P5, "p5"), seg("well index", "W" * 5, "cbc", placeholder=True), seg("P1B minus terminal T", S.P1B[:-1], "r1")], mods="Bio"),
            oligo("STRT-Tn5-Idx", S.tn5_indexed_top()),
            oligo("STRT-TN5-U", [seg("paired tail'", S.TN5_U, "me")], mods="5Phos"),
            oligo("P1/P2 second-PCR primers", [seg("P5 fragment", il.P5[:24], "p5"), seg(" / ", "", placeholder=True), seg("P7", il.P7, "p7")]),
        ]
    return "<h2>Key oligos</h2>\n<seq>\n" + "\n".join(rows) + "\n</seq>"


def rt_scenes(protocol: str) -> tuple[Scene, Scene]:
    body = seg("body", "XXXXXXXX...XXXXXXXX", placeholder=True)
    if protocol == S.STRT:
        poly = "A" * S.STRT_DT_LEN
        primer = [seg("handle", S.H, "tso"), seg("CGAC", S.STRT_DT_LINK),
                  seg("dT", "T" * S.STRT_DT_LEN)]
        tso = [seg("TSO handle", S.H, "tso"), seg("spacer", S.STRT_TSO_LINK, "r1"), seg("cell barcode", "B" * 6, "cbc", placeholder=True), seg("GGG", "GGG", "tso")]
    elif protocol == S.C1:
        poly = "A" * S.C1_DT_LEN
        primer = [seg("handle", S.C1_HANDLE, "p5"), seg("CG", S.C1_DT_LINK), seg("dT", "T" * S.C1_DT_LEN)]
        tso = [seg("TSO handle", S.C1_HANDLE, "p5"), seg("UMI", "N" * 5, "umi", placeholder=True), seg("GGG", "GGG", "tso")]
    else:
        poly = "A" * S.C1_DT_LEN
        primer = [seg("handle", S.C1_HANDLE, "p5"), seg("CG", S.C1_DT_LINK), seg("dT", "T" * S.C1_DT_LEN)]
        tso = [seg("P1B", S.P1B, "r1"), seg("UMI", "N" * 6, "umi", placeholder=True), seg("GGG", "GGG", "tso")]
    mrna = [body, seg("poly(A)", poly)]
    first = [*primer, *complement_segments([body]), seg("CCC", "CCC", "tso")]
    p = Scene(); p.strand("mRNA", mrna); p.anneal("oligo-dT", primer, to="mRNA", pair=("dT", "poly(A)"), mod5="biotin"); p.arrow("oligo-dT", "reverse transcriptase")
    sw = Scene(); sw.strand("mRNA", mrna); sw.anneal("first-strand cDNA", first, to="mRNA", pair=("dT", "poly(A)"), mod5="biotin"); sw.anneal("TSO", tso, to="first-strand cDNA", pair=("GGG", "CCC"), above=True); sw.arrow("first-strand cDNA", "template switch")
    return p, sw


def indexed_tn5_scene() -> Scene:
    top = S.tn5_indexed_top()
    sc = Scene(); sc.strand("indexed top", top, label="indexed strand"); sc.anneal("TN5-U", [seg("paired tail'", S.TN5_U, "me")], to="indexed top", pair=("paired tail'", "GCGTC"), label="TN5-U", mod5="p")
    sc.mark("indexed top", "subarray barcode", "8-nt Tn5 barcode")
    return sc


def strt_p2_scene() -> Scene:
    """The original STRT TSO-end fragment after P2-adapter ligation."""
    top = [seg("TSO end", S.H, "tso"),
           seg("cell barcode", "B" * 6, "cbc", placeholder=True),
           seg("5′ cDNA fragment", "XXXXXXXX...XXXXXXXX", placeholder=True),
           seg("P2 adapter'", revcomp(S.P2_TOP), "p7")]
    sc = Scene.duplex(top)
    sc.junction("top", "5′ cDNA fragment", "P2 adapter'", "P2 ligation")
    return sc


def strt_p2_adapter_scene() -> Scene:
    """Source-exact P2 fork before ligation to the dA-tailed cDNA fragment."""
    top = [seg("paired stem", S.P2_TOP[:-1], "r2"),
           seg("3′ T", S.P2_TOP[-1], "r2")]
    bottom = [seg("paired stem reverse complement", revcomp(S.P2_TOP[:-1]), "r2")]
    sc = Scene(); sc.strand("top", top, label="P2 adapter top")
    sc.anneal("bottom", bottom, to="top",
              pair=("paired stem reverse complement", "paired stem"),
              label="P2 adapter bottom", mod5="p")
    sc.overhang("top", "3′ T", polarity="3-prime", owner="adapter",
                text="adapter 3′-T overhang")
    return sc


def _front(graph: ProtocolGraph, protocol: str):
    prime, switched_state = rt_scenes(protocol)
    mrna = graph.input("mrna", "polyadenylated RNA from one cell",
                       MaterialKind.BIOLOGICAL, state=prime, primary=True)
    if protocol == S.STRT:
        dt_state = Construct([seg("SMART handle", S.H), seg("SalI junction", S.STRT_DT_LINK),
                              seg("T30", "T" * S.STRT_DT_LEN),
                              seg("VN", "VN", placeholder=True)])
        tso_state = Construct([seg("SMART handle", S.H),
                               seg("spacer", S.STRT_TSO_LINK),
                               seg("cell barcode", "B" * S.STRT_BC_LEN,
                                   placeholder=True), seg("GGG", "GGG")])
    elif protocol == S.C1:
        dt_state = Construct([seg("C1 handle", S.C1_HANDLE),
                              seg("PvuI junction", S.C1_DT_LINK),
                              seg("T31", "T" * S.C1_DT_LEN)])
        tso_state = Construct([seg("C1 handle", S.C1_HANDLE),
                               seg("UMI", "N" * S.C1_UMI_LEN, placeholder=True),
                               seg("GGG", "GGG")])
    else:
        dt_state = Construct([seg("C1 handle", S.C1_HANDLE),
                              seg("CG", S.C1_DT_LINK),
                              seg("T31", "T" * S.C1_DT_LEN)])
        tso_state = Construct([seg("P1B", S.P1B),
                               seg("UMI", "N" * S.TWO_I_UMI_LEN, placeholder=True),
                               seg("GGG", "GGG")])
    dt = graph.input("oligo_dt", "biotinylated oligo-dT primer",
                     MaterialKind.OLIGO, state=dt_state)
    rt = graph.input("reverse_transcriptase", "template-switching reverse transcriptase",
                     MaterialKind.CATALYST)
    first = graph.product("first_strand", "first-strand cDNA",
                          MaterialKind.MOLECULAR, state=prime)
    graph.react("reverse_transcription", "Reverse-transcribe polyadenylated RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((mrna, InputRole.TEMPLATE), (dt, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (first,))
    tso = graph.input("tso", "template-switching oligo", MaterialKind.OLIGO,
                      state=tso_state)
    switched = graph.product("template_switched_cdna", "template-switched cDNA",
                             MaterialKind.MOLECULAR, state=switched_state)
    graph.react("template_switch", "Copy the TSO at the transcript 5′ end",
                StepKind.TEMPLATE_SWITCHING,
                ((first, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (switched,))
    return switched, ["mrna", "oligo_dt", "reverse_transcriptase",
                      "reverse_transcription", "tso", "template_switch"]


def _pool(graph: ProtocolGraph, material, state, *, identifier: str, label: str):
    peers = graph.input(identifier + "_peers", label, MaterialKind.MOLECULAR, state=state)
    pooled = graph.product(identifier, label.replace("other ", "pooled "),
                           MaterialKind.MOLECULAR, state=state)
    graph.react(identifier + "_step", "Pool after molecular indexing", StepKind.POOL,
                ((material, InputRole.TEMPLATE), (peers, InputRole.TEMPLATE)),
                (pooled,), changes_molecule=False)
    return pooled, [identifier + "_peers", identifier + "_step"]


def _tagment(graph: ProtocolGraph, template, state, reviewed, *, inferred=False,
             source=""):
    transposome = graph.input(
        "indexed_transposome", "indexed STRT Tn5 transposome", MaterialKind.ADAPTER,
        state=indexed_tn5_scene(), inferred=inferred,
        source=(source + " (transposase source/selection details unavailable)"
                if inferred else None))
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST,
                      inferred=inferred,
                      source=(source + " (transposase source unavailable)"
                              if inferred else None))
    tagged = graph.product("tagmented_cdna", "indexed tagmented cDNA",
                           MaterialKind.MOLECULAR, state=state, inferred=inferred,
                           source=(source + " (selected tagmentation product reconstructed)"
                                   if inferred else None))
    waste = graph.product("tagment_waste", "other tagmentation products",
                          MaterialKind.WASTE)
    graph.react("tagment", "Tagment amplified cDNA with the indexed adapter",
                StepKind.TAGMENTATION,
                ((template, InputRole.TEMPLATE),
                 (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, waste))
    graph.discard(waste, "does not retain the selected transcript-5′ topology")
    reviewed.extend(("indexed_transposome", "tn5", "tagment"))
    return tagged


def graph_for(protocol: str) -> ProtocolGraph:
    if protocol != S.STRT:
        raise ValueError(f"unknown STRT protocol: {protocol!r}")
    source = "https://doi.org/" + PAPERS[protocol][0]
    graph = ProtocolGraph(protocol, authority=source)
    switched, reviewed = _front(graph, protocol)

    pooled, pool_reviewed = _pool(
        graph, switched, switched.state, identifier="pooled_cdna",
        label="other well-barcoded cDNA from the plate")
    reviewed.extend(pool_reviewed)
    beads = graph.input("cleanup_beads", "carboxylate cleanup beads",
                        MaterialKind.SELECTION)
    clean = graph.product("clean_cdna", "cleaned pooled cDNA",
                          MaterialKind.MOLECULAR, state=switched.state)
    clean_waste = graph.product("cleanup_waste", "unretained cleanup material",
                                MaterialKind.WASTE)
    graph.react("pool_cleanup", "Clean the pooled cDNA", StepKind.CLEANUP,
                ((pooled, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                (clean, clean_waste))
    graph.discard(clean_waste, "not retained by the post-pooling bead cleanup")
    primer = PCRPrimer("STRT-PCR", S.H)
    preamp_primer = graph.input("preamp_primer", "biotinylated STRT-PCR primer",
                                MaterialKind.OLIGO,
                                state=Construct([seg("STRT-PCR", S.H)]))
    polymerase = graph.input("preamp_polymerase", "Advantage2 polymerase",
                             MaterialKind.CATALYST)
    preamp_state = S.strt_preamp_product()
    amplified = graph.product("amplified_cdna", "biotinylated amplified cDNA",
                              MaterialKind.MOLECULAR, state=preamp_state)
    graph.react("preamp", "Amplify pooled cDNA: 20 cycles with the single STRT-PCR primer",
                StepKind.PCR,
                ((clean, InputRole.TEMPLATE),
                 (preamp_primer, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (amplified,),
                geometry=PCRGeometry.symmetric(
                    preamp_state, primer,
                    statement="STRT-PCR binds the shared SMART handle at both ends."))
    capture = graph.input("streptavidin_beads", "streptavidin beads",
                          MaterialKind.SELECTION)
    bound = graph.product("bead_bound_cdna", "streptavidin-bound amplified cDNA",
                          MaterialKind.MOLECULAR, state=preamp_state)
    unbound = graph.product("capture_waste", "unbound amplification material",
                            MaterialKind.WASTE)
    graph.react("immobilize", "Immobilize biotinylated cDNA", StepKind.AFFINITY_CAPTURE,
                ((amplified, InputRole.TEMPLATE),
                 (capture, InputRole.SELECTION)), (bound, unbound),
                changes_molecule=False)
    graph.discard(unbound, "not retained by streptavidin capture")
    dnase = graph.input("dnase_i", "DNase I", MaterialKind.CATALYST)
    fragments = graph.product("fragments", "bead-bound cDNA fragments",
                              MaterialKind.MOLECULAR,
                              state=Construct([seg("5′ cDNA fragment", "X" * 36,
                                                   placeholder=True)]))
    graph.react("fragment", "Fragment bead-bound cDNA", StepKind.FRAGMENTATION,
                ((bound, InputRole.TEMPLATE), (dnase, InputRole.CATALYST)),
                (fragments,))
    repair = graph.input("end_repair_enzymes", "end-repair enzyme mix",
                         MaterialKind.CATALYST)
    repaired = graph.product("repaired_fragments", "end-repaired cDNA fragments",
                             MaterialKind.MOLECULAR, state=fragments.state)
    graph.react("end_repair", "End-repair cDNA fragments", StepKind.END_REPAIR,
                ((fragments, InputRole.TEMPLATE),
                 (repair, InputRole.CATALYST)), (repaired,))
    da_enzyme = graph.input("da_tailing_enzyme", "dA-tailing polymerase",
                            MaterialKind.CATALYST)
    da_product = repair_and_dA_tail(fragments.state)
    tailed = graph.product("da_tailed_fragments", "3′-dA-tailed fragments",
                           MaterialKind.MOLECULAR,
                           state=dA_tailed_scene(da_product))
    graph.react("da_tail", "Add 3′ dA overhangs", StepKind.TAILING,
                ((repaired, InputRole.TEMPLATE),
                 (da_enzyme, InputRole.CATALYST)), (tailed,))
    adapter = graph.input("p2_adapter", "P2 forked adapter", MaterialKind.ADAPTER,
                          state=strt_p2_adapter_scene())
    ligase = graph.input("dna_ligase", "T4 DNA ligase", MaterialKind.CATALYST)
    ligated = graph.product("p2_ligated_fragments", "P2-ligated cDNA fragments",
                            MaterialKind.MOLECULAR, state=strt_p2_scene())
    geometry = AdapterLigationGeometry.known(
        "cohesive", dA_tailed_scene(da_product), strt_p2_adapter_scene(),
        statement="The printed P2 adapter exposes a 3′-T for the cDNA 3′-dA.")
    graph.react("ligate_p2", "Ligate the P2 adapter", StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (ligated,), geometry=geometry)
    sali = graph.input("sali", "SalI-HF", MaterialKind.CATALYST)
    selected = graph.product("tso_end_fragments", "bead-bound TSO-end fragments",
                             MaterialKind.MOLECULAR, state=strt_p2_scene())
    released = graph.product("released_three_prime_fragments",
                             "released oligo-dT-end fragments", MaterialKind.WASTE)
    graph.react("sali_select", "Release oligo-dT-end fragments with SalI-HF",
                StepKind.DIGESTION,
                ((ligated, InputRole.TEMPLATE),
                 (sali, InputRole.CATALYST)), (selected, released))
    graph.discard(released, "released from the beads at the oligo-dT SalI site")
    final = S.strt_library()
    f_primer = PCRPrimer("library PCR 1", S.STRT_LIB_PCR1)
    r_primer = PCRPrimer("library PCR 2", S.STRT_LIB_PCR2)
    forward = graph.input("library_pcr1", f_primer.name, MaterialKind.OLIGO,
                          state=Construct([seg(f_primer.name, f_primer.sequence)]))
    reverse = graph.input("library_pcr2", r_primer.name, MaterialKind.OLIGO,
                          state=Construct([seg(r_primer.name, r_primer.sequence)]))
    lib_polymerase = graph.input("library_polymerase", "Phusion polymerase",
                                 MaterialKind.CATALYST)
    unselected_library = graph.product("unselected_library", "amplified STRT library",
                                       MaterialKind.MOLECULAR, state=final)
    graph.react("library_pcr", "Amplify selected TSO-end fragments: 12 cycles",
                StepKind.PCR,
                ((selected, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL),
                 (lib_polymerase, InputRole.CATALYST)), (unselected_library,),
                geometry=PCRGeometry.exact(
                    final, f_primer, r_primer,
                    statement="The defining paper prints both library-PCR primers."))
    size_reagent = graph.input("size_selection", "200–400 bp size-selection matrix",
                               MaterialKind.SELECTION)
    library = graph.product("library", "size-selected STRT-seq library",
                            MaterialKind.LIBRARY, state=final)
    size_waste = graph.product("size_waste", "out-of-range library molecules",
                               MaterialKind.WASTE)
    graph.react("size_select", "Select 200–400 bp library molecules",
                StepKind.SELECTION,
                ((unselected_library, InputRole.TEMPLATE),
                 (size_reagent, InputRole.SELECTION)), (library, size_waste))
    graph.discard(size_waste, "outside the reported 200–400 bp library interval")
    graph.sequence(library, SequencingPlan(
        library=final, primers=(S.STRT_READ1,), acquired_roles=("Read 1",),
        read_lengths={"Read 1": 55}, limitation="Custom primer is a derived reference site; original supplement remains deferred."))
    reviewed.extend([
        "cleanup_beads", "pool_cleanup", "preamp_primer", "preamp_polymerase",
        "preamp", "streptavidin_beads", "immobilize", "dnase_i", "fragment",
        "end_repair_enzymes", "end_repair", "da_tailing_enzyme", "da_tail",
        "p2_adapter", "dna_ligase", "ligate_p2", "sali", "sali_select",
        "library_pcr1", "library_pcr2", "library_polymerase", "library_pcr",
        "size_selection", "size_select", "library",
    ])
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph() -> ProtocolGraph:
    """The defining STRT-seq chemistry for source and molecular audits."""
    return graph_for(S.STRT)


def sequencing(protocol: str) -> str:
    lib = S.library(protocol)
    if protocol == S.STRT:
        primers = (S.STRT_READ1,)
        hit = sp.locate(lib, S.STRT_READ1)
        rows = [("Read 1", "~55", "6-nt cell barcode &middot; GGG &middot; transcript 5' end")]
        note = "The custom primer assignment remains unverified; the drawn reference site ends immediately before the barcode. The original article reports 55-base reads and three to six C residues after the barcode; GGG in the opposite construct orientation is representative rather than a fixed-length tail."
    elif protocol == S.C1:
        primers = S.C1_PRIMERS
        hits = [sp.locate(lib, p) for p in S.C1_PRIMERS]
        if any(h is None for h in hits): raise ValueError("C1 read site missing")
        rows = [("Read 1", "50", "5-nt UMI &middot; GGG &middot; transcript 5' end"), ("Index 1", "8", "Tn5 cell barcode")]
        note = "Read 1 and index-primer sequences are not printed in the accessible paper materials; their landing sites follow from the published molecule."
    else:
        primers = S.TWO_I_INDEX_PRIMERS
        hits = [sp.locate(lib, p) for p in S.TWO_I_INDEX_PRIMERS]
        if any(h is None for h in hits): raise ValueError("2i index site missing")
        rows = [("Read 1", "45 or 48", "6-nt UMI &middot; GGG &middot; transcript 5' end"), ("Index 1", "8", "subarray index"), ("Index 2", "5", "well index")]
        note = "The paper reports 45 cycles in Methods and 48 in its supplement. Its printed Read-1 primer contains six N opposite a five-base well index; the constant 3' P1B portion still fixes the read start."
        # Scene enforces the printed degenerate primer against the published P5+5-index+P1B template.
        sc = Scene.duplex([seg("P5", il.P5, "p5"),
                           seg("well index", "W" * 5, "cbc", placeholder=True),
                           seg("P1B", S.P1B, "r1"),
                           seg("UMI", "N" * 6, "umi", placeholder=True)])
        sc.anneal("DI-Read1-Seq", [seg("P1B primer (paired 3' end)", S.P1B, "r1")],
                  to="bottom", pair=("P1B primer (paired 3' end)", "P1B'"))
        sc.arrow("DI-Read1-Seq", "UMI, then transcript")
        note += panel(sc.rows(), cls="small", caption="The printed 6-N Read-1 primer spans the last P5 base plus the five-base well index, then pairs exactly over P1B.")
    if protocol == S.STRT and hit is None:
        raise ValueError("STRT Read 1 site missing")
    specs = [sp.read_spec({"Index 1": "Index 1 (i7)",
                           "Index 2": "Index 2 (i5)"}.get(r, r),
                          c, content.replace("&middot;", "·"),
                          "extract_then_include" if r == "Read 1" else "extract")
             for r, c, content in rows]
    return ('<h2>Read layout</h2>' + sp.diagram(lib, primers)
            + table(("Read", "Cycles", "Content"), rows) + info(note)
            + sp.read_spec_table(specs, "Read roles and cycle counts are reported by the defining source.")
            + '</div>')


def render_page(protocol: str, out: Path, *, graph=None) -> None:
    title = PAPERS[protocol][1]
    graph = graph if graph is not None else graph_for(protocol)
    workflow = "<h2>Library generation</h2>\n" + protocol_graph_panels(graph)
    out.write_text("\n".join([head(f"{title} library chemistry"),
                              preamble(protocol, complete=True), oligos(protocol),
                              workflow, sequencing(protocol)]), encoding="utf-8")
    print(f"wrote {out}  ({out.stat().st_size:,} bytes)")
