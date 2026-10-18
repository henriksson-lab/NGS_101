"""Construct model for 10x Chromium Single Cell 3' Gene Expression v2 (CG000108)."""
from __future__ import annotations

import illumina as il
import rt
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, revcomp
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import Disposition, InputRole, MaterialKind, SequencingPlan, StepKind, Terminal
from protocol_model import molecule, role
from protocol_graph import ProtocolGraph, SourceClaim, SourceDiscrepancy

USER_GUIDE = ('<a href="https://assets.ctfassets.net/an68im79xiti/UhAMGmlaEMmYMaA4A4Uwa/'
              '274a813b81e42cba81345d49380432d7/CG00052_SingleCell3_ReagentKitv2UserGuide_RevD.pdf">'
              '10x Genomics User Guide CG00052 Rev D</a> (2018; cleanup and size-selection steps)')

PARTIAL_R1 = il.TRUSEQ_READ1[11:]
TSO_DNA = rt.SMART_HANDLE + "ACAT"
CELL_BARCODE = feature("cell_barcode", "cell_barcode", "whitelist",
                       whitelist="10x-chromium-3prime-v2")
UMI = feature("umi", "umi", "random")
SAMPLE_INDEX = feature("sample_index_i7", "sample_index", "unknown")
FEATURES = {"cell barcode": CELL_BARCODE, "UMI": UMI,
            "i7 reverse complement": SAMPLE_INDEX}


def seg(name, top, tag=None, **kw):
    kw.setdefault("feature", FEATURES.get(name))
    return Segment(name=name, top=top, tag=tag, **kw)


def bead_oligo():
    return [seg("Partial Read 1", PARTIAL_R1, "r1"),
            seg("cell barcode", "B" * 16, "cbc", placeholder=True),
            seg("UMI", "U" * 10, "umi", placeholder=True),
            seg("poly(dT)30", "T" * 30), seg("V", "V", placeholder=True),
            seg("N", "N", placeholder=True)]


def tso():
    return [seg("SMART handle", rt.SMART_HANDLE, "tso"), seg("ACAT", "ACAT", "tso"),
            seg("rGrGrG", "GGG", "tso")]


def capture_scene():
    mrna = [seg("mRNA", "X" * 32, placeholder=True), seg("poly(A)", "A" * 30)]
    sc = Scene(); sc.strand("mRNA", mrna, label="mRNA")
    sc.anneal("bead primer", bead_oligo(), to="mRNA", pair=("poly(dT)30", "poly(A)"),
              label="released gel-bead primer", unpaired=("V", "N"))
    sc.arrow("bead primer", "reverse transcriptase")
    return sc


def first_strand():
    return Construct([*bead_oligo(), seg("antisense cDNA", "X" * 34, placeholder=True),
                      seg("CCC", rt.UNTEMPLATED_TAIL, "tso")], name="v2 first strand")


def switch_scene():
    first = first_strand(); sc = Scene(); sc.strand("first", list(first), label="first strand")
    sc.anneal("TSO", tso(), to="first", pair=("rGrGrG", "CCC"), label="TSO", above=True)
    sc.arrow("first", "RT switches template")
    return sc


def amplified_cdna():
    return Construct([*bead_oligo(), seg("antisense cDNA", "X" * 34, placeholder=True),
                      seg("CCC", rt.UNTEMPLATED_TAIL, "tso"),
                      seg("copied TSO", revcomp(TSO_DNA), "tso")],
                     name="v2 amplified cDNA")


def adapter_top():
    return [seg("Read 2 adapter", il.INDEX1_PRIMER, "r2")]


def adapter_bottom():
    return [seg("stem complement", il.STEM_COMPLEMENT, "r2"), seg("3' T", "T", "r2")]


def adapter_scene():
    """Adapter PN-220026 as printed in CG000108 Rev A, Protocol 3.3, with its 3'-T end."""
    sc = Scene(); sc.strand("top", adapter_top(), label="Read 2 adapter")
    sc.anneal("bottom", adapter_bottom(), to="top", pair=("stem complement", "Read 2 adapter"),
              label="short strand")
    sc.ligation_end("bottom", "3' T", kind="cohesive", polarity="3-prime", owner="adapter",
                    text="3'-T ligation overhang")
    return sc


def read2_adapter_end():
    return AdapterLigationGeometry.known(
        "cohesive", adapter_scene(),
        statement=("CG000108 Rev A, Protocol 3.3 prints Adapter (Read 2) PN-220026 as the "
                   "33-nt top strand paired over 12 nt by 3'-TCTAGCCTTCTCG-5', leaving a "
                   "single 3'-T for the dA-tailed fragment; 5' phosphorylation is not stated."))


SEQ_PRIMERS = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])


def final_library(insert_nt=36):
    lib = Construct([seg("P5", il.P5, "p5"), seg("Read 1 arm", il.TRUSEQ_READ1[4:], "r1"),
        seg("cell barcode", "B" * 16, "cbc", placeholder=True),
        seg("UMI", "U" * 10, "umi", placeholder=True), seg("poly(dT)30", "T" * 30),
        seg("V", "V", placeholder=True), seg("N", "N", placeholder=True),
        seg("cDNA insert", "X" * insert_nt, placeholder=True), seg("dA junction", "A"),
        seg("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2"),
        seg("i7 reverse complement", "I" * 8, "cbc", placeholder=True),
        seg("P7 reverse complement", il.P7_RC, "p7")], name="10x 3' v2 library")
    problems = sp.verify(lib, SEQ_PRIMERS,
                         required_roles=("Read 1", "Index 1 (i7)", "Read 2"))
    if problems: raise ValueError("invalid v2 library: " + "; ".join(problems))
    return lib


def read1_layout():
    return [(x.cycles, x.feature.label) for x in
            sp.feature_spans(final_library(), SEQ_PRIMERS, {"Read 1": 26})]



def switched_cdna():
    """Single-stranded template-switched first strand; cDNA PCR makes the duplex later."""
    sc = Scene()
    sc.strand("switched", list(amplified_cdna()), label="v2 template-switched first-strand cDNA")
    return sc


def cdna_primers():
    return (PCRPrimer("cDNA forward primer", PARTIAL_R1),
            PCRPrimer("cDNA reverse primer", TSO_DNA))


def index_primers():
    # Source prints an N8 sample-index window; it represents one oligo in a well mix.
    return (PCRPrimer("sample-index forward primer", il.NEBNEXT_UNIVERSAL_PRIMER[:49],
                      anneal_nt=13),
            PCRPrimer("sample-index reverse primer",
                      il.P7+"N"*8+il.TRUSEQ_READ2[:21]))


def pre_index_library(insert_nt=36):
    return Construct([*bead_oligo(),
                      seg("cDNA insert", "X"*insert_nt, placeholder=True),
                      seg("dA junction", "A"), *adapter_top()],
                     name="v2 adapter-ligated cDNA fragment")


def _cleanup(g, reviewed, ident, label, material, reagent, wastes, *, state=None,
             kind=MaterialKind.MOLECULAR, kept_id=None, kept_label=None):
    """A CG00052 Rev D bead cleanup: the kept fraction and each removed fraction."""
    beads=g.input(ident+"_beads", reagent, MaterialKind.SELECTION, source=USER_GUIDE)
    kept=g.product(kept_id or ident+"_kept", kept_label or label+" (retained DNA)", kind,
                   state=material.state if state is None else state, source=USER_GUIDE)
    removed=[g.product(ident+"_waste%d" % i, text, MaterialKind.WASTE, source=USER_GUIDE)
             for i, (text, _) in enumerate(wastes, 1)]
    g.react(ident, label, StepKind.CLEANUP,
            ((material,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),
            (kept,*removed), source=USER_GUIDE, changes_molecule=False)
    for waste, (_, why) in zip(removed, wastes):
        g.discard(waste, why)
    reviewed.extend((beads.id, ident, kept.id, *(w.id for w in removed)))
    return kept


SUPERNATANT=("bead supernatant (unbound primers, reagents and fragments below the cut-off)",
             "supernatant removed from the magnet-held SPRIselect beads")
LARGE=("fragments bound at 0.6X SPRIselect (above the upper size cut-off)",
       "first-SPRI beads discarded after the supernatant is transferred")
SMALL=("supernatant after 0.8X SPRIselect (below the lower size cut-off)",
       "supernatant removed from the magnet-held SPRIselect beads")


def _source_pcr(g, reviewed, prefix, label, template, product_state, primers, *, library=False,
                site_template=None):
    """The source prints both primers; locate their sites rather than claiming absence.

    ``site_template`` is the duplex on which both primer sites are drawn when the input
    material is single-stranded (the switched first strand gains its complement in cycle 1).
    """
    forward, reverse=primers
    f=g.input(prefix+"_forward", forward.name, MaterialKind.OLIGO,
              state=Construct([Segment(forward.name, forward.sequence)]))
    r=g.input(prefix+"_reverse", reverse.name, MaterialKind.OLIGO,
              state=Construct([Segment(reverse.name, reverse.sequence)]))
    polymerase=g.input(prefix+"_polymerase", label+" polymerase", MaterialKind.CATALYST)
    product=g.product("index_pcr_product" if library else "amplified", label+" product",
                      MaterialKind.MOLECULAR, state=product_state)
    g.react(prefix, label, StepKind.PCR,
            ((template,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),
             (r,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)), (product,),
            geometry=PCRGeometry.exact(site_template or template.state,forward,reverse,
                statement="Both primer architectures are printed in CG000108 Rev A Fig. 4."))
    reviewed.extend((f.id,r.id,polymerase.id,prefix,product.id))
    return product


def _source_ligation(g, reviewed, tailed):
    adapter=g.input("read2_adapter", "3-prime-dT Read-2 adapter", MaterialKind.ADAPTER,
                    state=adapter_scene())
    ligase=g.input("ligase", "adapter ligase", MaterialKind.CATALYST)
    adapted=g.product("adapted", "Read-2-adapted barcode-bearing fragment",
                      MaterialKind.MOLECULAR,state=pre_index_library())
    g.react("ligate", "Ligate the Read-2 adapter", StepKind.LIGATION,
            ((tailed,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),
             (ligase,InputRole.CATALYST)),(adapted,),geometry=read2_adapter_end())
    reviewed.extend((adapter.id,ligase.id,"ligate"))
    return adapted

def protocol_graph():
    authority = ('Primary chemistry source: <a href="https://www.10xgenomics.com/support">'
                 '10x Genomics Technical Note CG000108 Rev A</a>.')
    graph = ProtocolGraph("10x Chromium Single Cell 3-prime Gene Expression v2", authority=authority)
    # Reaction choices and source coverage belong to this guide-specific declaration.
    reviewed = ["cells", "gems", "partition", "rt_primer", "rt", "reverse_transcribe",
                "tso", "template_switch", "gem_recovery", "recover", "fragmentase",
                "fragment", "end_mix", "end_prep"]
    final = final_library()
    cells = graph.input("cells", "cells containing poly(A) RNA", MaterialKind.BIOLOGICAL,
                        state=molecule("cells", (role("poly(A) RNA"),)), primary=True)
    gems = graph.input("gems", "Chromium GEM partitioning condition", MaterialKind.CONDITION)
    occupied = graph.product("occupied", "occupied GEMs containing cells", MaterialKind.MOLECULAR,
                             state=cells.state)
    failed = graph.product("failed", "empty, multiplet or unresolved GEM products", MaterialKind.MOLECULAR)
    graph.react("partition", "Partition cells with barcoded gel beads", StepKind.PARTITION,
                ((cells, InputRole.TEMPLATE), (gems, InputRole.CONDITION)),
                (occupied, failed), changes_molecule=False)
    graph.mark_terminal(Terminal(failed.id, Disposition.UNRESOLVED,
        "GEM contents are pooled after GEM break. Barcode resolution is a downstream analysis "
        "criterion; the guide does not establish physical removal of these GEM classes here."))
    primer = graph.input("rt_primer", "cell-barcode/UMI oligo(dT) gel-bead primer", MaterialKind.OLIGO,
                         state=Construct(bead_oligo(), name="gel-bead primer"))
    enzyme = graph.input("rt", "template-switching reverse transcriptase", MaterialKind.CATALYST)
    first_state = Scene()
    first_state.strand("first", list(first_strand()), label="first-strand transcript cDNA")
    first = graph.product("first_strand", "first-strand transcript cDNA", MaterialKind.MOLECULAR,
                          state=first_state)
    graph.react("reverse_transcribe", "Capture and reverse-transcribe poly(A) RNA", StepKind.REVERSE_TRANSCRIPTION,
                ((occupied, InputRole.TEMPLATE), (primer, InputRole.STRUCTURAL),
                 (enzyme, InputRole.CATALYST)), (first,))
    switch_oligo = graph.input("tso", "template-switch oligo", MaterialKind.OLIGO,
                              state=Construct(tso(), name="template-switch oligo"))
    switched = graph.product("switched", "template-switched cell-barcoded cDNA", MaterialKind.MOLECULAR,
                             state=switched_cdna())
    graph.react("template_switch", "Install the opposite cDNA amplification handle", StepKind.TEMPLATE_SWITCHING,
                ((first, InputRole.TEMPLATE), (switch_oligo, InputRole.STRUCTURAL),
                 (enzyme, InputRole.CATALYST)), (switched,))
    recovery = graph.input("gem_recovery", "GEM break and cDNA recovery condition", MaterialKind.CONDITION)
    recovered = graph.product("recovered", "recovered cell-barcoded cDNA", MaterialKind.MOLECULAR,
                              state=switched.state)
    graph.react("recover", "Break GEMs and recover pooled cDNA", StepKind.EXTRACTION,
                ((switched, InputRole.TEMPLATE), (recovery, InputRole.CONDITION)),
                (recovered,), changes_molecule=False)
    recovered = _cleanup(graph, reviewed, "silane_cleanup", "Post GEM-RT cleanup with Dynabeads MyOne Silane",
        recovered, "Dynabeads MyOne Silane beads",
        (("Silane-bead supernatant (leftover GEM-RT reagents and primers)",
          "removed from the magnet-held Silane beads"),))
    amplified = _source_pcr(graph, reviewed, "preamp", "Amplify recovered full-length cDNA",
                            recovered, amplified_cdna(), cdna_primers(),
                            site_template=amplified_cdna())
    amplified = _cleanup(graph, reviewed, "cdna_spri", "Post cDNA amplification SPRIselect cleanup (0.6X)",
                         amplified, "SPRIselect reagent, 0.6X", (SUPERNATANT,))
    fragmentase = graph.input("fragmentase", "CG00052 Rev D Step 3.1 Fragmentation Enzyme Blend", MaterialKind.CATALYST)
    fragments = graph.product("fragments", "fragmented barcode-bearing cDNA", MaterialKind.MOLECULAR)
    graph.react("fragment", "Fragment amplified cDNA", StepKind.FRAGMENTATION,
                ((amplified, InputRole.TEMPLATE), (fragmentase, InputRole.CATALYST)), (fragments,))
    end_mix = graph.input("end_mix", "CG00052 Rev D Step 3.1 end-repair and dA-tailing mix", MaterialKind.CATALYST)
    tailed = graph.product("tailed", "end-repaired dA-tailed cDNA fragments", MaterialKind.MOLECULAR)
    graph.react("end_prep", "End-repair and dA-tail cDNA fragments", StepKind.END_REPAIR,
                ((fragments, InputRole.TEMPLATE), (end_mix, InputRole.CATALYST)), (tailed,))
    tailed = _cleanup(graph, reviewed, "fragment_size_selection",
                      "Double-sided SPRIselect size selection (0.6X, then 0.8X)",
                      tailed, "SPRIselect reagent, 0.6X then 0.8X", (LARGE, SMALL))
    adapted = _source_ligation(graph, reviewed, tailed)
    adapted = _cleanup(graph, reviewed, "ligation_spri", "Post-ligation SPRIselect cleanup (0.8X)",
                       adapted, "SPRIselect reagent, 0.8X", (SUPERNATANT,))
    indexed = _source_pcr(graph, reviewed, "index_pcr", "Select barcode-bearing end and add sample index",
                          adapted, final, index_primers(), library=True)
    library = _cleanup(graph, reviewed, "index_size_selection",
                       "Double-sided SPRIselect size selection (0.6X, then 0.8X)",
                       indexed, "SPRIselect reagent, 0.6X then 0.8X", (LARGE, SMALL),
                       kind=MaterialKind.LIBRARY, kept_id="library",
                       kept_label="size-selected indexed v2 library")
    graph.sequence(library, SequencingPlan(library=final, primers=SEQ_PRIMERS,
        acquired_roles=tuple(p.role for p in SEQ_PRIMERS), read_lengths={"Read 1": 26, "Index 1 (i7)": 8, "Read 2": 98}))
    graph.declare_reviewed_scope(reviewed, authority)
    graph.require_publishable()
    graph.annotate(SourceDiscrepancy(
        "preamp_forward", "cDNA Primer Mix part number",
        (SourceClaim("CG000108 Rev A, Fig. 4, Protocol Step 2.2", "cDNA Primer Mix is labelled PN-22016."),
         SourceClaim("CG00052 Rev D, kit contents and Step 2.2 reaction mix",
                     "cDNA Primer Mix is listed as PN-220106.")),
        "Only the printed primer sequences are modelled; no part number is used in the construct.",
        "PN-22016 has one digit fewer than the six-digit part numbers used elsewhere in both documents.",
        "Whether the technical note omits a digit or names a different reagent lot is not stated."))
    return graph
