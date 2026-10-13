"""Shared exact chemistry for Spatial-CUT&Tag and Spatial-ATAC-seq.

The two assays differ in how Tn5 is targeted, but use the same deterministic
barcoding and sequencing-library architecture.  The oligos are from Farzad et
al., Nature Protocols 2024, Supplementary Tables 2 and 3.  The read-facing
orientation is derived here from those as-ordered oligos, so protocol pages do
not independently transcribe or reverse them.
"""
from __future__ import annotations

import illumina as il
import nextera as nx
import seqprimers as sp
from batch_ngs import seg
from chemdraw import (Construct, MolecularState, Scene, Workflow, feature,
                      revcomp)
from adapter_ligation import reported_adapter_end
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)


# As ordered, 5' -> 3'.  Asterisks in NEW_P5_SMART_WRITTEN are phosphorothioate
# bonds and therefore remain display metadata rather than bases in the model.
TN5_ME_REV_WRITTEN = "/5Phos/CTGTCTCTTATACACATCT"
TN5_ME_B = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
LIGATION_LINKER_1 = "AGTCGTACGCCGATGCTGTCTCTTATACACATCT"
LIGATION_LINKER_2 = "CGAATGCTCTGGCCTCTCAAGCACGTGGAT"
CUSTOM_READ_1 = "GCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAC"
NEW_P5_SMART = il.P5 + CUSTOM_READ_1
NEW_P5_SMART_WRITTEN = NEW_P5_SMART[:-3] + "T*A*C"
AD2_1 = "CAAGCAGAAGACGGCATACGAGATTCGCCTTAGTCTCGTGGGCTCGGAGATGTG"

# Representative members of the two 50-oligo spatial plates.
BARCODE_A_1 = "AACGTGAT"
BARCODE_B_1 = "AACGTGAT"
BARCODE_A_OLIGO_1 = "/5Phos/AGGCCAGAGCATTCG" + BARCODE_A_1 + nx.ME
BARCODE_B_OLIGO_1 = CUSTOM_READ_1 + BARCODE_B_1 + "ATCCACGTGCTTGAG"

# The sequenced strand contains the reverse complements of the ordered splints.
# Split them at the experimentally ligated bonds.  Concatenation is checked at
# import time, making the exact junctions part of construction rather than prose.
LINKER_2_READ = revcomp(LIGATION_LINKER_2)
LINKER_2_B_SIDE = BARCODE_B_OLIGO_1[-15:]
LINKER_2_A_SIDE = BARCODE_A_OLIGO_1.removeprefix("/5Phos/")[:15]
LINKER_1_READ = revcomp(LIGATION_LINKER_1)
LINKER_1_A_SIDE = BARCODE_A_OLIGO_1.removeprefix("/5Phos/")[-len(nx.ME):]
LINKER_1_GENOMIC_SIDE = LINKER_1_READ[len(nx.ME):]
if LINKER_2_B_SIDE + LINKER_2_A_SIDE != LINKER_2_READ:
    raise ValueError("spatial linker 2 oligos do not reconstruct the read strand")
if LINKER_1_A_SIDE + LINKER_1_GENOMIC_SIDE != LINKER_1_READ:
    raise ValueError("spatial linker 1 oligos do not reconstruct the read strand")


SPATIAL_A = feature(
    "spatial_a", "spatial_barcode", "whitelist", group="spatial_coordinate",
    part="A", whitelist="Farzad et al. Supplementary Table 3, first barcodes")
SPATIAL_B = feature(
    "spatial_b", "spatial_barcode", "whitelist", group="spatial_coordinate",
    part="B", whitelist="Farzad et al. Supplementary Table 3, second barcodes")
SAMPLE_I7 = feature("sample_index_i7", "sample_index", "whitelist",
                    whitelist="Farzad et al. Supplementary Table 2, Ad2 panel")


def _right_arm():
    return [
        seg("right mosaic end reverse complement", nx.ME_RC, "me"),
        seg("S7 reverse complement", nx.S7_RC, "s7"),
    ]


def tagged_fragment(insert_name: str) -> Construct:
    """Representative linker-1/S7 product after Tn5 strand transfer."""
    return Construct([
        seg("ligation linker 1", LINKER_1_READ, "me"),
        seg(insert_name, "X" * 46, placeholder=True),
        *_right_arm(),
    ], name="linker-bearing tagmented fragment")


def a_barcoded_fragment(insert_name: str) -> Construct:
    """Product after barcode A is ligated at the linker-1 junction."""
    return Construct([
        seg("ligation linker 2 A-side", LINKER_2_A_SIDE),
        seg("spatial barcode A", "A" * 8, "cbc", placeholder=True,
            feature=SPATIAL_A),
        seg("ligation linker 1 A-side", LINKER_1_A_SIDE, "me"),
        seg("ligation linker 1 genomic-side", LINKER_1_GENOMIC_SIDE, "me"),
        seg(insert_name, "X" * 46, placeholder=True),
        *_right_arm(),
    ], name="A-barcoded tagmented fragment")


def ab_barcoded_fragment(insert_name: str) -> Construct:
    """Product after barcode B is ligated inside linker 2."""
    return Construct([
        seg("custom Read 1 site", CUSTOM_READ_1, "r1"),
        seg("spatial barcode B", "B" * 8, "cbc", placeholder=True,
            feature=SPATIAL_B),
        seg("ligation linker 2 B-side", LINKER_2_B_SIDE),
        *list(a_barcoded_fragment(insert_name)),
    ], name="A/B-barcoded tagmented fragment")


def library(insert_name: str):
    """PCR-completed, single-index library with every primer site validated."""
    core = list(ab_barcoded_fragment(insert_name))
    # The intermediate already carries the S7 end; PCR adds the P5 and indexed P7
    # flow-cell arms without changing the internal barcode/linker order.
    lib = Construct([
        seg("P5", il.P5, "p5"),
        *core,
        seg("i7 reverse complement", "I" * 8, "cbc", placeholder=True,
            feature=SAMPLE_I7),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name=insert_name + " spatial epigenome library")
    primers = (
        sp.custom("Read 1", "Spatial epigenome custom Read 1", CUSTOM_READ_1,
                  "Farzad et al. 2024 Supplementary Tables 1-2"),
        sp.NEXTERA["I1"],
        sp.NEXTERA["R2"],
    )
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError("invalid spatial epigenome library: " + "; ".join(problems))
    return lib, primers


def _scene(con: Construct, label: str, junctions=()):
    sc = Scene.duplex(list(con), label=label)
    for left, right, text in junctions:
        sc.junction("top", left, right, text)
    sc.labels("top")
    return sc


def tagged_rows(insert_name: str):
    return _scene(tagged_fragment(insert_name), "tagmented fragment").rows()


def a_barcoded_rows(insert_name: str):
    return _scene(a_barcoded_fragment(insert_name), "after barcode A", (
        ("ligation linker 1 A-side", "ligation linker 1 genomic-side", "ligation"),
    )).rows()


def ab_barcoded_rows(insert_name: str):
    return _scene(ab_barcoded_fragment(insert_name), "after barcode B", (
        ("ligation linker 2 B-side", "ligation linker 2 A-side", "ligation"),
        ("ligation linker 1 A-side", "ligation linker 1 genomic-side", "ligation"),
    )).rows()


def final_rows(lib: Construct):
    return _scene(lib, "library", (
        ("ligation linker 2 B-side", "ligation linker 2 A-side", "ligation"),
        ("ligation linker 1 A-side", "ligation linker 1 genomic-side", "ligation"),
    )).rows()


def workflow(insert_name: str, *, first_action: str, first_note: str) -> Workflow:
    """Shared deterministic-barcoding workflow with assay-specific targeting."""
    start = Scene.duplex(
        [seg("fixed tissue chromatin", "X" * 46, placeholder=True)],
        label="fixed tissue chromatin")
    wf = Workflow(MolecularState("Fixed tissue chromatin", tuple(start.rows())))
    wf.react(first_action, tagged_rows(insert_name), name="Linker-bearing fragment",
             note=first_note)
    wf.react("Flow barcode A through the first microchannel array",
             a_barcoded_rows(insert_name), name="A-barcoded fragment",
             note="The first 8-nt whitelist identifies one tissue coordinate; ligation occurs inside linker 1.")
    wf.react("Rotate the device and flow barcode B through perpendicular channels",
             ab_barcoded_rows(insert_name), name="A/B-barcoded fragment",
             note="The second 8-nt whitelist identifies the orthogonal coordinate; ligation occurs inside linker 2.")
    lib, _ = library(insert_name)
    wf.react("Reverse crosslinks and PCR-complete P5, indexed P7 and flow-cell arms",
             final_rows(lib), name="Sequencing library",
             note="The ordered B+A pair identifies the channel intersection. The construct uses one i7 sample index.")
    return wf


def oligos():
    """Copyable representative oligos, exactly as written in the protocol tables."""
    return (
        f"<b>Tn5MErev</b> — 5′-{TN5_ME_REV_WRITTEN}-3′",
        f"<b>Tn5ME-B</b> — 5′-{TN5_ME_B}-3′",
        f"<b>Ligation linker 1</b> — 5′-{LIGATION_LINKER_1}-3′",
        f"<b>Ligation linker 2</b> — 5′-{LIGATION_LINKER_2}-3′",
        f"<b>Tn5ME-A_BC1_1</b> — 5′-{BARCODE_A_OLIGO_1}-3′",
        f"<b>BC2_R1_1</b> — 5′-{BARCODE_B_OLIGO_1}-3′",
        f"<b>New-P5-SMART PCR hybrid oligo</b> — 5′-{NEW_P5_SMART_WRITTEN}-3′",
        f"<b>Ad2.1_TAAGGCGA</b> — 5′-{AD2_1}-3′",
        f"<b>Custom Read 1 primer</b> — 5′-{CUSTOM_READ_1}-3′",
    )


def protocol_graph(insert_name: str, *, title: str, authority: str,
                   antibody_targeted: bool) -> ProtocolGraph:
    graph = ProtocolGraph(title, authority=authority)
    tissue = graph.input("tissue", "fixed tissue chromatin", MaterialKind.BIOLOGICAL,
                         primary=True)
    substrate = tissue
    scope = ["tissue"]
    if antibody_targeted:
        antibody = graph.input("antibody", "histone-mark primary/secondary antibodies",
                               MaterialKind.MOLECULAR)
        targeted = graph.product("targeted_chromatin", "antibody-marked tissue chromatin",
                                 MaterialKind.BIOLOGICAL)
        unbound = graph.product("unbound_antibody", "unbound antibody wash fraction",
                                MaterialKind.WASTE)
        graph.react("antibody_binding", "Bind antibodies at the selected histone mark",
                    StepKind.BINDING,
                    ((tissue, InputRole.TEMPLATE), (antibody, InputRole.ANALYTE)),
                    (targeted, unbound))
        graph.discard(unbound, "unbound antibody is removed before pA-Tn5 recruitment")
        substrate = targeted
        scope.extend(("antibody", "antibody_binding"))
    transposome = graph.input("linker_transposome", "linker-1-loaded Tn5 transposome",
                              MaterialKind.ADAPTER)
    tn5_label = "pA-Tn5 transposase" if antibody_targeted else "Tn5 transposase"
    tn5 = graph.input("tn5", tn5_label, MaterialKind.CATALYST)
    tagged = graph.product("tagged_fragment", "linker-bearing tagmented fragment",
                           MaterialKind.MOLECULAR, state=tagged_fragment(insert_name))
    t_waste = graph.product("tagmentation_waste", "non-library tagmentation products",
                            MaterialKind.WASTE)
    graph.react("tagment", "Transfer linker 1 and the S7-side adapter into chromatin",
                StepKind.TAGMENTATION,
                ((substrate, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, t_waste))
    graph.discard(t_waste, "fragments without a usable linker/S7 topology do not reach PCR")
    barcode_a = graph.input("barcode_a", "first-channel barcode-A oligo",
                            MaterialKind.ADAPTER)
    ligase = graph.input("ligase", "in-tissue DNA ligase", MaterialKind.CATALYST)
    a_product = graph.product("a_barcoded", "A-barcoded tagmented fragment",
                              MaterialKind.MOLECULAR, state=a_barcoded_fragment(insert_name))
    a_geometry = reported_adapter_end(
        "barcode-A oligo", kind="splinted", polarity="5-prime",
        terminal_sequence=LINKER_1_A_SIDE,
        statement="The published barcode-A oligo is 5-prime phosphorylated and its mosaic-end segment reconstructs the exact linker-1 ligation boundary.")
    graph.react("barcode_a_ligation", "Flow and ligate barcode A through the first channels",
                StepKind.LIGATION,
                ((tagged, InputRole.TEMPLATE), (barcode_a, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (a_product,), geometry=a_geometry)
    barcode_b = graph.input("barcode_b", "perpendicular-channel barcode-B oligo",
                            MaterialKind.ADAPTER)
    ab_product = graph.product("ab_barcoded", "A/B-barcoded spatial fragment",
                               MaterialKind.MOLECULAR, state=ab_barcoded_fragment(insert_name))
    b_geometry = reported_adapter_end(
        "barcode-B oligo", kind="splinted", polarity="5-prime",
        terminal_sequence=LINKER_2_B_SIDE,
        statement="The published barcode-B and barcode-A oligos reconstruct the exact linker-2 ligation boundary.")
    graph.react("barcode_b_ligation", "Rotate the device and ligate barcode B",
                StepKind.LIGATION,
                ((a_product, InputRole.TEMPLATE), (barcode_b, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (ab_product,), geometry=b_geometry)
    proteinase = graph.input("proteinase_k", "reverse-crosslink proteinase",
                             MaterialKind.CATALYST)
    released = graph.product("released_fragments", "reverse-crosslinked spatial fragments",
                             MaterialKind.MOLECULAR, state=ab_barcoded_fragment(insert_name))
    graph.react("reverse_crosslinks", "Reverse crosslinks and release spatial fragments",
                StepKind.DEPROTEINIZATION,
                ((ab_product, InputRole.TEMPLATE), (proteinase, InputRole.CATALYST)),
                (released,))
    p5 = graph.input("p5_primer", "New-P5-SMART PCR primer", MaterialKind.OLIGO)
    p7 = graph.input("indexed_p7_primer", "indexed Ad2 PCR primer", MaterialKind.OLIGO)
    polymerase = graph.input("pcr_polymerase", "library PCR polymerase",
                             MaterialKind.CATALYST)
    final, primers = library(insert_name)
    out = graph.product("library", final.name, MaterialKind.LIBRARY, state=final)
    graph.react("library_pcr", "Complete P5 and indexed P7 flow-cell arms",
                StepKind.PCR,
                ((released, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL),
                 (p7, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (out,), geometry=PCRGeometry.unavailable(
                    final, reason="The A/B-barcoded intermediate is exact, while PCR adds the verified outer flow-cell arms."))
    graph.sequence(out, SequencingPlan(
        library=final, primers=primers,
        acquired_roles=tuple(p.role for p in primers),
        read_lengths={"Read 1": 150, "Index 1 (i7)": 8, "Read 2": 150},
        read_uses={"Read 1": "extract_then_include",
                   "Index 1 (i7)": "demultiplex", "Read 2": "include"}))
    scope.extend(("linker_transposome", "tn5", "tagment", "barcode_a", "ligase",
                  "barcode_a_ligation", "barcode_b", "barcode_b_ligation",
                  "proteinase_k", "reverse_crosslinks", "p5_primer",
                  "indexed_p7_primer", "pcr_polymerase", "library_pcr", "library"))
    graph.declare_reviewed_scope(tuple(scope), authority)
    return graph.require_publishable()
