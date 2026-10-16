"""Molecular construct model for the paper-described PIP-seq workflow."""
from __future__ import annotations

import illumina as il
import nextera as nx
import rt
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene, Segment, feature
from pcr import PCRGeometry
from protocol_graph import (Disposition, Terminal, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind)

CELL_CASSETTE = feature("cell_barcode_cassette", "cell_barcode", "combinatorial",
                        note="four barcode blocks; internal boundaries unpublished")
UMI_FEATURE = feature("umi", "umi", "random")
I7_FEATURE = feature("sample_index_i7", "sample_index", "fixed")
FEATURES = {"ligated barcode/linker cassette": CELL_CASSETTE, "UMI": UMI_FEATURE,
            "i7 index read": I7_FEATURE}


def bead_ligation_geometry():
    return AdapterLigationGeometry.unavailable(
        "The defining PIP-seq paper states that the hydrogel-bead barcode cassette is assembled by split-pool ligation, but it does not publish the barcode-block ends or splint oligos.")


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    kw.setdefault("feature", FEATURES.get(name))
    return Segment(name=name, top=top, tag=tag, **kw)


BEAD_PREFIX = "T" * 7
BEAD_JUNCTION = "ACGACTC"
BEAD_READ1_PART = il.TRUSEQ_READ1[3:28]
P5_SPACER = "GCCTGTCCGCGG"
TSO = rt.SMART_HANDLE + "GAAT" + "GGG"
P5_LIBRARY = il.P5 + P5_SPACER + rt.SMART_HANDLE + "AC"
I7_INDEX_READ = "TAAGGCGA"


def inferred_barcode_cassette() -> Segment:
    """Unpublished split-pool block: keep it visibly unresolved for every caller."""
    return seg("ligated barcode/linker cassette", "X" * 36, "cbc",
               placeholder=True, inferred=True,
               note="The paper states four ligated barcodes but does not print their blocks or splints.")


def inferred_read1_completion() -> Segment:
    """Unprinted ligation bases needed to complete the named Read-1 primer site."""
    return seg("Read-1-site completion", il.TRUSEQ_READ1[28:], "r1", inferred=True,
               note="These five bases are absent from the printed bead primer.")


def bead_oligo() -> Construct:
    return Construct([
        seg("bead spacer", BEAD_PREFIX), seg("SMART handle", rt.SMART_HANDLE, "tso"),
        seg("bead junction", BEAD_JUNCTION),
        seg("partial TruSeq Read 1", BEAD_READ1_PART, "r1"),
        inferred_read1_completion(), inferred_barcode_cassette(),
        seg("UMI", "N" * 12, "umi", placeholder=True),
        seg("poly(T)", "T" * 19), seg("V anchor", "V", placeholder=True),
    ], name="PIP-seq bead oligo")


def tso_oligo() -> list[Segment]:
    return [seg("SMART handle", rt.SMART_HANDLE, "tso"), seg("GAAT", "GAAT"),
            seg("rGrGrG", "rGrGrG", placeholder=True)]


def rt_scene(insert_nt: int = 28) -> Scene:
    mrna = [seg("transcript", "X" * (insert_nt - 1), placeholder=True),
            seg("anchor partner", "V", placeholder=True), seg("poly(A)", "A" * 19)]
    sc = Scene()
    sc.strand("mRNA", mrna, label="mRNA", mod3="poly(A)")
    sc.anneal("bead", list(bead_oligo()), to="mRNA", pair=("poly(T)", "poly(A)"),
              label="released bead primer")
    sc.arrow("bead", "bulk reverse transcription and template switching")
    return sc


SEQ_PRIMERS = (
    sp.custom("Read 1", "TruSeq Read 1", il.TRUSEQ_READ1,
              "Illumina; the library omits the site's first three bases",
              "Its 5' three bases form an unpaired flap; the 3' 30 nt land exactly."),
    sp.NEXTERA["I1"], sp.NEXTERA["R2"],
)
RUN_ROLES = ("Read 1", "Index 1 (i7)", "Read 2")


def final_library(insert_nt: int = 28) -> Construct:
    bead = list(bead_oligo())
    lib = Construct([
        seg("P5", il.P5, "p5"), seg("P5 custom spacer", P5_SPACER),
        *bead[1:], seg("cDNA", "X" * insert_nt, placeholder=True),
        seg("mosaic end reverse complement", nx.ME_RC, "me"),
        seg("s7 reverse complement", nx.S7_RC, "s7"),
        seg("i7 index read", I7_INDEX_READ, "cbc"),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="paper-described PIP-seq library")
    problems = sp.verify(lib, SEQ_PRIMERS, required_roles=RUN_ROLES)
    if problems:
        raise ValueError("invalid PIP-seq library: " + "; ".join(problems))
    return lib


def final_scene() -> Scene:
    return Scene.duplex(list(final_library()), label="library")


def protocol_graph() -> ProtocolGraph:
    source = "Defining source: https://doi.org/10.1038/s41587-023-01685-z"
    unresolved = source + "; barcode-block and splint sequences are not published"
    graph = ProtocolGraph("PIP-seq", authority=source)
    bead_seed = graph.input("bead_seed", "hydrogel bead primer seed",
                            MaterialKind.MOLECULAR, primary=True)
    barcode_cassettes = graph.input("barcode_cassettes",
                                    "four split-pool barcode/linker cassette rounds",
                                    MaterialKind.ADAPTER, inferred=True, source=unresolved)
    ligase = graph.input("ligase", "bead barcode ligase", MaterialKind.CATALYST)
    bead = graph.product("barcoded_bead", bead_oligo().name, MaterialKind.MOLECULAR,
                         state=bead_oligo(), inferred=True, source=unresolved)
    graph.react("bead_ligation", "Build split-pool-barcoded hydrogel beads",
                StepKind.LIGATION,
                ((bead_seed, InputRole.TEMPLATE), (barcode_cassettes, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (bead,), geometry=bead_ligation_geometry(),
                source=unresolved)
    cells = graph.input("cells", "single-cell suspension", MaterialKind.BIOLOGICAL,
                        primary=True)
    loaded = graph.product("loaded_droplets", "cell-and-bead particle-templated droplets",
                           MaterialKind.BIOLOGICAL)
    excluded = graph.product("partition_loss", "empty or multiply loaded droplets",
                             MaterialKind.WASTE)
    graph.react("partition", "Vortex cells and beads into particle-templated droplets",
                StepKind.PARTITION,
                ((cells, InputRole.TEMPLATE), (bead, InputRole.STRUCTURAL)),
                (loaded, excluded), changes_molecule=False)
    graph.discard(excluded, "droplets outside the usable cell-and-bead occupancy class are excluded")
    lysis = graph.input("lysis_condition", "droplet lysis and emulsion-break condition",
                        MaterialKind.CONDITION)
    rna = graph.product("released_rna", "released cellular mRNA", MaterialKind.MOLECULAR)
    primer = graph.product("released_bead_primer", "released bead capture primer",
                           MaterialKind.OLIGO, state=bead_oligo(), inferred=True,
                           source=unresolved)
    graph.react("lysis", "Lyse cells at65°C for35min, then break the emulsion", StepKind.EXTRACTION,
                ((loaded, InputRole.TEMPLATE), (lysis, InputRole.CONDITION)), (rna, primer))
    reverse_transcriptase = graph.input("reverse_transcriptase",
                                       "template-switching reverse transcriptase",
                                       MaterialKind.CATALYST)
    first_strand = graph.product("first_strand", "bead-primed cDNA", MaterialKind.MOLECULAR,
                                 state=rt_scene(), inferred=True, source=unresolved)
    graph.react("reverse_transcription", "Reverse-transcribe captured mRNA in bulk",
                StepKind.REVERSE_TRANSCRIPTION,
                ((rna, InputRole.TEMPLATE), (primer, InputRole.STRUCTURAL),
                 (reverse_transcriptase, InputRole.CATALYST)), (first_strand,))
    tso = graph.input("tso", "PIPS template-switching oligo", MaterialKind.OLIGO)
    switched = graph.product("switched_cdna", "SMART-handle-flanked cDNA",
                             MaterialKind.MOLECULAR, state=rt_scene(), inferred=True,
                             source=unresolved)
    graph.react("template_switch", "Template-switch at the cDNA 5′ end",
                StepKind.TEMPLATE_SWITCHING,
                ((first_strand, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL),
                 (reverse_transcriptase, InputRole.CATALYST)), (switched,))
    wta_primer = graph.input("wta_primer", "single PIPS WTA primer", MaterialKind.OLIGO)
    polymerase = graph.input("polymerase", "WTA and library PCR polymerase",
                             MaterialKind.CATALYST)
    wta = graph.product("wta_cdna", "whole-transcriptome-amplified cDNA",
                        MaterialKind.MOLECULAR)
    graph.react("wta_pcr", "Amplify handle-flanked cDNA: 16 WTA cycles", StepKind.PCR,
                ((switched, InputRole.TEMPLATE), (wta_primer, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (wta,),
                geometry=PCRGeometry.unavailable(
                    final_library(), reason="The WTA substrate is symmetric before tagmentation.",
                    symmetric_species=True))
    tn5_adapter = graph.input("tn5_adapter", "s7 tagmentation adapter", MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagged_cdna", "s7-tagmented WTA fragments", MaterialKind.MOLECULAR)
    tag_waste = graph.product("tagmentation_waste", "non-library tagmentation products",
                              MaterialKind.WASTE)
    graph.react("tagmentation", "Tagment the amplified cDNA", StepKind.TAGMENTATION,
                ((wta, InputRole.TEMPLATE), (tn5_adapter, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, tag_waste))
    graph.discard(tag_waste, "fragments without a usable bead-end/s7 pair do not amplify")
    p5 = graph.input("p5_primer", "published bead-end-selective P5 primer", MaterialKind.OLIGO)
    p7 = graph.input("p7_primer", "named N7xx indexed P7 primer", MaterialKind.OLIGO)
    final = final_library()
    library = graph.product("library", final.name, MaterialKind.LIBRARY, state=final,
                            inferred=True, source=unresolved)
    graph.react("library_pcr", "Select the bead end and PCR-complete the library", StepKind.PCR,
                ((tagged, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL),
                 (p7, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)), (library,),
                geometry=PCRGeometry.unavailable(
                    final, reason="Unpublished barcode/linker bases prevent an exact full PCR landing model."),
                source=unresolved)
    graph.sequence(library, SequencingPlan(
        library=final, primers=SEQ_PRIMERS, acquired_roles=RUN_ROLES,
        read_uses={"Read 1": "extract", "Index 1 (i7)": "demultiplex",
                   "Read 2": "include"},
        limitation="Exact cell-barcode cycle boundaries remain unavailable because the four barcode blocks and splints are unpublished."))
    graph.declare_reviewed_scope((
        "bead_seed", "barcode_cassettes", "ligase", "bead_ligation", "cells",
        "partition", "lysis_condition", "lysis", "reverse_transcriptase",
        "reverse_transcription", "tso", "template_switch", "wta_primer",
        "polymerase", "wta_pcr", "tn5_adapter", "tn5", "tagmentation",
        "p5_primer", "p7_primer", "library_pcr", "library"), source)
    return graph.require_publishable()
