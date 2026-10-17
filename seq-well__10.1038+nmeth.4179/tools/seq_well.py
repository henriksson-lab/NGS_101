"""Molecular construct model for Seq-Well (Gierahn et al., 2017)."""
from __future__ import annotations

import illumina as il
import nextera as nx
import rt
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, revcomp
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (Disposition, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind, Terminal)

CELL_BARCODE = feature("cell_barcode", "cell_barcode", "random")
UMI_FEATURE = feature("umi", "umi", "random")
I7_FEATURE = feature("sample_index_i7", "sample_index", "unknown")
FEATURES = {"cell barcode": CELL_BARCODE, "UMI": UMI_FEATURE,
            "i7 reverse complement": I7_FEATURE}


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    kw.setdefault("feature", FEATURES.get(name))
    return Segment(name=name, top=top, tag=tag, **kw)


LINKER = "T" * 7
SPACER = "GCCTGTCCGCGG"
BEAD_TEMPLATE = LINKER + rt.SMART_HANDLE + "AC" + "J" * 12 + "N" * 8 + "T" * 30
TSO = rt.SMART_HANDLE + "GAAT" + "GGG"
SMART_PCR = rt.SMART_HANDLE
P5_HYBRID = il.P5 + SPACER + rt.SMART_HANDLE + "AC"
CUSTOM_R1 = SPACER + rt.SMART_HANDLE + "AC"


def bead_primer() -> list[Segment]:
    return [seg("bead linker", LINKER), seg("SMART handle", rt.SMART_HANDLE, "tso"),
            seg("AC", "AC"), seg("cell barcode", "B" * 12, "cbc", placeholder=True),
            seg("UMI", "U" * 8, "umi", placeholder=True), seg("dT30", "T" * 30)]


def tso() -> list[Segment]:
    return [seg("SMART handle", rt.SMART_HANDLE, "tso"), seg("GAAT", "GAAT"),
            seg("rGrGrG", "GGG", "tso")]


def p5_hybrid() -> list[Segment]:
    return [seg("P5", il.P5, "p5"), seg("spacer", SPACER, "r1"),
            seg("SMART handle", rt.SMART_HANDLE, "tso"), seg("AC", "AC")]


def capture_scene() -> Scene:
    mrna=[seg("RNA body","X"*28,placeholder=True),seg("poly(A)","A"*30)]
    sc=Scene(); sc.strand("mRNA",mrna,label="captured mRNA")
    sc.anneal("bead primer",bead_primer(),to="mRNA",pair=("dT30","poly(A)"),label="bead primer",mod5="bead")
    sc.arrow("bead primer","reverse transcriptase"); return sc


def first_strand() -> Construct:
    return Construct([*bead_primer(),seg("antisense cDNA","X"*28,placeholder=True),
                      seg("CCC","CCC","tso")],name="Seq-Well first-strand cDNA")


def template_switch_scene() -> Scene:
    first=first_strand(); sc=Scene(); sc.strand("first strand",list(first),label="first strand",mod5="bead")
    sc.anneal("TSO",tso(),to="first strand",pair=("rGrGrG","CCC"),label="TSO",above=True)
    sc.arrow("first strand","RT switches template"); return sc


def wta_product() -> Construct:
    return Construct([*first_strand(),seg("copied GAAT",revcomp("GAAT")),
                      seg("copied SMART handle",revcomp(rt.SMART_HANDLE),"tso")],
                     name="Seq-Well template-switched cDNA")


def selected_tagmented_strand(insert_nt: int=34) -> Construct:
    return Construct([
        seg("SMART handle",rt.SMART_HANDLE,"tso"),seg("AC","AC"),
        seg("cell barcode","B"*12,"cbc",placeholder=True),seg("UMI","U"*8,"umi",placeholder=True),
        seg("poly(T)","T"*30),seg("antisense cDNA","X"*insert_nt,placeholder=True),
        seg("mosaic end reverse complement",nx.ME_RC,"me"),seg("s7 reverse complement",nx.S7_RC,"s7")],
        name="Seq-Well bead-end tagmented strand")


READ1=sp.custom("Read 1","Seq-Well custom Read 1",CUSTOM_R1,
                "Gierahn et al. 2017, Supplementary Table 1")
SEQ_PRIMERS=(READ1,sp.NEXTERA["I1"],sp.NEXTERA["R2"])


def final_library(insert_nt: int=34) -> Construct:
    selected=selected_tagmented_strand(insert_nt)
    lib=Construct([seg("P5",il.P5,"p5"),seg("spacer",SPACER,"r1"),*list(selected),
                   seg("i7 reverse complement","I"*8,"cbc",placeholder=True),
                   seg("P7 reverse complement",il.P7_RC,"p7")],name="Seq-Well sequencing library")
    for p in SEQ_PRIMERS:
        if sp.locate(lib,p) is None: raise ValueError(f"{p.name} does not land on {lib.name}")
    return lib


def primer_landings(lib: Construct|None=None):
    lib=lib or final_library(); out=[]
    for p in SEQ_PRIMERS:
        h=sp.locate(lib,p)
        if h is None: raise ValueError(f"{p.name} does not land on {lib.name}")
        out.append((p,h))
    return out


def read1_layout() -> list[tuple[str,str]]:
    return [(x.cycles, x.feature.label) for x in
            sp.feature_spans(final_library(), SEQ_PRIMERS, {"Read 1": 20})]


MASTER = "Seq-Well Master Protocol V1.0 (2/13/2017), Nature Methods MOESM238"


def _cleanup(graph, reviewed, ident, label, template, matrix_label, kept_label, *,
             material=MaterialKind.MOLECULAR, waste_label="SPRI supernatant",
             reason="not retained by the SPRI beads"):
    matrix = graph.input(ident + "_matrix", matrix_label, MaterialKind.SELECTION)
    kept = graph.product(ident + "_kept", kept_label, material, state=template.state)
    waste = graph.product(ident + "_waste", waste_label, MaterialKind.WASTE)
    graph.react(ident, label, StepKind.CLEANUP,
                ((template, InputRole.TEMPLATE), (matrix, InputRole.SELECTION)), (kept, waste))
    graph.discard(waste, reason)
    reviewed.extend((ident + "_matrix", ident, ident + "_kept"))
    return kept


def protocol_graph() -> ProtocolGraph:
    source = ("Gierahn et al., Nat Methods 14:395 (2017), doi:10.1038/nmeth.4179; "
              "Supplementary Table 1 (MOESM237) and " + MASTER)
    graph = ProtocolGraph("Seq-Well", authority=source)
    reviewed = []
    cells = graph.input("cells", "single-cell suspension loaded onto bead-loaded arrays",
                        MaterialKind.BIOLOGICAL, primary=True)
    beads = graph.input("barcoded_beads", "barcoded capture beads (Supplementary Table 1, oligo 1)",
                        MaterialKind.ADAPTER, state=tuple(bead_primer()))
    membrane = graph.input("membrane", "plasma-treated polycarbonate membrane, clamped 37 °C 30 min",
                           MaterialKind.CONDITION)
    sealed = graph.product("sealed_wells", "membrane-sealed microwells holding cells over beads",
                           MaterialKind.MOLECULAR)
    empty_wells = graph.product("empty_wells", "unproductive cell/bead loading outcomes",
                                MaterialKind.WASTE)
    graph.react("seal_wells", "Load cells onto bead-loaded arrays and seal with the membrane",
                StepKind.PARTITION,
                ((cells, InputRole.TEMPLATE), (membrane, InputRole.CONDITION)),
                (sealed, empty_wells))
    graph.discard(empty_wells, "wells without a productive cell-bead pairing do not yield library")
    lysis = graph.input("lysis", "guanidine thiocyanate lysis buffer with Sarkosyl and "
                        "beta-mercaptoethanol, 20 min rotation", MaterialKind.CONDITION)
    lysate = graph.product("well_lysate", "membrane-retained cellular RNA in each well",
                           MaterialKind.MOLECULAR)
    graph.react("lyse", "Lyse cells through the semipermeable membrane",
                StepKind.PERMEABILIZATION,
                ((sealed, InputRole.TEMPLATE), (lysis, InputRole.CONDITION)), (lysate,))
    hyb = graph.input("hybridization", "hybridization buffer, 40 min rotation", MaterialKind.CONDITION)
    captured = graph.product("captured_rna", "microwell-confined bead-captured poly(A) RNA",
                             MaterialKind.MOLECULAR, state=capture_scene())
    graph.react("hybridize", "Hybridize released poly(A) RNA to the bead oligo-dT",
                StepKind.HYBRIDIZATION,
                ((lysate, InputRole.TEMPLATE), (beads, InputRole.STRUCTURAL),
                 (hyb, InputRole.CONDITION)), (captured,))
    recover = graph.input("bead_removal", "membrane removal, 1000 × g bead spin-off and collection",
                          MaterialKind.SELECTION)
    recovered = graph.product("recovered_beads", "recovered beads carrying hybridized mRNA",
                              MaterialKind.MOLECULAR, state=captured.state)
    array_waste = graph.product("array_waste", "array, membrane and lysate/wash buffers",
                                MaterialKind.WASTE)
    graph.react("recover_beads", "Remove the membrane and spin beads out of the array",
                StepKind.CLEANUP,
                ((captured, InputRole.TEMPLATE), (recover, InputRole.SELECTION)),
                (recovered, array_waste))
    graph.discard(array_waste, "removed with the array and wash buffer")
    tso_input = graph.input("tso", "template-switching oligo (Supplementary Table 1, oligo 2)",
                            MaterialKind.OLIGO, state=tuple(tso()))
    reverse_transcriptase = graph.input(
        "reverse_transcriptase", "Maxima H− RT in Maxima RT buffer with 4% Ficoll PM-400, "
        "RT 30 min then 52 °C 90 min", MaterialKind.CATALYST)
    first = graph.product("first_strand", first_strand().name, MaterialKind.MOLECULAR,
                          state=first_strand())
    graph.react("reverse_transcribe", "Extend the bead oligo-dT into barcoded first-strand cDNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((recovered, InputRole.TEMPLATE), (recovered, InputRole.SELF_PRIMER),
                 (reverse_transcriptase, InputRole.CATALYST)), (first,))
    switched = graph.product("template_switched", wta_product().name, MaterialKind.MOLECULAR,
                             state=wta_product())
    graph.react("template_switch", "Switch template onto the TSO in the same RT mix",
                StepKind.TEMPLATE_SWITCHING,
                ((first, InputRole.TEMPLATE), (tso_input, InputRole.STRUCTURAL),
                 (reverse_transcriptase, InputRole.CATALYST)), (switched,))
    washed = _cleanup(graph, reviewed, "rt_wash", "Wash beads in TE-TW / TE-SDS", switched,
                      "TE-TW, TE-SDS, TE-TW washes", "washed bead-bound cDNA",
                      waste_label="RT mix and free TSO", reason="removed in bead washes")
    exo = graph.input("exo1", "Exonuclease I (NEB M0293S), 37 °C 5 min then 45 min with rotation", MaterialKind.CATALYST)
    trimmed = graph.product("exo_trimmed", "bead cDNA after Exo I", MaterialKind.MOLECULAR,
                            state=wta_product())
    exo_waste = graph.product("exo_waste", "unextended single-stranded bead primers",
                              MaterialKind.WASTE)
    graph.react("exo_digest", "Exo I removes unextended bead oligos", StepKind.DIGESTION,
                ((washed, InputRole.TEMPLATE), (exo, InputRole.CATALYST)), (trimmed, exo_waste))
    graph.discard(exo_waste, "degraded by Exo I; beads then washed in TE-SDS / TE-TW")
    smart = graph.input("smart_primer", "SMART PCR primer (Supplementary Table 1, oligo 3)",
                        MaterialKind.OLIGO, state=(seg("SMART handle", SMART_PCR, "tso"),))
    preamp_pol = graph.input("preamp_polymerase",
                             "KAPA HiFi HotStart ReadyMix, 1,500–2,000 beads per 50 µL; "
                             "4 cycles (65 °C) + 9–12 cycles (67 °C)", MaterialKind.CATALYST)
    amplified = graph.product("amplified_cdna", "amplified full-length barcoded cDNA",
                              MaterialKind.MOLECULAR, state=wta_product())
    graph.react("preamp", "Whole-transcriptome amplification", StepKind.PCR,
                ((trimmed, InputRole.TEMPLATE), (smart, InputRole.STRUCTURAL),
                 (preamp_pol, InputRole.CATALYST)), (amplified,),
                geometry=PCRGeometry.symmetric(
                    wta_product(), PCRPrimer("SMART PCR primer", SMART_PCR),
                    statement="Supplementary Table 1 prints the SMART PCR primer and both copied handle sites."))
    wta_clean = _cleanup(graph, reviewed, "wta_spri", "Pool and purify WTA product (0.6× SPRI)",
                         amplified, "AMPure SPRI 0.6×, elution in 10% of the starting volume",
                         "purified WTA cDNA (1–2 kb average)")
    transposome = graph.input("transposome", "Nextera XT TD buffer + ATM, 750 pg cDNA, 55 °C 5 min, "
                              "NT neutralization", MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Nextera XT Tn5 transposase (within ATM)", MaterialKind.CATALYST)
    tagged = graph.product("selected_tagmented", selected_tagmented_strand().name,
                           MaterialKind.MOLECULAR, state=selected_tagmented_strand())
    other = graph.product("other_tagmented", "tagmented fragments without the bead-end SMART handle "
                          "(internal fragments and TSO-end fragments)", MaterialKind.MOLECULAR)
    graph.react("tagment", "Tagment the WTA cDNA with Nextera XT", StepKind.TAGMENTATION,
                ((wta_clean, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, other))
    graph.mark_terminal(Terminal(
        other.id, Disposition.UNRESOLVED,
        "Not removed by a represented separation: these fragments enter the library PCR, "
        "where the P5-SMART primer is described as selectively priming the bead-specific "
        "SMART sequence (Supplementary Table 1, note 4), so they are not exponentially "
        "amplified; the final 0.6× SPRI does not select by end identity."))
    p5 = graph.input("p5_primer", "New-P5-SMART PCR hybrid oligo (Supplementary Table 1, oligo 4)",
                     MaterialKind.OLIGO, state=tuple(p5_hybrid()))
    p7_sequence = il.P7 + "N" * 8 + nx.S7
    p7 = graph.input("p7_primer", "Nextera N70X index oligo (Illumina)", MaterialKind.OLIGO)
    library_pol = graph.input("library_polymerase", "Nextera PCR mix, 12 cycles (95/55/72 °C)",
                              MaterialKind.CATALYST)
    final = final_library()
    library = graph.product("library", final.name, MaterialKind.LIBRARY, state=final)
    graph.react("library_pcr", "Select the bead-end fragment and add P5/P7",
                StepKind.PCR,
                ((tagged, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL),
                 (p7, InputRole.STRUCTURAL), (library_pol, InputRole.CATALYST)), (library,),
                geometry=PCRGeometry.exact(
                    selected_tagmented_strand(), PCRPrimer("P5-SMART", P5_HYBRID, anneal_nt=len(SMART_PCR)+2),
                    PCRPrimer("indexed N7", p7_sequence, anneal_nt=len(nx.S7)),
                    statement="Supplementary Table 1 prints the selective P5-SMART primer; the N70X family is Illumina's Nextera index primer."))
    clean = _cleanup(graph, reviewed, "library_spri", "Purify the library (0.6× SPRI)", library,
                     "AMPure SPRI 0.6×, 13 µL elution", "Seq-Well sequencing library (650–750 bp average)",
                     material=MaterialKind.LIBRARY)
    lengths = {"Read 1": 20, "Index 1 (i7)": 8, "Read 2": 50}
    graph.sequence(clean, SequencingPlan(
        library=final, primers=SEQ_PRIMERS, acquired_roles=tuple(lengths), read_lengths=lengths,
        read_uses={"Read 1": "extract", "Index 1 (i7)": "demultiplex", "Read 2": "include"}))
    reviewed.extend((
        "cells", "barcoded_beads", "membrane", "seal_wells", "lysis", "lyse",
        "well_lysate", "hybridization", "hybridize", "bead_removal", "recover_beads", "reverse_transcriptase",
        "reverse_transcribe", "tso", "template_switch", "exo1", "exo_digest", "smart_primer",
        "preamp_polymerase", "preamp", "transposome", "tn5", "tagment", "p5_primer", "p7_primer",
        "library_polymerase", "library_pcr", "library"))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
