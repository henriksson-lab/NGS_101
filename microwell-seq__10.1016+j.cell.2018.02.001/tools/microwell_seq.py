"""Molecular construct model for Microwell-seq (Han et al., 2018)."""
from __future__ import annotations

import illumina as il
import nextera as nx
import rt
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, revcomp
from pcr import PCRGeometry
from protocol_graph import (Disposition, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind, Terminal)

UMI_FEATURE = feature("umi", "umi", "random")
I7_FEATURE = feature("sample_index_i7", "sample_index", "fixed")


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    if name == "UMI":
        kw.setdefault("feature", UMI_FEATURE)
    elif name == "i7 index read":
        kw.setdefault("feature", I7_FEATURE)
    return Segment(name=name, top=top, tag=tag, **kw)


BEAD_HEAD = "TTTAGGGATAACAGGGTAAT"
BEAD_END = "ACGT"
LINK1 = "CGACTCACTACAGGG"
LINK2 = "TCGGTGACACGATCG"
P5_SPACER = "GCCTGTCCGCGG"
I7_OLIGO_INDEX = "TAAGGCGA"
I7_INDEX_READ = revcomp(I7_OLIGO_INDEX)
TSO = rt.SMART_HANDLE + "GAAT" + "GGG"
P5_PRIMER = il.P5 + P5_SPACER + rt.SMART_HANDLE + "AC"
CUSTOM_R1 = P5_SPACER + rt.SMART_HANDLE + BEAD_END


def barcode(name: str) -> Segment:
    part = name.replace("barcode ", "round ")
    return seg(name, "B" * 6, "cbc", placeholder=True,
               feature=feature(name.replace(" ", "_"), "cell_barcode", "combinatorial",
                               group="cell_barcode", part=part))


def bead_oligo() -> Construct:
    """Published seqA bead architecture after three extension rounds."""
    return Construct([
        seg("bead head", BEAD_HEAD), seg("SMART handle", rt.SMART_HANDLE, "tso"),
        seg("seqA constant", BEAD_END), barcode("barcode 1"), seg("linker 1", LINK1),
        barcode("barcode 2"), seg("linker 2", LINK2), barcode("barcode 3"),
        seg("UMI", "N" * 6, "umi", placeholder=True), seg("poly(T)", "T" * 30),
    ], name="Microwell-seq seqA bead oligo")


def tso_oligo() -> list[Segment]:
    """Table S2 "TSO LNA": the three terminal G residues kept in source notation."""
    return [seg("SMART handle", rt.SMART_HANDLE, "tso"), seg("GAAT", "GAAT"),
            seg("source modified terminal G residues", "/rG/rG/iXNA_G", placeholder=True,
                length_bp=3)]


def barcode_oligo_sets() -> list[Construct]:
    """Table S1 split-pool oligo families A, B and C, one representative strand each.

    A is the bead-coupled round-1 strand; B and C are the extension templates, so
    their barcode windows are complements of the bead's barcode 2 and 3 and C's N6
    templates the UMI. Every variable window is an identifier placeholder.
    """
    def window(name: str, complement: bool = True) -> Segment:
        block = barcode(name)
        return Segment(name=f"{name} complement" if complement else name, top=block.top,
                       tag="cbc", placeholder=True, feature=block.feature)

    return [
        Construct([seg("bead head", BEAD_HEAD), seg("SMART handle", rt.SMART_HANDLE, "tso"),
                   seg("seqA constant", BEAD_END), window("barcode 1", complement=False),
                   seg("linker 1", LINK1)], name="Table S1 round A oligo"),
        Construct([seg("linker 2 complement", revcomp(LINK2)), window("barcode 2"),
                   seg("linker 1 complement", revcomp(LINK1))], name="Table S1 round B oligo"),
        Construct([seg("poly(A)", "A" * 30),
                   Segment(name="UMI complement", top="N" * 6, tag="umi", placeholder=True,
                           feature=UMI_FEATURE),
                   window("barcode 3"), seg("linker 2 complement", revcomp(LINK2))],
                  name="Table S1 round C oligo"),
    ]


def p7_primer(index: str | None = None) -> Construct:
    """Table S2 N7xx primer; a given index draws that row, otherwise an 8-nt window."""
    i7 = (Segment(name="i7 index", top=index, tag="cbc", feature=I7_FEATURE) if index
          else Segment(name="i7 index", top="N" * 8, tag="cbc", placeholder=True,
                       feature=I7_FEATURE))
    return Construct([seg("P7", il.P7, "p7"), i7, seg("s7", nx.S7, "s7")],
                     name="Microwell-seq N7xx library-PCR primer")


def rt_scene(insert_nt: int = 28) -> Scene:
    mrna = [seg("transcript", "X" * insert_nt, placeholder=True), seg("poly(A)", "A" * 30)]
    sc = Scene()
    sc.strand("mRNA", mrna, label="mRNA", mod3="poly(A)")
    sc.anneal("bead", list(bead_oligo()), to="mRNA", pair=("poly(T)", "poly(A)"),
              label="barcoded bead primer")
    sc.arrow("bead", "reverse transcription and template switching")
    return sc


def inferred_s7_end() -> list[Segment]:
    """The N7xx-selected s7-ME end.

    STAR Methods: a customized TruePrep V2 (Vazyme) transposase "carries two identical
    insertion sequences" and N7xx primers are used; the insertion-sequence bases are not
    printed, so s7-ME is inferred from the N7xx primer.
    """
    return [seg("mosaic end reverse complement", nx.ME_RC, "me", inferred=True),
            seg("s7 reverse complement", nx.S7_RC, "s7", inferred=True),
            seg("i7 index read", I7_INDEX_READ, "cbc", inferred=True),
            seg("P7 reverse complement", il.P7_RC, "p7", inferred=True)]


SEQ_PRIMERS = (
    sp.custom("Read 1", "Microwell-seq custom Read 1", CUSTOM_R1,
              "Han et al. 2018, Supplementary Table S2"),
    sp.NEXTERA["I1"], sp.NEXTERA["R2"],
)
RUN_ROLES = ("Read 1", "Index 1 (i7)", "Read 2")


def final_library(insert_nt: int = 28) -> Construct:
    bead = list(bead_oligo())
    lib = Construct([
        seg("P5", il.P5, "p5"), seg("custom spacer", P5_SPACER),
        *bead[1:], seg("cDNA", "X" * insert_nt, placeholder=True), *inferred_s7_end(),
    ], name="Microwell-seq sequencing library")
    problems = sp.verify(lib, SEQ_PRIMERS, required_roles=RUN_ROLES)
    if problems:
        raise ValueError("invalid final Microwell-seq construct: " + "; ".join(problems))
    return lib


def final_scene() -> Scene:
    return Scene.duplex(list(final_library()), label="library")


def read1_layout() -> list[tuple[str, str]]:
    spans = sp.feature_spans(final_library(), SEQ_PRIMERS, {"Read 1": 54})
    by_start = {x.cycle_start: (x.cycles, x.feature.label +
                (f" ({x.feature.part})" if x.feature.part else "")) for x in spans}
    return [by_start[1], ("7–21", "linker 1"), by_start[22],
            ("28–42", "linker 2"), by_start[43], by_start[49]]


def protocol_graph() -> ProtocolGraph:
    source = "Defining source: https://doi.org/10.1016/j.cell.2018.02.001"
    star = ("Han et al. 2018 Cell 172:1091, STAR Methods (Method Details: synthesis of barcoded "
            "beads, cell collection and lysis, reverse transcription, exonuclease I treatment, "
            "cDNA amplification, transposase fragmentation and selective PCR)")
    krt = "Han et al. 2018 Cell, STAR Methods Key Resources Table and Supplementary Table S2"
    inferred_source = (star + ": customized TruePrep DNA Library Prep Kit V2 (Vazyme) transposase "
                       "carrying two identical insertion sequences, N5xx replaced by the P5 "
                       "primer; the insertion-sequence bases are not printed, so the s7-ME end "
                       "is inferred from the N7xx primer")
    graph = ProtocolGraph("Microwell-seq", authority=source)
    bead_seed = graph.input("bead_seed", "seqA bead oligo seed", MaterialKind.MOLECULAR,
                            primary=True)
    barcode_oligos = graph.input("barcode_oligos", "three split-pool barcode-extension oligo sets",
                                 MaterialKind.OLIGO,
                                 source="Han et al. 2018, Supplementary Table S1 (A-, B-, C-1..96)",
                                 state=barcode_oligo_sets())
    extension_polymerase = graph.input(
        "extension_polymerase",
        "EDC coupling of amino-modified round-A oligos to carboxyl magnetic beads; rounds B/C "
        "extended with 1× Phanta Master Mix (Vazyme), 48.8 °C annealing; 95 °C strip of template strands",
        MaterialKind.CATALYST, source=star)
    bead = graph.product("barcoded_bead", bead_oligo().name, MaterialKind.MOLECULAR,
                         state=bead_oligo())
    graph.react("bead_assembly", "Build the bead barcode by three split-pool extension rounds",
                StepKind.EXTENSION,
                ((bead_seed, InputRole.TEMPLATE), (barcode_oligos, InputRole.STRUCTURAL),
                 (extension_polymerase, InputRole.CATALYST)), (bead,))
    rna = graph.input("rna", "RNA of single cells in agarose microwells, released by cold lysis buffer "
                      "(0.1 M Tris-HCl pH 7.5, 0.5 M LiCl, 1% SDS, 10 mM EDTA, 5 mM DTT; 12 min) "
                      "onto the co-loaded bead; beads then washed in 6× SSC and Tris",
                      MaterialKind.BIOLOGICAL, primary=True, source=star)
    reverse_transcriptase = graph.input(
        "reverse_transcriptase",
        "SuperScript II RT (200 U), 1 M betaine, 6 mM MgCl2, 2.5 mM DTT, 1 mM dNTP, 42 °C 90 min "
        "on pooled beads; TE-SDS wash inactivates the RT",
        MaterialKind.CATALYST, source=star)
    first_strand = graph.product("first_strand", "bead-primed first-strand cDNA",
                                 MaterialKind.MOLECULAR, state=rt_scene())
    graph.react("reverse_transcription", "Capture microwell RNA and reverse-transcribe",
                StepKind.REVERSE_TRANSCRIPTION,
                ((rna, InputRole.TEMPLATE), (bead, InputRole.STRUCTURAL),
                 (reverse_transcriptase, InputRole.CATALYST)), (first_strand,))
    # Han et al. 2018 Supplementary Table S2 "TSO LNA" (also the Guo-lab protocol
    # Table 1 TSO_RT primer); the SMART handle assignment rests on these printed rows.
    tso = graph.input("tso", "published template-switching oligo", MaterialKind.OLIGO,
                      source=krt + " (TSO LNA / 'Primer for Reverse transcription')",
                      state=Construct(tso_oligo(), name="Microwell-seq TSO"))
    switched = graph.product("switched_cdna", "SMART-handle-flanked cDNA",
                             MaterialKind.MOLECULAR, state=rt_scene())
    graph.react("template_switch", "Template-switch at the cDNA 5′ end",
                StepKind.TEMPLATE_SWITCHING,
                ((first_strand, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL),
                 (reverse_transcriptase, InputRole.CATALYST)), (switched,))
    exonuclease = graph.input("exonuclease_i", "Exonuclease I (NEB) 1 U/µl, 37 °C 60 min",
                              MaterialKind.CATALYST, source=star)
    trimmed = graph.product("exo_treated_cdna", "bead cDNA after Exo I",
                            MaterialKind.MOLECULAR, state=rt_scene())
    unused = graph.product("unused_bead_oligos", "digested bead oligos that captured no mRNA",
                           MaterialKind.WASTE, source=star)
    graph.react("exo_i", "Digest bead oligos that captured no mRNA", StepKind.DIGESTION,
                ((switched, InputRole.TEMPLATE), (exonuclease, InputRole.CATALYST)),
                (trimmed, unused), source=star)
    graph.discard(unused, "removed by the TE-SDS / TE-TW / Tris bead washes")
    smart_primer = graph.input("smart_primer", "single TSO-PCR primer", MaterialKind.OLIGO,
                               source=krt + " (TSO-PCR / 'Primer for cDNA amplification')",
                               state=Construct([seg("TSO_PCR primer", rt.SMART_HANDLE)]))
    polymerase = graph.input("pcr_polymerase",
                             "KAPA HiFi HotStart ReadyMix (WTA: 4 cycles at 65 °C + 10–14 cycles at 67 °C, then AMPure XP; "
                             "library PCR: 72 °C 3 min, 14 cycles at 60 °C)",
                             MaterialKind.CATALYST, source=star)
    wta = graph.product("wta_cdna", "whole-transcriptome-amplified cDNA",
                        MaterialKind.MOLECULAR)
    graph.react("wta_pcr", "Amplify with the single SMART-handle primer", StepKind.PCR,
                ((trimmed, InputRole.TEMPLATE), (smart_primer, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (wta,),
                geometry=PCRGeometry.unavailable(
                    final_library(), reason="The WTA substrate is symmetric before the final tagmentation shell.",
                    symmetric_species=True))
    tn5_adapter = graph.input("tn5_adapter", "identical Tn5 insertion sequence at both ends (s7-ME inferred from N7xx)",
                              MaterialKind.ADAPTER, inferred=True, source=inferred_source)
    tn5 = graph.input("tn5", "customized transposase, TruePrep DNA Library Prep Kit V2 for Illumina (Vazyme)",
                      MaterialKind.CATALYST, source=star)
    tagged = graph.product("tagged_cdna", "tagmented WTA fragments", MaterialKind.MOLECULAR,
                           inferred=True, source=inferred_source)
    graph.react("tagmentation", "INFERRED — Tagment WTA with a single-insertion-sequence TruePrep V2 transposase; the insertion-sequence bases are not printed (s7-ME inferred)",
                StepKind.TAGMENTATION,
                ((wta, InputRole.TEMPLATE), (tn5_adapter, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged,), source=inferred_source)
    p5 = graph.input("p5_primer", "published bead-end-selective P5 primer (3′ *A*C)",
                     MaterialKind.OLIGO, source="Han et al. 2018, Supplementary Table S2 P5",
                     state=Construct([seg("P5 primer", P5_PRIMER)]))
    p7 = graph.input("p7_primer", "published N7xx indexed P7 primer (N701 drawn)",
                     MaterialKind.OLIGO, source="Han et al. 2018, Supplementary Table S2 N7xx",
                     state=p7_primer())
    final = final_library()
    library = graph.product("library", final.name, MaterialKind.LIBRARY, state=final,
                            inferred=True, source=inferred_source)
    looped = graph.product("self_looped_fragments",
                           "internal cDNA fragments with identical Tn5 ends (self-loop, not amplified)",
                           MaterialKind.MOLECULAR, source=star)
    graph.react("library_pcr", "INFERRED — Select the bead end with P5 + N7xx and PCR-complete the library; AMPure XP removes primer dimers and large fragments (the s7/ME and P7-side regions follow the printed primers)", StepKind.PCR,
                ((tagged, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL),
                 (p7, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)), (library, looped),
                geometry=PCRGeometry.unavailable(
                    final, reason="The final construct is reconstructed from exact PCR primers around an unpublished Tn5 adapter."),
                source=inferred_source)
    graph.mark_terminal(Terminal(looped.id, Disposition.UNRESOLVED,
        "STAR Methods: fragments without the bead end 'form self-loops, which impede their "
        "binding to PCR primers'; they are not exponentially amplified but no separation "
        "specific to them is described"))
    graph.sequence(library, SequencingPlan(
        library=final, primers=SEQ_PRIMERS, acquired_roles=RUN_ROLES,
        read_uses={"Read 1": "extract", "Index 1 (i7)": "demultiplex",
                   "Read 2": "include"},
        limitation=("The s7/ME and P7-side construct regions are inferred from the source "
                    "boundary. Acquired cycle counts are not declared: the 2018 STAR Methods "
                    "say only 'sequencing on an Illumina Hiseq system' and Tables S1/S2 print "
                    "none (the GEO GSE108097 samples name HiSeq 2500; the undated lab web "
                    "protocol names NextSeq 500). Read 1 needs at least 54 cycles to reach the UMI end and the "
                    "Table S2 i7 indices are 8 nt; these are informative minima computed "
                    "from the oligos, not published run lengths.")))
    graph.declare_reviewed_scope((
        "bead_seed", "barcode_oligos", "extension_polymerase", "bead_assembly",
        "rna", "reverse_transcriptase", "reverse_transcription", "tso",
        "template_switch", "exonuclease_i", "exo_i", "smart_primer", "pcr_polymerase", "wta_pcr",
        "tn5_adapter", "tn5", "tagmentation", "p5_primer", "p7_primer",
        "library_pcr", "library", "self_looped_fragments"), source)
    return graph.require_publishable()
