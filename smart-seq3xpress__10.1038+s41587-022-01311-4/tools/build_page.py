#!/usr/bin/env python3
"""Build the SMART-seq3xpress protocol schematic."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "lib"))
from protocol_graph import ProtocolGraph  # noqa: E402

OUT = HERE.parent / "smart-seq3xpress.html"

# The original xpress Methods specify 10-base dual indexes. The protocol's
# attached 384set_Nextera_BC_order.xls supplies complete oligos in two sheets.
# A1 on each sheet is an illustrative pair, not a recovered experimental well.
import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Segment, Scene, feature, oligo, panel, protocol_graph_panels, revcomp
from page import head, info
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (Disposition, Terminal, InputRole, MaterialKind, SequencingPlan, StepKind,
                            AnnotationSlot, SourceClaim, SourceDiscrepancy)
from protocol_views import assemble_page

I5_A1 = "TTGGTACGCG"
I7_A1 = "ATAGGCGCTC"
INDEX_FORWARD = il.P5 + I5_A1 + nx.S5
INDEX_REVERSE = il.P7 + I7_A1 + nx.S7
SEQ_PRIMERS = tuple(sp.NEXTERA[k] for k in ("R1", "I1", "I2", "R2"))
TAG = "ATTGCGCAATG"
DT_HANDLE = "ACGAGCATCAGCAGCATACGA"
TSO_PREFIX = "AGAGACAG"
DT_LENGTH = 30
UMI_RANDOM_LENGTH = 8
SPACER = "WW"
FORWARD_PCR = nx.S5 + nx.ME + TAG
REVERSE_PCR = DT_HANDLE
SOURCE = "Final2022 Methods and author protocol DOI10.17504/protocols.io.yxmvmk1yng3p/v1"
UMI = feature("umi", "umi", "random")
I5 = feature("sample_index_i5", "sample_index", "fixed")
I7 = feature("sample_index_i7", "sample_index", "fixed")


def tso_segments():
    """Final published improved TSO; RNA terminal Gs drawn in DNA alphabet."""
    return [Segment("ME3", TSO_PREFIX, "me"), Segment("5-prime tag", TAG, "r1"),
            Segment("Eight-N UMI region", "N" * UMI_RANDOM_LENGTH, "umi", placeholder=True, feature=UMI),
            Segment("WW spacer / UMI positions 9-10", SPACER, "r1", placeholder=True, feature=UMI),
            Segment("rGrGrG", "GGG", "tso")]


def rt_primer():
    return Construct([Segment("dT handle", DT_HANDLE, "r3"),
                      Segment("dT", "T" * DT_LENGTH),
                      Segment("VN anchor", "VN", placeholder=True)],
                     name="Biotinylated Smartseq3_OligodT30VN")


def switched_cdna(*, amplified=False):
    # The forward primer adds s5 and the rest of ME during amplification;
    # those bases are not already present in the ordered TSO.
    prefix = ([Segment("forward-primer handle", FORWARD_PCR, "s5"),
               *tso_segments()[2:]] if amplified else tso_segments())
    return Construct([*prefix, Segment("transcript", "X" * 32, placeholder=True),
                      Segment("anchor complement", "NB", placeholder=True),
                      Segment("polyA", "A" * DT_LENGTH),
                      Segment("dT handle reverse complement", revcomp(DT_HANDLE), "r3")],
                     name="Amplified xpress cDNA" if amplified else "Template-switched xpress cDNA")


def require_preamp_sites(product):
    return PCRGeometry.exact(product,
        PCRPrimer("Forward PCR", FORWARD_PCR, anneal_nt=len(FORWARD_PCR)),
        PCRPrimer("Reverse PCR", REVERSE_PCR, anneal_nt=len(REVERSE_PCR)),
        statement="Complete source-known primers face inward on the amplified product, including the introduced forward-primer overhang; this site check does not prove the complete input-to-product transition.")


def ordered_oligos():
    rows = [oligo("Smartseq3_OligodT30VN", rt_primer().segments, mods="/5BiosG/"),
            oligo("Smartseq3xpress_TSO", [*tso_segments()[:-1],
                  Segment("RNA G tail", "rGrGrG", placeholder=True)], mods="/5BiosG/"),
            oligo("Fwd_PCR_primer", [Segment("", FORWARD_PCR[:-3]),
                  Segment("terminal phosphorothioates", "*".join(FORWARD_PCR[-3:]), placeholder=True)]),
            oligo("Rev_PCR_primer", [Segment("", REVERSE_PCR[:-3]),
                  Segment("terminal phosphorothioates", "*".join(REVERSE_PCR[-3:]), placeholder=True)]),
            oligo("Illustrative i5 A1", [Segment("", INDEX_FORWARD)]),
            oligo("Illustrative i7 A1", [Segment("", INDEX_REVERSE)])]
    return "<h2>Adapter and primer sequences</h2><seq>" + "".join(rows) + "</seq>"


def final_library(*, five_prime=True):
    middle = ([Segment("5-prime tag", TAG, "r1"),
               Segment("Eight-N UMI region", "N" * UMI_RANDOM_LENGTH, "umi", placeholder=True, feature=UMI),
               Segment("WW spacer / UMI positions 9-10", SPACER, "r1", placeholder=True, feature=UMI),
               Segment("template-switch GGG", "GGG", "tso")]
              if five_prime else [])
    library = Construct([
        Segment("P5", il.P5, "p5"),
        Segment("illustrative i5 A1", I5_A1, "cbc", feature=I5),
        Segment("s5", nx.S5, "s5"), Segment("ME", nx.ME, "me"),
        *middle, Segment("cDNA", "X" * 32, placeholder=True),
        Segment("ME reverse complement", nx.ME_RC, "me"),
        Segment("s7 reverse complement", nx.S7_RC, "s7"),
        Segment("illustrative i7 A1 reverse complement", revcomp(I7_A1), "cbc", feature=I7),
        Segment("P7 reverse complement", il.P7_RC, "p7")],
        name="SMART-seq3xpress " + ("5-prime" if five_prime else "internal") + " library; illustrative A1/A1")
    require_index_sites(library)
    problems = sp.verify(library, SEQ_PRIMERS,
                         required_roles=tuple(p.role for p in SEQ_PRIMERS))
    if problems:
        raise ValueError("; ".join(problems))
    return library


def require_index_sites(library):
    return PCRGeometry.exact(library,
        PCRPrimer("UD_i5_indices A1", INDEX_FORWARD, anneal_nt=len(INDEX_FORWARD)),
        PCRPrimer("UD_i7_indices A1", INDEX_REVERSE, anneal_nt=len(INDEX_REVERSE)),
        statement="Complete original workbook A1 primers match the displayed illustrative 10-base dual-index endpoint.")


def index_pcr(graph, filled, polymerase):
    forward = graph.input("index_forward", "Custom i5 A1 primer (illustrative)", MaterialKind.OLIGO,
                          state=Construct([Segment("ordered i5 A1", INDEX_FORWARD)]))
    reverse = graph.input("index_reverse", "Custom i7 A1 primer (illustrative)", MaterialKind.OLIGO,
                          state=Construct([Segment("ordered i7 A1", INDEX_REVERSE)]))
    five, internal = final_library(), final_library(five_prime=False)
    library = graph.product("indexed_library", "5-prime library, illustrative A1/A1", MaterialKind.LIBRARY, state=five)
    other = graph.product("indexed_internal", "Internal library, illustrative A1/A1", MaterialKind.LIBRARY, state=internal)
    graph.react("index_pcr", "Index PCR, 10-base dual indexes (illustrative A1/A1): 95 C 30 s, "
                "12-14x(95 C 10 s, 55 C 30 s, 72 C 30-60 s), 72 C 5 min", StepKind.PCR,
                ((filled, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library, other), geometry=require_index_sites(five))
    # Final Methods: plates spun out into a reservoir; pool purified with 22% PEG beads 1:0.7.
    beads = graph.input("library_beads", "Homemade 22% PEG bead suspension, 0.7:1", MaterialKind.SELECTION)
    pooled = graph.product("library", "Pooled, bead-purified 5-prime library", MaterialKind.LIBRARY, state=five)
    pooled_internal = graph.product("internal_library", "Pooled, bead-purified internal library",
                                    MaterialKind.LIBRARY, state=internal)
    supernatant = graph.product("library_supernatant", "Bead supernatant and washes", MaterialKind.WASTE)
    graph.react("library_cleanup", "Spin-out pooling of plates and 0.7:1 PEG-bead cleanup", StepKind.CLEANUP,
                ((library, InputRole.TEMPLATE), (other, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                (pooled, pooled_internal, supernatant), changes_molecule=False)
    graph.discard(supernatant, "Removed on the magnet; removal of primers and short products is not quantified.")
    for material, construct in ((pooled, five), (pooled_internal, internal)):
        graph.sequence(material, SequencingPlan(library=construct, primers=SEQ_PRIMERS,
                       acquired_roles=tuple(p.role for p in SEQ_PRIMERS)))


def render_page():
    graph = protocol_graph()
    ordered = ordered_oligos() + info("Terminal phosphorothioate linkages are marked *. Both A1 index oligos come from the author-linked workbook; their pairing is illustrative and does not assign a barcode to a published cell.")
    finals, sequencing = [], []
    for five_prime in (True, False):
        con = final_library(five_prime=five_prime)
        heading = "5-prime fragments" if five_prime else "Internal fragments"
        finals.extend(["<h3>" + heading + "</h3>", panel(Scene.duplex(con.segments).rows(), cls="long")])
        sequencing.append(sp.section(con, SEQ_PRIMERS, heading="Sequencing: " + heading,
                                     show_read_content=True))
    return assemble_page(graph, [head("SMART-seq3xpress library chemistry"),
        '<div class="wrap" data-protocol-complete="true"><h1>SMART-seq3xpress — full-length single-cell mRNA</h1>',
        info('Defining source: <a href="https://doi.org/10.1038/s41587-022-01311-4">Hagemann-Jensen et al.2022</a>; original <a href="https://doi.org/10.17504/protocols.io.yxmvmk1yng3p/v1">author protocol v1</a>.'),
        info("This declaration selects the final published improved WW TSO, Maxima H Minus RT, the SeqAmp pre-amplification option, dilution of the cDNA without bead cleanup, the commercial ATM (Nextera XT) or TDE1 tagmentation route with SDS stop, Phusion index PCR, and the pooled 0.7:1 PEG-bead cleanup. The paper also tests other pre-amplification enzymes and in-house Tn5; their choices are not silently inherited. Tagmented-end and extension intermediates remain role-level states here: their exact covalent geometry is not established by drawing the final library. Stock Nextera sequencing-primer sites are computed on that endpoint; no experiment-specific run cycle recipe is selected."),
        info("This page represents the Illumina-compatible dual-index library before any MGI conversion. The defining paper also reports MGI adapter conversion/circularization and custom sequencing primers; that branch is outside this schematic."),
        AnnotationSlot.INTRO, ordered,
        "<h2>Step-by-step library generation</h2>", protocol_graph_panels(graph, cls="long"),
        AnnotationSlot.AFTER_WORKFLOW, "<h2>Final libraries</h2>", *finals,
        "<h2>Library sequencing</h2>", *sequencing,
        info("Index source: author protocol DOI yxmvmk1yng3p/v1, attached 384set_Nextera_BC_order.xls, sheets UD_i5_indices and UD_i7_indices, A1. The workbook defines available oligos; it does not assign this illustrative pair to an experimental well."),
        info("The ordered TSO has eight N positions followed by the two-base WW spacer. The deposited annotation tables call NNNNNNNNWW a ten-position UMI, matching UMI(12-21) after the eleven-base tag in the author protocol. Both regions therefore carry the UMI feature here; only eight positions are fully degenerate N, and WW retains its spacer function."),
        '<details class="sources"><summary>Sources for the declaration</summary><ul>'
        '<li>Final paper Methods: low-volume RT, improved TSO, pre-amplification options, library preparation and sequencing.</li>'
        '<li>Original author protocol v1: modified ordered oligos and linked complete index workbook.</li>'
        '<li>Workbook sheets UD_i5_indices and UD_i7_indices: complete illustrative A1 primers.</li>'
        '</ul></details>', "</div>"])


def protocol_graph() -> ProtocolGraph:
    graph = ProtocolGraph("SMART-seq3xpress", authority=SOURCE)
    mrna = graph.input("mrna", "single-cell poly(A) RNA", MaterialKind.MOLECULAR,
        state=Construct([Segment("transcript", "X" * 32, placeholder=True),
                         Segment("polyA", "A" * DT_LENGTH)]), primary=True)
    primer = graph.input("rt_primer", "Biotinylated anchored oligo(dT)30VN", MaterialKind.OLIGO,
                         state=rt_primer())
    enzyme = graph.input("reverse_transcriptase", "Maxima H Minus reverse transcriptase (2 U/uL)", MaterialKind.CATALYST)
    # Final Methods, "Low-volume Smart-seq3 and Smart-seq3xpress library preparation".
    lysis = graph.input("lysis", "Lysis under 3 uL overlay: 0.1% Triton X-100, 5% PEG8000, 0.5 mM dNTP each, "
                        "RNase inhibitor; 72 C 10 min", MaterialKind.CONDITION)
    rt_buffer = graph.input("rt_buffer", "RT mix: 25 mM Tris-HCl pH 8.3, 30 mM NaCl, 0.5 mM GTP, 2.5 mM MgCl2, "
                            "8 mM DTT; 42 C 90 min, 10x(50 C 2 min, 42 C 2 min)", MaterialKind.CONDITION)
    first = graph.product("first_strand", "First-strand cDNA; terminal addition not resolved base by base", MaterialKind.MOLECULAR)
    graph.react("reverse_transcription", "Reverse-transcribe from anchored oligo(dT)", StepKind.REVERSE_TRANSCRIPTION,
                ((mrna, InputRole.TEMPLATE), (primer, InputRole.STRUCTURAL), (lysis, InputRole.CONDITION),
                 (rt_buffer, InputRole.CONDITION), (enzyme, InputRole.CATALYST)), (first,))
    tso = graph.input("tso", "Biotinylated improved WW template-switch oligo", MaterialKind.OLIGO,
                      state=Construct(tso_segments()))
    switched = graph.product("switched", "Template-switched full-length cDNA", MaterialKind.MOLECULAR,
                            state=switched_cdna())
    graph.react("template_switch", "Template-switch using the final published improved TSO", StepKind.TEMPLATE_SWITCHING,
                ((first, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL), (enzyme, InputRole.CATALYST)), (switched,))
    forward = graph.input("preamp_forward", "Phosphorothioate-protected forward PCR primer", MaterialKind.OLIGO,
                          state=Construct([Segment("forward PCR", FORWARD_PCR)]))
    reverse = graph.input("preamp_reverse", "Phosphorothioate-protected reverse PCR primer", MaterialKind.OLIGO,
                          state=Construct([Segment("reverse PCR", REVERSE_PCR)]))
    polymerase = graph.input("preamp_polymerase", "SeqAmp polymerase, 0.025 U/uL in 0.6 uL (represented published option); 95 C 1 min, "
                             "6-18 cycles of 98 C 10 s / 65 C 30 s / 68 C 2-6 min, 72 C 10 min", MaterialKind.CATALYST)
    preamp = switched_cdna(amplified=True)
    amplified = graph.product("amplified", "Amplified cDNA including the forward-primer s5/ME overhang", MaterialKind.MOLECULAR, state=preamp)
    graph.react("preamp", "Pre-amplify with the published asymmetric primer pair", StepKind.PCR,
                ((switched, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)), (amplified,),
                geometry=require_preamp_sites(preamp))
    diluted = graph.product("diluted_cdna", "Diluted cDNA; 1 uL of the well transferred", MaterialKind.MOLECULAR,
                            state=preamp)
    water = graph.input("dilution_water", "9 uL water added per well (no bead cleanup)", MaterialKind.CONDITION)
    remainder = graph.product("cdna_remainder", "Diluted cDNA remaining in the pre-amplification well",
                              MaterialKind.MOLECULAR, state=preamp)
    graph.react("dilute", "Dilute pre-amplified cDNA and transfer 1 uL into tagmentation", StepKind.SAMPLING,
                ((amplified, InputRole.TEMPLATE), (water, InputRole.CONDITION)), (diluted, remainder),
                changes_molecule=False)
    graph.mark_terminal(Terminal(remainder.id, Disposition.STORED,
                                 "Not carried into this library; its later use is not represented."))
    transposome = graph.input("transposome", "Commercial ATM (Nextera XT) or TDE1 loaded transposome (published options)", MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Tn5 activity in the commercial transposome", MaterialKind.CATALYST)
    td = graph.input("td_buffer", "1x tagmentation buffer: 10 mM Tris pH 7.5, 5 mM MgCl2, 5% DMF; 55 C 10 min",
                     MaterialKind.CONDITION)
    tagged = graph.product("tagmented", "Tagmented cDNA fragments; exact end geometry unassigned", MaterialKind.MOLECULAR)
    other = graph.product("tagment_waste", "Other tagmentation products outside the represented endpoint", MaterialKind.WASTE)
    graph.react("tagment", "Tagment amplified cDNA using the commercial Nextera option", StepKind.TAGMENTATION,
                ((diluted, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL), (td, InputRole.CONDITION),
                 (tn5, InputRole.CATALYST)), (tagged, other))
    graph.discard(other, "Outside the represented final-library classes; no absolute amplification exclusion or product frequency is asserted.")
    sds = graph.input("sds", "0.5 uL 0.2% SDS per well", MaterialKind.CONDITION)
    released = graph.product("tn5_released", "Tn5-stripped tagmented fragments; exact end geometry unassigned",
                             MaterialKind.MOLECULAR)
    graph.react("strip_tn5", "Stop tagmentation and strip Tn5 with SDS", StepKind.DEPROTEINIZATION,
                ((tagged, InputRole.TEMPLATE), (sds, InputRole.CONDITION)), (released,), changes_molecule=False)
    index_polymerase = graph.input("index_polymerase", "Phusion DNA polymerase, 0.01 U/uL in the index-PCR mix", MaterialKind.CATALYST)
    filled = graph.product("gap_filled", "Extension intermediate before index PCR; exact end state unassigned", MaterialKind.MOLECULAR)
    graph.react("gap_fill", "Initial 72 C 3 min extension hold before index PCR", StepKind.EXTENSION,
                ((released, InputRole.TEMPLATE), (index_polymerase, InputRole.CATALYST)), (filled,))
    index_pcr(graph, filled, index_polymerase)
    graph.declare_reviewed_scope(("mrna", "rt_primer", "reverse_transcriptase", "reverse_transcription",
        "tso", "template_switch", "preamp_forward", "preamp_reverse", "preamp_polymerase", "preamp",
        "transposome", "tn5", "tagment", "gap_fill", "index_forward", "index_reverse",
        "index_polymerase", "index_pcr", "library", "internal_library", "lysis", "rt_buffer",
        "diluted_cdna", "cdna_remainder", "dilution_water", "dilute", "td_buffer", "sds", "tn5_released", "strip_tn5",
        "indexed_library", "indexed_internal", "library_beads", "library_cleanup"), SOURCE)
    graph.annotate(SourceDiscrepancy(
        "index_pcr", "Index length and the published processing example",
        (SourceClaim("Author protocol v1, p33, bcl2fastq example",
                     "The example uses Y150N,I8,I8,Y150N, selecting eight cycles from each index read."),
         SourceClaim("Author protocol v1, p34, configuration; final paper library-preparation Methods; linked order workbook",
                     "The configuration requests BC(1-10) from both index files; the article specifies ten-base dual indexes and the workbook supplies ten-base index primers.")),
        "The molecular endpoint uses source-known ten-base indexes, with A1/A1 explicitly illustrative; no run mask is selected.",
        "Complete ordered primers and the final article independently support the ten-base molecular architecture.",
        "The eight-cycle example does not supply ten index bases for its adjacent configuration. The sources do not establish which example or run mask applied to every experiment; the page does not silently repair the command."))
    graph.annotate(SourceDiscrepancy(
        "preamp", "SeqAmp primer concentrations",
        (SourceClaim("Final paper Methods, SeqAmp paragraph",
                     "PCR mastermix contains '0.5 uM/1 uM Smartseq3 forward and reverse primer'."),
         SourceClaim("Author protocol v1, step 9 PCR table (p14)",
                     "Forward and reverse primers are each listed at 0.5 uM reaction concentration.")),
        "The graph names the asymmetric primer pair and does not select a concentration.",
        "Both sources agree on primer identity, SeqAmp polymerase at 0.025 U/uL and a 0.6 uL mix.",
        "The article's '0.5 uM/1 uM' may denote forward/reverse values or tested alternatives; the sources do not establish which concentration each experiment used."))
    return graph.require_publishable()


def main() -> None:
    OUT.write_text(render_page(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
