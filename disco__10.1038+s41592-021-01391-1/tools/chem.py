"""Protocol-owned source, molecular workflow and reaction graph for DisCo.

Source scope (bounded pass 2026-10-17). The defining Nature Methods article
(doi:10.1038/s41592-021-01391-1) is subscription-gated and its Methods remain
unread. The library stages, reagents and printed oligos below are taken from the
authors' "DisCo step-by-step laboratory protocol" (EPFL Deplancke lab, PDF dated
14 Sep 2021): p. 2 "Required primers"; "Reverse transcription"; "Exonuclease I
treatment and bead elution"; "cDNA amplification"; "Tagmentation of cDNA with
in-house Tn5"; "Purification of the tagmented library"; "Sequencing Strategy".
That protocol describes itself as the authors' most up-to-date protocol; its
equivalence to the published experiments is not established here.

Not in the inspected originals: the capture-bead oligo sequence (cell-barcode and
UMI lengths, constant bases), the Tn5 p7|ME adapter strands, the i7 index-primer
bases (cited only as Buenrostro et al. 2013) and the Read 2 / index primers.
Those stay explicit roles; PCR geometry is recorded as unavailable.
"""
from chemdraw import Construct, Segment, feature
import illumina as il
import rt
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import chain, molecule, role
from extended_protocols import Definition
import seqprimers as sp

DEFINITION = Definition('DisCo', '10.1038/s41592-021-01391-1', 'Bues et al., Nature Methods (2022)', 'Machine vision and multilayer valves deterministically co-encapsulate one low-input cell and one mRNA-capture bead; captured RNA is reverse-transcribed with template switching on barcoded beads, amplified by SMART PCR and tagmented into a 3′ RNA-seq library.', 'Read 1 (custom Drop-seq SeqB primer, 20 cycles) reports the bead-derived cell barcode and UMI; Read 2 (50 cycles) reports cDNA from the Tn5 end.', '01_disco.html')

AUTHOR_PROTOCOL = ("Bues, Biocanin, Pezoldt et al., DisCo step-by-step laboratory protocol "
                   "(EPFL Deplancke lab, PDF 14 Sep 2021; "
                   "https://www.epfl.ch/labs/deplanckelab/wp-content/uploads/2021/09/DisCo_protocol.pdf)")

# --- Oligos as printed in the author protocol, p. 2 "Required primers" ---------------
# Source transcriptions, kept with their written modifications.
TSO_AS_WRITTEN = "AAGCAGTGGTATCAACGCAGAGTGAATrGrGrG"
SMART_PCR_AS_WRITTEN = "AAGCAGTGGTATCAACGCAGAGT"
P5_SMART_AS_WRITTEN = "AATGATACGGCGACCACCGAGATCTACACGCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGT*A*C"
SEQB_AS_WRITTEN = "GCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAC"

# Decomposition into canonical parts (computed and checked below, never retyped
# into a sequence column).
SPACER = "GCCTGTCCGCGG"          # between P5 and the SMART handle in P5 SMART PCR / SeqB
PRIMER_3P = "AC"                 # 3'-terminal bases of P5 SMART PCR (written *A*C) and SeqB
TSO_LINKER = "GAAT"


def _letters(written: str) -> str:
    return written.replace("r", "").replace("*", "")


def _decompose() -> None:
    """Fail construction if the printed oligos do not decompose as modelled."""
    expected = {
        "TSO": (TSO_AS_WRITTEN, rt.SMART_HANDLE + TSO_LINKER + "GGG"),
        "SMART PCR": (SMART_PCR_AS_WRITTEN, rt.SMART_HANDLE),
        "P5 SMART PCR": (P5_SMART_AS_WRITTEN, il.P5 + SPACER + rt.SMART_HANDLE + PRIMER_3P),
        "SeqB": (SEQB_AS_WRITTEN, SPACER + rt.SMART_HANDLE + PRIMER_3P),
    }
    for name, (written, parts) in expected.items():
        if _letters(written) != parts:
            raise ValueError(f"{name} printed sequence does not decompose as modelled")


_decompose()

CELL_BARCODE = feature("cell_barcode", "cell_barcode", "unknown")
UMI_FEATURE = feature("umi", "umi", "unknown")
I7_FEATURE = feature("sample_index_i7", "sample_index", "unknown")


def tso() -> Construct:
    return Construct([Segment("SMART handle", rt.SMART_HANDLE, "tso"),
                      Segment("GAAT", TSO_LINKER),
                      Segment("rGrGrG (RNA)", "GGG", "tso")], name="Template Switch Oligo")


def smart_pcr_primer() -> Construct:
    return Construct([Segment("SMART handle", rt.SMART_HANDLE, "tso")], name="SMART PCR primer")


def p5_smart_pcr() -> Construct:
    return Construct([Segment("P5", il.P5, "p5"), Segment("spacer", SPACER, "r1"),
                      Segment("SMART handle", rt.SMART_HANDLE, "tso"),
                      Segment("A*C (two phosphorothioates)", PRIMER_3P)],
                     name="P5 SMART PCR")


SEQB = SPACER + rt.SMART_HANDLE + PRIMER_3P


def seq_primers():
    return (sp.custom("Read 1", "Custom Drop-seq SeqB sequencing primer", SEQB,
                      AUTHOR_PROTOCOL + ", p. 2 Required primers and Sequencing Strategy"),)


READ_LENGTHS = {"Read 1": 20}


def final_library() -> Construct:
    """Selected bead-end product of the P5 SMART PCR / i7 index-primer PCR.

    The left end is the P5 SMART PCR primer (exact, printed). Its 3' AC must be
    followed by bead-derived barcode/UMI/poly(T) sequence; the bead oligo itself is
    not in the inspected originals, so those regions are roles. The right end is the
    Tn5 p7|ME adapter extended by a Buenrostro i7 index primer; no bases are printed.
    """
    lib = Construct([
        Segment("P5", il.P5, "p5"), Segment("spacer", SPACER, "r1"),
        Segment("SMART handle", rt.SMART_HANDLE, "tso"), Segment("primer 3′ AC", PRIMER_3P),
        role("cell barcode", tag="cbc", molecular_feature=CELL_BARCODE),
        role("UMI", tag="umi", molecular_feature=UMI_FEATURE),
        role("bead oligo(dT)-derived poly(T)"),
        role("3-prime cDNA insert"),
        role("Tn5 mosaic end + p7 adapter (bases not printed)"),
        role("i7 sample index", tag="cbc", molecular_feature=I7_FEATURE),
        role("P7 (from Buenrostro i7 index primer)"),
    ], name="DisCo bead-end library")
    if sp.locate(lib, seq_primers()[0]) is None:
        raise ValueError("SeqB does not land on the DisCo library")
    return lib


TITLE = DEFINITION.title
NOTES = DEFINITION.notes
SOURCE = (f'Defining source: <a href="https://doi.org/{DEFINITION.doi}">{DEFINITION.citation}</a> '
          f'(full text not read: subscription-gated). Library chemistry from the authors\' '
          f'<a href="https://www.epfl.ch/labs/deplanckelab/wp-content/uploads/2021/09/DisCo_protocol.pdf">'
          f'DisCo step-by-step laboratory protocol</a> (EPFL, Sep 2021).')
SUMMARY = DEFINITION.summary
CAVEAT = ("<b>Source scope.</b> The reaction steps, printed oligos and read lengths come from the "
          "authors' step-by-step protocol, which describes itself as their most up-to-date version; "
          "the defining Nature Methods Methods section was not available, so equivalence with the "
          "published experiments (for example cDNA and library PCR cycle numbers) is not established. "
          "The capture-bead oligo, the Tn5 p7|ME adapter, the i7 index primer and the Read 2/index "
          "primers are not printed in the inspected sources; they are drawn as named roles, and the "
          "two PCR geometries are recorded as unavailable.")
FINAL_LIBRARY = final_library()
SEQ_PRIMERS = seq_primers()
FINAL_CAPTION = ("Selected bead-end product: P5 SMART PCR primer (printed) — bead-derived cell barcode, "
                 "UMI and poly(T) (bases not printed) — cDNA — Tn5 p7|ME end and i7 index (bases not printed).")
SEQUENCING_INTRO = ("NextSeq, following Drop-seq specifications (author protocol, Sequencing Strategy): "
                    "Read 1 = 20 bp primed by the custom Drop-seq SeqB primer; Read 2 = 50 bp. "
                    "The SeqB site is computed on the printed P5 SMART PCR end.")
SEQUENCING_ENDING = DEFINITION.endpoint
SEQUENCING_UNAVAILABLE = ("Read 2 (50 cycles) and the i7 index read are named but their primers, the "
                          "Tn5 p7|ME adapter bases and the i7 index-primer bases are not printed in the "
                          "inspected sources; only the Read 1 SeqB landing is computed.")


def oligos():
    def row(name, written, note):
        return f"<b>{name}</b>: 5′-{written}-3′ <small>({note})</small>"
    return (
        row("Template Switch Oligo", TSO_AS_WRITTEN, "rGrGrG = three 3′ riboguanosines; 2.5 µM final in RT"),
        row("SMART PCR primer", SMART_PCR_AS_WRITTEN, "cDNA amplification"),
        row("P5 SMART PCR", P5_SMART_AS_WRITTEN, "* = phosphorothioate; library PCR"),
        row("Custom Drop-seq SeqB sequencing primer", SEQB_AS_WRITTEN, "Read 1"),
        "<small>All four as printed on p. 2 (Required primers) of the author protocol. The i7 index "
        "primers are cited there only as Buenrostro et al., Nat. Methods 10, 1213 (2013); no bases "
        "are printed.</small>",
    )


def workflow():
    return chain(molecule("low-input cell suspension and mRNA-capture beads", (role("single cells"), role("barcoded capture bead"))), (
        ("Machine-vision position exactly one cell and one bead at separate valves",
         molecule("deterministically paired cell and bead", (role("one cell"), role("one bead"))),
         "Active positioning replaces Poisson co-loading for low-input samples."),
        ("Co-eject the pair into one lysis droplet and route it to the sample outlet",
         molecule("occupied DisCo droplet", (role("lysed cell RNA"), role("bead oligo(dT)"))),
         "Droplets are collected under stuffer droplets and transferred to the bead-capture chip."),
        ("RT with template switching on the captured beads, Exo I, SMART PCR",
         molecule("amplified full-length cDNA", (role("SMART handle"), role("cell barcode"), role("UMI"),
                                                 role("poly(T)"), role("cDNA"), role("SMART handle"))),
         "Maxima H− RT with the TSO; Exonuclease I; SMART PCR on eluted beads."),
        ("Tagment with p7|ME-loaded Tn5 and amplify with P5 SMART PCR + i7 index primer",
         molecule("DisCo bead-end library", tuple(FINAL_LIBRARY.segments)),
         "Only bead-end fragments carry the SMART handle followed by AC at the P5 side."),
    ))


def protocol_graph() -> ProtocolGraph:
    d = DEFINITION
    final = FINAL_LIBRARY
    article = f"{d.citation}; doi:{d.doi}"
    proto = AUTHOR_PROTOCOL
    graph = ProtocolGraph(d.title, authority=article)
    scope: dict[str, str] = {}

    def from_protocol(*ids):
        for x in ids:
            scope["src_" + x] = x

    inputs = graph.input(
        "cells_and_beads", "low-input cells and barcoded mRNA-capture beads in lysis buffer",
        MaterialKind.BIOLOGICAL,
        state=molecule("low-input suspension and capture beads", (
            role("single cells"), role("barcoded capture bead"))), primary=True, source=proto)
    positioning = graph.input("vision_positioning", "machine-vision valve positioning condition",
                              MaterialKind.CONDITION)
    paired = graph.product("paired_cell_bead", "deterministically paired cell and bead",
                           MaterialKind.MOLECULAR,
                           state=molecule("deterministically paired cell and bead", (
                               role("one cell"), role("one bead"))))
    rejected = graph.product("positioning_rejects", "empty or mis-paired valve events",
                             MaterialKind.WASTE)
    graph.react("position_pair", "Position exactly one cell and one bead",
                StepKind.PARTITION,
                ((inputs, InputRole.TEMPLATE), (positioning, InputRole.CONDITION)),
                (paired, rejected))
    graph.discard(rejected, "fails deterministic one-cell/one-bead positioning")
    droplet_condition = graph.input("droplet_condition", "co-ejection and droplet-sorting condition",
                                    MaterialKind.CONDITION)
    occupied = graph.product("occupied_droplet", "occupied DisCo droplet with lysed cell RNA",
                             MaterialKind.MOLECULAR,
                             state=molecule("occupied DisCo droplet", (
                                 role("lysed cell RNA"), role("bead oligo(dT)"))))
    failed_droplets = graph.product("failed_droplets", "droplets routed to the waste outlet",
                                    MaterialKind.WASTE)
    graph.react("make_droplet", "Co-eject the pair and sort the droplet to the sample outlet",
                StepKind.PARTITION,
                ((paired, InputRole.TEMPLATE), (droplet_condition, InputRole.CONDITION)),
                (occupied, failed_droplets))
    graph.discard(failed_droplets, "does not reach the sample outlet")
    scope.update({"src_cells_and_beads": "cells_and_beads",
                  "src_vision_positioning": "vision_positioning",
                  "src_position_pair": "position_pair",
                  "src_droplet_condition": "droplet_condition",
                  "src_make_droplet": "make_droplet"})

    # Bead capture chip: droplets transferred, beads retained on pillars, washed.
    cp_chip = graph.input("cp_chip_wash", "bead-capture (cp) chip pillars; 2x 6X SSC and 1x RT-buffer washes",
                          MaterialKind.SELECTION, source=proto)
    captured = graph.product("captured_beads", "beads with hybridized poly(A) RNA retained on the cp-chip",
                             MaterialKind.MOLECULAR, source=proto,
                             state=molecule("captured beads", (role("bead oligo(dT)"), role("poly(A) RNA"))))
    emulsion_waste = graph.product("emulsion_waste", "oil, stuffer droplets and wash flow-through",
                                   MaterialKind.WASTE, source=proto)
    graph.react("capture_beads", "Recover beads on the cp-chip and wash (6X SSC, RT buffer)",
                StepKind.CLEANUP,
                ((occupied, InputRole.TEMPLATE), (cp_chip, InputRole.SELECTION)),
                (captured, emulsion_waste), source=proto)
    graph.discard(emulsion_waste, "flows through the bead-capture chip")

    bead_primer = graph.input("bead_oligo_dt", "bead oligo(dT) primer carrying cell barcode and UMI (sequence not in inspected sources)",
                              MaterialKind.OLIGO, source=proto,
                              state=molecule("bead oligo(dT)", (role("bead PCR handle (not printed)", inferred=True), role("cell barcode"), role("UMI"), role("oligo(dT)"))))
    rtase = graph.input("maxima_h", "Maxima H− reverse transcriptase (Maxima 5x RT buffer, Ficoll PM-400, dNTPs, RNase inhibitor)",
                        MaterialKind.CATALYST, source=proto)
    first = graph.product("first_strand", "bead-bound first-strand cDNA with 3′ non-templated C overhang",
                          MaterialKind.MOLECULAR, source=proto, inferred=True,
                          state=molecule("first-strand cDNA", (role("bead oligo(dT)"), role("antisense cDNA"), role("CCC", inferred=True))))
    graph.react("rt", "Reverse-transcribe on the cp-chip, 42 °C 90 min",
                StepKind.REVERSE_TRANSCRIPTION,
                ((captured, InputRole.TEMPLATE), (bead_primer, InputRole.STRUCTURAL),
                 (rtase, InputRole.CATALYST)), (first,), source=proto,
                reason="RT mix contains the TSO; the C-overhang intermediate is the standard template-switching inference, not printed.")
    tso_in = graph.input("tso", "Template Switch Oligo (rGrGrG)", MaterialKind.OLIGO,
                         state=tso(), source=proto + ", p. 2")
    switched = graph.product("switched_cdna", "bead-bound cDNA ending in the copied TSO",
                             MaterialKind.MOLECULAR, source=proto,
                             state=molecule("template-switched cDNA", (
                                 role("bead oligo(dT)"), role("antisense cDNA"), role("TSO complement"))))
    graph.react("template_switch", "Template switch onto the TSO (same RT reaction)",
                StepKind.TEMPLATE_SWITCHING,
                ((first, InputRole.TEMPLATE), (tso_in, InputRole.STRUCTURAL),
                 (rtase, InputRole.CATALYST)), (switched,), source=proto)
    rt_wash = graph.input("rt_wash", "TE-SDS once, TE-TW twice, 10 mM Tris once", MaterialKind.CONDITION, source=proto)
    washed = graph.product("washed_beads", "washed beads carrying cDNA", MaterialKind.MOLECULAR,
                           state=switched.state, source=proto)
    rt_waste = graph.product("rt_waste", "RT mix including free TSO, RNA and enzyme", MaterialKind.WASTE, source=proto)
    graph.react("rt_wash_step", "Wash beads on chip after RT", StepKind.CLEANUP,
                ((switched, InputRole.TEMPLATE), (rt_wash, InputRole.CONDITION)), (washed, rt_waste), source=proto)
    graph.discard(rt_waste, "washed off the cp-chip")

    exo = graph.input("exo_i", "Exonuclease I (10X Exo I buffer), 37 °C 45 min", MaterialKind.CATALYST, source=proto)
    exo_beads = graph.product("exo_beads", "beads with cDNA after Exo I", MaterialKind.MOLECULAR,
                              state=switched.state, source=proto)
    digested = graph.product("exo_digest", "nucleotides from unextended single-stranded bead oligos (expected Exo I substrate; purpose not stated in protocol)",
                             MaterialKind.WASTE, source=proto, inferred=True)
    graph.react("exo_digest_step", "Exonuclease I treatment on chip", StepKind.DIGESTION,
                ((washed, InputRole.TEMPLATE), (exo, InputRole.CATALYST)), (exo_beads, digested), source=proto)
    graph.discard(digested, "removed by the TE-SDS / 2x TE-TW / water washes that follow Exo I on the chip")
    elution = graph.input("bead_elution", "TE-SDS, 2x TE-TW, water washes; elute beads from cp-chip with water; pellet",
                          MaterialKind.CONDITION, source=proto)
    eluted = graph.product("eluted_beads", "eluted cDNA beads in PCR tube", MaterialKind.MOLECULAR,
                           state=switched.state, source=proto)
    exo_wash = graph.product("exo_wash", "Exo I mix and wash supernatant", MaterialKind.WASTE, source=proto)
    graph.react("elute", "Wash and elute beads from the cp-chip", StepKind.CLEANUP,
                ((exo_beads, InputRole.TEMPLATE), (elution, InputRole.CONDITION)), (eluted, exo_wash), source=proto)
    graph.discard(exo_wash, "supernatant removed before PCR")

    smart = graph.input("smart_pcr_primer", "SMART PCR primer", MaterialKind.OLIGO,
                        state=smart_pcr_primer(), source=proto + ", p. 2")
    kapa_rm = graph.input("kapa_readymix", "2x KAPA HiFi HotStart ReadyMix", MaterialKind.CATALYST, source=proto)
    amplified = graph.product("amplified_cdna", "SMART-amplified cDNA with beads", MaterialKind.MOLECULAR, source=proto,
                              state=molecule("amplified cDNA", (role("bead PCR handle (SMART-primed; inferred)", inferred=True), role("cell barcode"), role("UMI"),
                                                                role("poly(T)"), role("cDNA"), role("SMART handle"))))
    graph.react("cdna_pcr", "SMART PCR: 4 cycles (65 °C anneal) + 15–19 cycles (67 °C)", StepKind.PCR,
                ((eluted, InputRole.TEMPLATE), (smart, InputRole.STRUCTURAL), (kapa_rm, InputRole.CATALYST)),
                (amplified,), source=proto,
                geometry=PCRGeometry.unavailable(
                    final, symmetric_species=True,
                    reason="Single SMART PCR primer; the bead-end handle is not printed in the inspected sources, so its site is not computed."))
    bead_removal = graph.input("bead_removal", "spin down capture beads and transfer the PCR supernatant",
                               MaterialKind.CONDITION, source=proto)
    soluble = graph.product("soluble_cdna", "bead-free amplified cDNA", MaterialKind.MOLECULAR,
                            state=amplified.state, source=proto)
    spent_beads = graph.product("spent_beads", "capture beads retaining bead-bound strands", MaterialKind.WASTE, source=proto)
    graph.react("remove_beads", "Separate the PCR supernatant from the capture beads", StepKind.CLEANUP,
                ((amplified, InputRole.TEMPLATE), (bead_removal, InputRole.CONDITION)), (soluble, spent_beads), source=proto)
    graph.discard(spent_beads, "left behind when the supernatant is transferred")
    ampure1 = graph.input("ampure_06x", "0.6x AMPure beads", MaterialKind.SELECTION, source=proto)
    cdna = graph.product("purified_cdna", "purified cDNA (Qubit)", MaterialKind.MOLECULAR,
                         state=amplified.state, source=proto)
    small1 = graph.product("small_cdna_waste", "fragments below ~200 bp and primers", MaterialKind.WASTE, source=proto)
    graph.react("cdna_cleanup", "0.6x AMPure cleanup", StepKind.CLEANUP,
                ((soluble, InputRole.TEMPLATE), (ampure1, InputRole.SELECTION)), (cdna, small1), source=proto)
    graph.discard(small1, "left in the 0.6x supernatant")

    tn5_adapter = graph.input("tn5_p7_me", "Tn5 loaded with p7|ME (bases not printed)", MaterialKind.ADAPTER, source=proto)
    tn5 = graph.input("tn5", "in-house Tn5 (Picelli et al. 2014), 5x TAPS, 55 °C 6 min, 2 ng cDNA", MaterialKind.CATALYST, source=proto)
    tagged = graph.product("tagmented", "tagmented cDNA fragments carrying p7|ME ends", MaterialKind.MOLECULAR, source=proto,
                           state=molecule("tagmented cDNA", (role("bead PCR handle (inferred)", inferred=True), role("cell barcode"), role("UMI"),
                                                            role("poly(T)"), role("cDNA"), role("Tn5 p7|ME end"))))
    graph.react("tagment", "Tagment cDNA with p7|ME-loaded Tn5", StepKind.TAGMENTATION,
                ((cdna, InputRole.TEMPLATE), (tn5_adapter, InputRole.STRUCTURAL), (tn5, InputRole.CATALYST)),
                (tagged,), source=proto,
                reason="Drawn state is the bead-end fragment selected by the later P5 SMART PCR; TSO-end and internal fragments also form, their fate is not stated.")
    sds = graph.input("sds_stop", "1.5 µL 0.2% SDS, RT 5 min", MaterialKind.CONDITION, source=proto)
    stopped = graph.product("tagmented_stopped", "SDS-stopped tagmentation", MaterialKind.MOLECULAR,
                            state=tagged.state, source=proto)
    graph.react("stop_tn5", "Stop tagmentation with SDS", StepKind.QUENCHING,
                ((tagged, InputRole.TEMPLATE), (sds, InputRole.CONDITION)), (stopped,), source=proto)

    p5s = graph.input("p5_smart_pcr", "P5 SMART PCR (two 3′ phosphorothioates)", MaterialKind.OLIGO,
                      state=p5_smart_pcr(), source=proto + ", p. 2")
    i7 = graph.input("i7_index_primer", "i7 index primer (Buenrostro et al. 2013 set; bases not printed)",
                     MaterialKind.OLIGO, source=proto)
    kapa = graph.input("kapa_hifi", "KAPA HiFi polymerase, 5x buffer, dNTP", MaterialKind.CATALYST, source=proto)
    pcr_product = graph.product("indexed_pcr", "indexed bead-end library PCR product", MaterialKind.MOLECULAR,
                                state=final, source=proto)
    graph.react("library_pcr", "72 °C 3 min, 98 °C 30 s, then 16 cycles (63 °C anneal) with P5 SMART PCR + i7 index primer",
                StepKind.PCR,
                ((stopped, InputRole.TEMPLATE), (p5s, InputRole.STRUCTURAL), (i7, InputRole.STRUCTURAL),
                 (kapa, InputRole.CATALYST)), (pcr_product,), source=proto,
                geometry=PCRGeometry.unavailable(
                    final,
                    reason="P5 SMART PCR is printed, but the bead-end template, Tn5 p7|ME adapter and i7 index-primer bases are not printed in the inspected sources."))
    ampure2 = graph.input("ampure_07x", "0.7x AMPure XP", MaterialKind.SELECTION, source=proto)
    library = graph.product("library", "DisCo sequencing library", MaterialKind.LIBRARY, state=final, source=proto)
    small2 = graph.product("small_lib_waste", "short products and primers", MaterialKind.WASTE, source=proto)
    graph.react("library_cleanup", "0.7x AMPure XP cleanup; BioAnalyzer/Fragment Analyzer", StepKind.CLEANUP,
                ((pcr_product, InputRole.TEMPLATE), (ampure2, InputRole.SELECTION)), (library, small2), source=proto)
    graph.discard(small2, "left in the 0.7x supernatant")
    graph.sequence(library, SequencingPlan(
        library=final, primers=SEQ_PRIMERS, acquired_roles=("Read 1",), read_lengths=READ_LENGTHS,
        limitation=SEQUENCING_UNAVAILABLE))

    from_protocol("cp_chip_wash", "capture_beads", "bead_oligo_dt", "maxima_h", "rt", "tso",
                  "template_switch", "rt_wash", "rt_wash_step", "exo_i", "exo_digest_step",
                  "bead_elution", "elute", "smart_pcr_primer", "kapa_readymix", "cdna_pcr",
                  "bead_removal", "remove_beads", "ampure_06x", "cdna_cleanup", "tn5_p7_me", "tn5",
                  "tagment", "sds_stop", "stop_tn5", "p5_smart_pcr", "i7_index_primer", "kapa_hifi",
                  "library_pcr", "ampure_07x", "library_cleanup", "library")
    graph.declare_reviewed_scope(scope, proto + "; positioning/droplet stages also " + article)
    return graph.require_publishable()
