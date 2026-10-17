"""CEL-Seq2 source-specific primers, read recipe and graph invocation."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] /
                       "cel-seq__10.1016+j.celrep.2012.08.003" / "tools"))

import celseq as C
import seqprimers as sp
import illumina as il
from chemdraw import Construct, Segment, oligo, panel, protocol_graph_panels, revcomp, strand_row
from page import head, info, table
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            SourceClaim, SourceDiscrepancy, StepKind)

PROTOCOL = "CEL-Seq2"
SOURCE = ("10.1186/s13059-016-0938-8", "Genome Biology", "2016")
DESIGN = C.PrimerDesign(
    spacer="GCCGG", ra5="AGTTCTACAGTCCGACGATC",
    umi_nt=6, barcode_nt=6,
    randomhex_tail="GCCTTGGCACCCGAGAATTCCA", randomhex_nt=6,
)
# The author protocol specifies Illumina RP1 and RPIX but refers the reader to
# Illumina for their bases. Illumina Adapter Sequences, doc 1000000002694 v22,
# "TruSeq Small RNA" prints both oligos and the RPI1 index. The shortened
# CEL-Seq2 RA5 tract leaves a 16-base RP1 3' landing site; RPI1 has a 21-base
# 3' landing site on the randomhexRT tail (source text and geometry below).
VENDOR_SOURCE = "https://support-docs.illumina.com/SHARE/AdapterSequences/Content/TruSeq-SmallRNA.htm"
RP1 = "AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA"
RPI1_INDEX = "CGTGAT"
RPI1 = il.P7 + RPI1_INDEX + "GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
SOURCE_OLIGO_EXCEPTIONS = {
    "pcr_reverse": ("Illumina doc 1000000002694 v22 prints the RPI family with a "
                    "[6 bases] slot and RPI1 CGTGAT in its index table, rather than "
                    "printing the assembled RPI1 oligo on one source line."),
}


def library_pcr_template() -> Construct:
    """Representative plus-strand sequence after randomhexRT; insert stays unknown."""
    return Construct([
        Segment("T7-start/RA5 tag", C.T7_TRANSCRIBED_G + DESIGN.ra5),
        Segment("UMI", "N" * DESIGN.umi_nt, "umi", placeholder=True,
                feature=C.UMI_FEATURE),
        Segment("cell barcode", "N" * DESIGN.barcode_nt, "cbc",
                placeholder=True, feature=C.CELL_BARCODE),
        Segment("poly(T)/cDNA", "TTTXXXXXXXX...XXXXXXXX", placeholder=True),
        Segment("randomhexRT tail reverse complement", revcomp(DESIGN.randomhex_tail)),
    ], name="CEL-Seq2 random-primed library PCR template")


def library_pcr_geometry(template: Construct | None = None) -> PCRGeometry:
    template = template or library_pcr_template()
    return PCRGeometry.exact(
        template, PCRPrimer("Illumina RP1", RP1, anneal_nt=16),
        PCRPrimer("Illumina RPI1", RPI1, anneal_nt=21),
        statement=("CEL-Seq2 detailed protocol calls for RP1 and an indexed RPIX; "
                   "Illumina doc 1000000002694 v22 prints the primer bases. "
                   "RPI1 is a representative of the indexed family."))


def final_library() -> Construct:
    """RP1/RPI1-defined ends around the source-declared CEL-Seq2 insert roles."""
    geom = library_pcr_geometry()
    template = geom.template.top()
    assert geom.forward_hit is not None and geom.reverse_hit is not None
    middle = template[geom.forward_hit.end:geom.reverse_hit.start]
    prefix = RP1 + middle[:len(il.TRUSEQ_SMALL_RNA_RA5) - len(RP1.removeprefix(il.P5))]
    if not prefix.endswith(il.TRUSEQ_SMALL_RNA_RA5):
        raise ValueError("RP1 extension fails to recover the CEL-Seq2 RA5 tract")
    if not RPI1.endswith(DESIGN.randomhex_tail[1:]):
        raise ValueError("RPI1 no longer lands on the randomhexRT tail")
    return Construct([
        Segment("Illumina P5", il.P5, "p5"),
        Segment("small-RNA RA5", il.TRUSEQ_SMALL_RNA_RA5, "r1"),
        Segment("UMI", "N" * DESIGN.umi_nt, "umi", placeholder=True,
                feature=C.UMI_FEATURE),
        Segment("cell barcode", "N" * DESIGN.barcode_nt, "cbc",
                placeholder=True, feature=C.CELL_BARCODE),
        Segment("poly(T)/cDNA", "TTTXXXXXXXX...XXXXXXXX", placeholder=True),
        Segment("RPI1 reverse complement", revcomp(RPI1), "r2"),
    ], name="CEL-Seq2 RPI1 library")


FINAL_LIBRARY = final_library()
READ1_T_NT = 3
READ1_CYCLES = DESIGN.umi_nt + DESIGN.barcode_nt + READ1_T_NT
INDEX1_CYCLES = 7
READ2_CYCLES = 36


PROTOCOL_FILE = "Additional file 4 (MOESM4, CEL-Seq2 Protocol)"


def protocol_graph() -> ProtocolGraph:
    """Manual tube/plate CEL-Seq2 workflow from the detailed author protocol.

    Steps, reagents and cleanups follow Additional file 4 (MOESM4 docx text lines
    166-328).  The graph is assembled locally rather than through the shared
    barcoded IVT helper so that the protocol's own heat inactivation, bead cleanups
    and post-IVT ExoSAP-IT treatment stay in their published order.
    """
    protocol, design, doi = PROTOCOL, DESIGN, SOURCE[0]
    source = f"https://doi.org/{doi}"
    graph = ProtocolGraph(protocol, authority=source)
    reviewed: list[str] = []

    def note(*ids: str) -> None:
        reviewed.extend(ids)

    primary = graph.input("mrna", "polyadenylated RNA from one frozen, 65 °C-lysed cell",
                          MaterialKind.BIOLOGICAL, state=C.rt_scene(protocol, design),
                          primary=True)
    rt_primer = graph.input("rt_primer", f"{protocol} RT primer (one barcode per cell, 25 ng/µL)",
                            MaterialKind.OLIGO, state=C.rt_scene(protocol, design))
    ercc = graph.input("ercc_spike_in", "ERCC RNA spike-in (1 µL of 1:1,000,000 per cell)",
                       MaterialKind.BIOLOGICAL)
    rt_enzyme = graph.input("reverse_transcriptase", "SuperScript II (42 °C, 1 h)",
                            MaterialKind.CATALYST)
    first = graph.product("first_strand", "well-barcoded first-strand cDNA:RNA hybrid",
                          MaterialKind.MOLECULAR, state=C.rt_scene(protocol, design))
    free_primers = graph.product("free_rt_primers", "unincorporated RT primers",
                                 MaterialKind.MOLECULAR)
    graph.react("reverse_transcription",
                "Reverse-transcribe with the well-barcoded anchored oligo(dT) primer",
                StepKind.REVERSE_TRANSCRIPTION,
                ((primary, InputRole.TEMPLATE), (ercc, InputRole.SPIKE_IN),
                 (rt_primer, InputRole.STRUCTURAL), (rt_enzyme, InputRole.CATALYST)),
                (first, free_primers))
    note("mrna", "rt_primer", "ercc_spike_in", "reverse_transcriptase",
         "reverse_transcription", "free_rt_primers")

    heat = graph.input("rt_heat_inactivation", "70 °C, 10 min", MaterialKind.CONDITION)
    inactivated = graph.product("inactivated_first_strand",
                                "first-strand cDNA after RT heat inactivation",
                                MaterialKind.MOLECULAR, state=first.state)
    graph.react("inactivate_rt", "Heat-inactivate the reverse transcriptase",
                StepKind.THERMAL_TREATMENT,
                ((first, InputRole.TEMPLATE), (heat, InputRole.CONDITION)),
                (inactivated,), changes_molecule=False)
    note("rt_heat_inactivation", "inactivate_rt")

    second_mix = graph.input(
        "second_strand_enzymes",
        "E. coli DNA polymerase I, E. coli DNA ligase and RNase H (16 °C, 2 h)",
        MaterialKind.CATALYST)
    duplex = graph.product("duplex", "double-stranded barcoded T7 template",
                           MaterialKind.MOLECULAR,
                           state=C.double_stranded_cdna(protocol, design))
    graph.react("second_strand", "Synthesize the second cDNA strand",
                StepKind.EXTENSION,
                ((inactivated, InputRole.TEMPLATE), (second_mix, InputRole.CATALYST)),
                (duplex,))
    note("second_strand_enzymes", "second_strand")

    peers = graph.input("peer_wells", "barcoded duplex cDNA from the other cells of one IVT",
                        MaterialKind.MOLECULAR,
                        state=C.double_stranded_cdna(protocol, design))
    pooled = graph.product("pooled", "pooled well-barcoded duplex cDNA",
                           MaterialKind.MOLECULAR,
                           state=C.double_stranded_cdna(protocol, design))
    graph.react("pool", "Pool cells that go to the same IVT", StepKind.POOL,
                ((duplex, InputRole.TEMPLATE), (peers, InputRole.TEMPLATE)),
                (pooled,), changes_molecule=False)
    note("peer_wells", "pool")

    ampure = graph.input("cdna_ampure",
                         "AMPure XP beads plus 20% PEG8000/2.5 M NaCl bead buffer",
                         MaterialKind.SELECTION)
    clean_duplex = graph.product("clean_duplex", "bead-purified pooled duplex cDNA",
                                 MaterialKind.MOLECULAR,
                                 state=C.double_stranded_cdna(protocol, design))
    cdna_sup = graph.product("cdna_cleanup_supernatant",
                             "bead supernatant and 80% ethanol washes",
                             MaterialKind.WASTE)
    graph.react("cdna_cleanup", "Bead-purify pooled cDNA (elute 6.4 µL)",
                StepKind.CLEANUP,
                ((pooled, InputRole.TEMPLATE), (ampure, InputRole.SELECTION)),
                (clean_duplex, cdna_sup))
    graph.discard(cdna_sup, "unbound material removed with the bead supernatant")
    note("cdna_ampure", "cdna_cleanup", "cdna_cleanup_supernatant")

    t7 = graph.input("t7_polymerase", "MEGAscript T7 enzyme mix and NTPs (37 °C, 13 h)",
                     MaterialKind.CATALYST)
    arna = graph.product("arna", "linearly amplified antisense RNA",
                         MaterialKind.MOLECULAR, state=C.arna(protocol, design))
    graph.react("ivt", "Linearly amplify pooled barcoded cDNA by T7 IVT",
                StepKind.TRANSCRIPTION,
                ((clean_duplex, InputRole.TEMPLATE), (t7, InputRole.CATALYST)),
                (arna,))
    note("t7_polymerase", "ivt")

    exosap = graph.input("exosap_it", "ExoSAP-IT (37 °C, 15 min)", MaterialKind.CATALYST)
    treated = graph.product("exosap_arna", "aRNA after ExoSAP-IT",
                            MaterialKind.MOLECULAR, state=C.arna(protocol, design))
    digested = graph.product("digested_primers",
                             "ExoSAP-IT digestion products of carried-over primers",
                             MaterialKind.MOLECULAR)
    graph.react("exosap", "Treat with ExoSAP-IT \u201cto remove primers\u201d",
                StepKind.DIGESTION,
                ((arna, InputRole.TEMPLATE), (free_primers, InputRole.TEMPLATE),
                 (exosap, InputRole.CATALYST)), (treated, digested))
    note("exosap_it", "exosap", "exosap_arna", "digested_primers")

    shear = graph.input("fragmentation_condition",
                        "200 mM Tris-acetate pH 8.1, 500 mM KOAc, 150 mM MgOAc; 94 °C, 3 min",
                        MaterialKind.CONDITION)
    stop = graph.input("fragmentation_stop", "0.5 M EDTA stop buffer",
                       MaterialKind.CONDITION)
    fragmented = graph.product("fragmented_arna", "fragmented amplified RNA",
                               MaterialKind.MOLECULAR, state=C.arna(protocol, design))
    graph.react("fragment_arna", "Fragment amplified RNA, then stop with EDTA",
                StepKind.FRAGMENTATION,
                ((treated, InputRole.TEMPLATE), (shear, InputRole.CONDITION),
                 (stop, InputRole.CONDITION)), (fragmented,))
    note("fragmentation_condition", "fragmentation_stop", "fragment_arna")

    rnaclean = graph.input("rnaclean_xp", "RNAClean XP beads, 1.8 volumes",
                           MaterialKind.SELECTION)
    clean_fragments = graph.product("clean_fragments", "bead-purified fragmented aRNA",
                                    MaterialKind.MOLECULAR, state=C.arna(protocol, design))
    rna_sup = graph.product("arna_cleanup_supernatant",
                            "bead supernatant, unbound digestion products and 70% ethanol washes",
                            MaterialKind.WASTE)
    graph.react("arna_cleanup", "Bead-purify fragmented aRNA", StepKind.CLEANUP,
                ((fragmented, InputRole.TEMPLATE), (digested, InputRole.TEMPLATE),
                 (rnaclean, InputRole.SELECTION)), (clean_fragments, rna_sup))
    graph.discard(rna_sup, "unbound material removed with the bead supernatant")
    note("rnaclean_xp", "arna_cleanup", "arna_cleanup_supernatant")

    randomhex = graph.input("randomhex_rt", "randomhexRT primer",
                            MaterialKind.OLIGO, state=C.random_rt_scene(design))
    library_rt = graph.input("library_reverse_transcriptase",
                             "SuperScript II (25 °C 10 min, then 42 °C 1 h)",
                             MaterialKind.CATALYST)
    pcr_template = graph.product(
        "pcr_template", "random-primed library-PCR template",
        MaterialKind.MOLECULAR, state=library_pcr_template())
    graph.react("library_reverse_transcription",
                "Reverse-transcribe fragmented aRNA with randomhexRT",
                StepKind.REVERSE_TRANSCRIPTION,
                ((clean_fragments, InputRole.TEMPLATE),
                 (randomhex, InputRole.STRUCTURAL),
                 (library_rt, InputRole.CATALYST)), (pcr_template,))
    note("randomhex_rt", "library_reverse_transcriptase",
         "library_reverse_transcription")
    forward = graph.input("pcr_forward", "Illumina RP1 library-PCR primer",
                          MaterialKind.OLIGO, state=RP1, source=VENDOR_SOURCE)
    reverse = graph.input("pcr_reverse", "Illumina indexed RPI1 library-PCR primer",
                          MaterialKind.OLIGO, state=RPI1, source=VENDOR_SOURCE)
    polymerase = graph.input(
        "pcr_polymerase",
        "Phusion High-Fidelity PCR Master Mix with HF buffer (11 cycles; up to 15)",
        MaterialKind.CATALYST)
    amplified = graph.product("amplified_library", "amplified CEL-Seq2 library",
                              MaterialKind.MOLECULAR, state=FINAL_LIBRARY)
    graph.react("library_pcr", "Amplify and index the sequencing library",
                StepKind.PCR,
                ((pcr_template, InputRole.TEMPLATE),
                 (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (amplified,),
                geometry=library_pcr_geometry(pcr_template.state))
    note("pcr_forward", "pcr_reverse", "pcr_polymerase", "library_pcr",
         "amplified_library")

    beads1 = graph.input("pcr_ampure_1", "AMPure XP beads, 1 volume (first cleanup)",
                         MaterialKind.SELECTION)
    once = graph.product("library_once_cleaned", "library after the first bead cleanup",
                         MaterialKind.MOLECULAR, state=FINAL_LIBRARY)
    sup1 = graph.product("pcr_cleanup_1_supernatant",
                         "first-cleanup bead supernatant and ethanol washes",
                         MaterialKind.WASTE)
    graph.react("pcr_cleanup_1", "Bead-purify the PCR product (elute 25 µL)",
                StepKind.CLEANUP,
                ((amplified, InputRole.TEMPLATE), (beads1, InputRole.SELECTION)),
                (once, sup1))
    graph.discard(sup1, "unbound material removed with the bead supernatant")
    beads2 = graph.input("pcr_ampure_2", "AMPure XP beads, 1 volume (second cleanup)",
                         MaterialKind.SELECTION)
    library = graph.product(
        "library", f"{protocol} sequencing library", MaterialKind.LIBRARY,
        state=FINAL_LIBRARY)
    sup2 = graph.product("pcr_cleanup_2_supernatant",
                         "second-cleanup bead supernatant and ethanol washes",
                         MaterialKind.WASTE)
    graph.react("pcr_cleanup_2", "Repeat the bead cleanup (elute 10 µL)",
                StepKind.CLEANUP,
                ((once, InputRole.TEMPLATE), (beads2, InputRole.SELECTION)),
                (library, sup2))
    graph.discard(sup2, "unbound material removed with the bead supernatant")
    note("pcr_ampure_1", "pcr_cleanup_1", "library_once_cleaned",
         "pcr_cleanup_1_supernatant", "pcr_ampure_2", "pcr_cleanup_2",
         "pcr_cleanup_2_supernatant", "library")
    graph.sequence(library, SequencingPlan(
        limitation=("The published read roles are retained below, but the defining "
                    "sources do not disclose exact sequencing-primer sites.")))
    graph.annotate(SourceDiscrepancy(
        "rt_primer", "Size of the CEL-Seq2 barcode set",
        (SourceClaim(f"{PROTOCOL_FILE}, Primers paragraph (docx text line 48)",
                     "\u201cIf not using all 48 primers, it is recommended to pool at least 10\u201d"),
         SourceClaim(f"{PROTOCOL_FILE}, primer list (docx text lines 51\u2013146)",
                     "96 barcoded primers, 1s\u201396s"),
         SourceClaim("Article Methods, Barcode design (PMC4848782 JATS)",
                     "168 unique barcodes matching the constraints were selected")),
        "The schematic draws one representative variable 6-nt barcode; it does not "
        "choose a set size.",
        "All three statements are in the published article package.",
        "The protocol does not say whether 48 refers to an earlier primer set or a "
        "typical IVT pool; Table S1 (C1) uses barcodes 1\u201348 twice."))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def readout() -> str:
    r1 = [Segment("UMI", "N" * DESIGN.umi_nt, "umi", placeholder=True),
          Segment("cell barcode", "N" * DESIGN.barcode_nt, "cbc", placeholder=True),
          Segment("poly(T)", "T" * READ1_T_NT)]
    r2 = [Segment("transcript", "XXXXXXXX...", placeholder=True)]
    reads = [
        ("Read 1", f"{READ1_CYCLES} nt: {DESIGN.umi_nt}-nt UMI + "
         f"{DESIGN.barcode_nt}-nt cell barcode + three T"),
        ("Index", f"{INDEX1_CYCLES} cycles: indexed library identity"),
        ("Read 2", f"{READ2_CYCLES} nt: transcript sequence from the random-primed end"),
    ]
    return f"""<h2>Read layout</h2>
{sp.unavailable_diagram(FINAL_LIBRARY, "PCR primer positions are source-backed; exact sequencing-primer placement is not yet declared")}
{panel([strand_row(Construct(r1), prefix="Read 1  ", suffix=""),
        strand_row(Construct(r2), prefix="Read 2  ", suffix="")],
       cls="small", caption="Paired-end sequencing keeps cell identity in Read 1 and "
                            "transcript sequence in Read 2.")}
{table(("Read", "Published layout"), reads)}
{sp.read_spec_table([
    sp.read_spec("Read 1", READ1_CYCLES,
                 f"{DESIGN.umi_nt}-nt UMI, {DESIGN.barcode_nt}-nt cell barcode, then three T", "extract"),
    sp.read_spec("Index 1 (i7)", INDEX1_CYCLES, "indexed library identity", "demultiplex"),
    sp.read_spec("Read 2", READ2_CYCLES,
                 "transcript sequence from the random-primed end", "include"),
], "Read roles, cycle counts and contents are source-declared; sequencing-primer placement remains unverified.")}"""


def render_page() -> str:
    from page_parts import oligos as family_oligos
    from protocol_graph import AnnotationSlot
    from protocol_views import render_annotations

    graph = protocol_graph()
    source = (f'Defining source: <a href="https://doi.org/{SOURCE[0]}">'
              f'Hashimshony et al., <i>{SOURCE[1]}</i> ({SOURCE[2]})</a>; '
              f'RP1/RPI1 bases: <a href="{VENDOR_SOURCE}">Illumina TruSeq Small RNA '
              'adapter sequences</a>, doc 1000000002694 v22.')
    caveat = ("<div class='caveat'><b>Source and variant boundary.</b> The CEL-Seq2 author "
              "protocol calls for Illumina RP1 and an indexed RPIX but refers their "
              "sequences to Illumina. The vendor table supplies the RP1/RPI family; "
              "RPI1 is a representative member, not an identified index used in the "
              "2016 experiments. This schematic follows the manual well workflow; "
              "the paper also supplies a distinct Fluidigm C1 protocol in Additional file 6. "
              "The transcript insert and cell UMI/barcode remain variable; "
              "sequencing-primer placement is not yet declared.</div>")
    pcr_oligos = ("<h2>Library PCR primers</h2><seq>" +
                  oligo("Illumina RP1", [Segment("RP1", RP1)]) +
                  oligo("Illumina RPI1", [Segment("RPI1", RPI1)]) + "</seq>")
    return "\n".join((head(f"{PROTOCOL} library chemistry"),
                      '<div class="wrap" data-protocol-complete="true">',
                      f"<h1>{PROTOCOL} &mdash; barcoded 3' RNA-seq by T7 linear amplification</h1>",
                      info(source), caveat,
                      render_annotations(graph, AnnotationSlot.INTRO), family_oligos(PROTOCOL, DESIGN), pcr_oligos,
                      "<h2>Library construction</h2>", protocol_graph_panels(graph),
                      readout(), "</div>"))
