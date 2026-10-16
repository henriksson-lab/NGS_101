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
                            StepKind)
from rna_amplification import barcoded_ivt_front

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


def protocol_graph() -> ProtocolGraph:
    protocol, design, doi = PROTOCOL, DESIGN, SOURCE[0]
    source = f"https://doi.org/{doi}"
    graph = ProtocolGraph(protocol, authority=source)

    primary = graph.input("mrna", "polyadenylated RNA from one well",
                          MaterialKind.BIOLOGICAL, state=C.rt_scene(protocol, design),
                          primary=True)
    arna, reviewed = barcoded_ivt_front(
        graph, primary, rt_primer_label=f"{protocol} RT primer",
        rt_state=C.rt_scene(protocol, design), first_state=C.rt_scene(protocol, design),
        duplex_state=C.double_stranded_cdna(protocol, design),
        arna_state=C.arna(protocol, design))
    reviewed.insert(0, "mrna")

    shear = graph.input("fragmentation_condition", "RNA fragmentation condition",
                        MaterialKind.CONDITION)
    fragmented = graph.product("fragmented_arna", "fragmented amplified RNA",
                               MaterialKind.MOLECULAR, state=C.arna(protocol, design))
    graph.react("fragment_arna", "Fragment amplified RNA", StepKind.FRAGMENTATION,
                ((arna, InputRole.TEMPLATE), (shear, InputRole.CONDITION)),
                (fragmented,))
    reviewed.extend(("fragmentation_condition", "fragment_arna"))

    final_state = FINAL_LIBRARY
    randomhex = graph.input("randomhex_rt", "randomhexRT primer",
                            MaterialKind.OLIGO,
                            state=C.random_rt_scene(design))
    library_rt = graph.input("library_reverse_transcriptase",
                             "library reverse transcriptase",
                             MaterialKind.CATALYST)
    pcr_template = graph.product(
        "pcr_template", "random-primed library-PCR template",
        MaterialKind.MOLECULAR, state=library_pcr_template())
    graph.react("library_reverse_transcription",
                "Reverse-transcribe fragmented aRNA with randomhexRT",
                StepKind.REVERSE_TRANSCRIPTION,
                ((fragmented, InputRole.TEMPLATE),
                 (randomhex, InputRole.STRUCTURAL),
                 (library_rt, InputRole.CATALYST)), (pcr_template,))
    reviewed.extend(("randomhex_rt", "library_reverse_transcriptase",
                     "library_reverse_transcription"))
    forward = graph.input("pcr_forward", "Illumina RP1 library-PCR primer",
                          MaterialKind.OLIGO, state=RP1, source=VENDOR_SOURCE)
    reverse = graph.input("pcr_reverse", "Illumina indexed RPI1 library-PCR primer",
                          MaterialKind.OLIGO, state=RPI1, source=VENDOR_SOURCE)
    polymerase = graph.input("pcr_polymerase", "PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product(
        "library", f"{protocol} sequencing library", MaterialKind.LIBRARY,
        state=final_state)
    graph.react("library_pcr", "Amplify and index the sequencing library",
                StepKind.PCR,
                ((pcr_template, InputRole.TEMPLATE),
                 (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (library,),
                geometry=library_pcr_geometry(pcr_template.state))
    graph.sequence(library, SequencingPlan(
        limitation=("The published read roles are retained below, but the defining "
                    "sources do not disclose exact sequencing-primer sites.")))
    reviewed.extend(("pcr_forward", "pcr_reverse", "pcr_polymerase",
                     "library_pcr", "library"))
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
                      info(source), caveat, family_oligos(PROTOCOL, DESIGN), pcr_oligos,
                      "<h2>Library construction</h2>", protocol_graph_panels(graph),
                      readout(), "</div>"))
