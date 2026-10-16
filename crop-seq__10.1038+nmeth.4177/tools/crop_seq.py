"""Molecular model for the original Drop-seq implementation of CROP-seq."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import nextera as nx
import seqprimers as sp
from crispr import SCAFFOLD_V1
from chemdraw import Construct, Segment, feature, revcomp
from lentivirus import TransferGenome
from pcr import PCRGeometry
from protocol_graph import (Disposition, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind, Terminal)

TITLE = "CROP-seq — pooled CRISPR screening with single-cell RNA readout"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/nmeth.4177">Datlinger et al., <i>Nature Methods</i> (2017)</a>.'
DROPSEQ_BEAD_SOURCE = ('<a href="https://mccarrolllab.org/wp-content/uploads/2015/05/'
                       'Online-Dropseq-Protocol-v.-3.1-Dec-2015.pdf">Macosko and Goldman, '
                       'Drop-seq lab protocol v3.1 (2015), primer table</a>')

GUIDE_FEATURE = feature("guide", "guide_barcode", "unknown")
CELL_BARCODE = feature("cell_barcode", "cell_barcode", "random")
UMI_FEATURE = feature("umi", "umi", "random")
I7_FEATURE = feature("sample_index_i7", "sample_index", "unknown")
FEATURES = {"guide": GUIDE_FEATURE, "guide cDNA": GUIDE_FEATURE,
            "cell barcode": CELL_BARCODE, "UMI": UMI_FEATURE,
            "i7 reverse complement": I7_FEATURE}


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    kw.setdefault("feature", FEATURES.get(name))
    return Segment(name=name, top=top, tag=tag, **kw)


def role(name: str, text: str, tag: str | None = None) -> Segment:
    return seg(name, f"[{text}]", tag, placeholder=True)


GIBSON_5 = "TGGAAAGGACGAAACACCG"
GIBSON_3 = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGC"
GUIDE_LEN = 20
SMART_HANDLE = "AAGCAGTGGTATCAACGCAGAGT"
P5_SPACER = "GCCTGTCCGCGG"
CUSTOM_R1 = P5_SPACER + SMART_HANDLE + "AC"


def pooled_guide_oligo() -> list[Segment]:
    return [seg("hU6 homology + G", GIBSON_5, "r3"),
            seg("guide", "N" * GUIDE_LEN, "cbc", placeholder=True),
            seg("scaffold homology", GIBSON_3, "r2")]


def guide_expression_cassette() -> tuple[Segment, ...]:
    return (seg("hU6 promoter tail + G", GIBSON_5, "r3"),
            seg("guide", "N" * GUIDE_LEN, "cbc", placeholder=True),
            seg("gRNA scaffold", SCAFFOLD_V1, "r2"))


def transfer_genome() -> TransferGenome:
    return TransferGenome(
        u3=(*guide_expression_cassette(), role("U3 remainder", "U3 remainder", "me")),
        r=(role("R", "R", "me"),), u5=(role("U5", "U5", "me"),),
        internal=(role("vector body", "psi ... EF1a-Puro ... WPRE", "w1"),),
        name="CROP-seq")


def pol3_guide_transcript() -> Construct:
    return Construct([seg("guide", "N" * GUIDE_LEN, "cbc", placeholder=True),
                      seg("gRNA scaffold", SCAFFOLD_V1, "r2"),
                      role("Pol III termination", "U6 termination")],
                     name="functional Pol III guide RNA")


def pol2_reporter_transcript() -> Construct:
    return Construct([role("puromycin transcript", "EF1a--Puro--WPRE", "w1"),
                      *guide_expression_cassette(), role("3' LTR", "3' LTR", "me"),
                      seg("poly(A)", "A" * 24)],
                     name="polyadenylated Pol II reporter RNA")


def guide_bearing_cdna_fragment() -> Construct:
    """Representative bead-end fragment cut at the hU6/guide boundary."""
    return Construct([
        seg("SMART handle", SMART_HANDLE, "tso"), seg("batch-B constant", "AC"),
        seg("cell barcode", "B" * 12, "cbc", placeholder=True),
        seg("UMI", "U" * 8, "umi", placeholder=True), seg("poly(T)", "T" * 30),
        seg("3' transcript cDNA", "X" * 16, "me", placeholder=True),
        seg("scaffold cDNA", revcomp(SCAFFOLD_V1), "r2"),
        seg("guide cDNA", "N" * GUIDE_LEN, "cbc", placeholder=True),
        seg("mosaic end reverse complement", nx.ME_RC, "me"),
        seg("s7 reverse complement", nx.S7_RC, "s7")],
        name="guide-bearing bead-end fragment")


def sequencing_primers():
    return (sp.custom("Read 1", "Drop-seq Custom Read1", CUSTOM_R1,
                      "CROP-seq Supplementary Protocol oligo-order table"),
            sp.NEXTERA["I1"], sp.NEXTERA["R2"])


def final_library() -> Construct:
    lib = Construct([seg("P5", il.P5, "p5"), seg("Read-1 spacer", P5_SPACER, "r1"),
                     *list(guide_bearing_cdna_fragment()),
                     seg("i7 reverse complement", "I" * 8, "cbc", placeholder=True),
                     seg("P7 reverse complement", il.P7_RC, "p7")],
                    name="CROP-seq guide-bearing Drop-seq library")
    problems = sp.verify(lib, sequencing_primers(),
                         required_roles=("Read 1", "Index 1 (i7)", "Read 2"))
    if problems:
        raise ValueError("invalid CROP-seq final library: " + "; ".join(problems))
    return lib


def read1_layout() -> list[tuple[str, str]]:
    return [(x.cycles, x.feature.label) for x in
            sp.feature_spans(final_library(), sequencing_primers(), {"Read 1": 20})]


def protocol_graph():
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    provirus = graph.input("provirus", "cell carrying an integrated CROP-seq provirus",
                           MaterialKind.BIOLOGICAL, state=transfer_genome().provirus(),
                           primary=True)
    pol3 = graph.input("pol3", "RNA polymerase III", MaterialKind.CATALYST)
    guide = graph.product("functional_guide", "functional non-polyadenylated guide RNA",
                          MaterialKind.MOLECULAR, state=pol3_guide_transcript())
    graph.react("guide_transcription", "Transcribe the functional guide from hU6",
                StepKind.TRANSCRIPTION,
                ((provirus, InputRole.TEMPLATE), (pol3, InputRole.CATALYST)), (guide,))
    graph.mark_terminal(Terminal(
        guide.id, Disposition.STORED,
        "the Pol III guide executes the perturbation but is not the sequenced reporter branch"))

    pol2 = graph.input("pol2", "RNA polymerase II", MaterialKind.CATALYST)
    reporter = graph.product("reporter_rna", "polyadenylated guide-identifying reporter RNA",
                             MaterialKind.MOLECULAR, state=pol2_reporter_transcript())
    graph.react("reporter_transcription", "Transcribe the polyadenylated reporter",
                StepKind.TRANSCRIPTION,
                ((provirus, InputRole.TEMPLATE), (pol2, InputRole.CATALYST)), (reporter,))
    bead_oligo = graph.input("bead_oligo", "Drop-seq bead oligo-dT with cell barcode and UMI",
                             MaterialKind.OLIGO, source=DROPSEQ_BEAD_SOURCE)
    rt = graph.input("reverse_transcriptase", "template-switching reverse transcriptase",
                     MaterialKind.CATALYST)
    cdna = graph.product("barcoded_cdna", "cell-barcoded reporter cDNA",
                         MaterialKind.MOLECULAR)
    graph.react("reverse_transcribe", "Capture and reverse-transcribe the reporter on a bead",
                StepKind.REVERSE_TRANSCRIPTION,
                ((reporter, InputRole.TEMPLATE), (bead_oligo, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (cdna,))
    tso = graph.input("tso", "Drop-seq template-switch oligo", MaterialKind.OLIGO)
    switched = graph.product("switched_cdna", "template-switched guide-reporter cDNA",
                             MaterialKind.MOLECULAR)
    graph.react("template_switch", "Install the distal SMART handle",
                StepKind.TEMPLATE_SWITCHING,
                ((cdna, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (switched,))
    smart_f = graph.input("smart_forward", "SMART cDNA PCR forward primer", MaterialKind.OLIGO)
    smart_r = graph.input("smart_reverse", "SMART cDNA PCR reverse primer", MaterialKind.OLIGO)
    preamp_pol = graph.input("preamp_polymerase", "cDNA amplification polymerase",
                             MaterialKind.CATALYST)
    amplified = graph.product("amplified_cdna", "amplified guide-reporter cDNA",
                              MaterialKind.MOLECULAR)
    graph.react("preamp", "Amplify Drop-seq cDNA", StepKind.PCR,
                ((switched, InputRole.TEMPLATE), (smart_f, InputRole.STRUCTURAL),
                 (smart_r, InputRole.STRUCTURAL), (preamp_pol, InputRole.CATALYST)),
                (amplified,), geometry=PCRGeometry.unavailable(
                    final_library(), reason="The intermediate SMART amplicon is not reconstructed as a full duplex."))
    transposome = graph.input("nextera_transposome", "Nextera XT transposome",
                              MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagged_cdna", "Nextera-tagged reporter cDNA",
                           MaterialKind.MOLECULAR, state=guide_bearing_cdna_fragment())
    tag_waste = graph.product("tagmentation_waste", "non-amplifiable Tn5 products",
                              MaterialKind.WASTE)
    graph.react("tagment", "Tagment the amplified cDNA", StepKind.TAGMENTATION,
                ((amplified, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, tag_waste))
    graph.discard(tag_waste, "homotypic and otherwise non-amplifiable Tn5 products")
    selective = graph.input("selective_p5", "P5-SMART selective enrichment primer",
                            MaterialKind.OLIGO)
    n70x = graph.input("n70x", "indexed N70x enrichment primer", MaterialKind.OLIGO)
    library_pol = graph.input("library_polymerase", "library PCR polymerase",
                              MaterialKind.CATALYST)
    final = final_library()
    library = graph.product("library", final.name, MaterialKind.LIBRARY, state=final)
    graph.react("selective_pcr", "Select guide-bearing bead-end fragments",
                StepKind.PCR,
                ((tagged, InputRole.TEMPLATE), (selective, InputRole.STRUCTURAL),
                 (n70x, InputRole.STRUCTURAL), (library_pol, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    final, reason="The selective primer identities are known, but the tagmented precursor is heterogeneous."))
    graph.sequence(library, SequencingPlan(
        library=final, primers=tuple(sequencing_primers()),
        read_lengths={"Read 1": 20, "Index 1 (i7)": 8, "Read 2": 64}))
    graph.declare_reviewed_scope((
        "provirus", "pol3", "guide_transcription", "pol2", "reporter_transcription",
        "bead_oligo", "reverse_transcriptase", "reverse_transcribe", "tso",
        "template_switch", "smart_forward", "smart_reverse", "preamp_polymerase",
        "preamp", "nextera_transposome", "tn5", "tagment", "selective_p5", "n70x",
        "library_polymerase", "selective_pcr", "library"), SOURCE)
    return graph.require_publishable()
