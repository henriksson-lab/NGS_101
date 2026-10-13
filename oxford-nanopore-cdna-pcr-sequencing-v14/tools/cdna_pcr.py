"""Oxford Nanopore SQK-PCS114 cDNA-PCR workflow."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Construct, Row, Segment, feature
from nanopore_library import motor_sequencing_plan
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from protocol_model import molecule, role
from rna_special import nanopore_entry_rows, template_switch_scene

TITLE = "Oxford Nanopore cDNA-PCR Sequencing V14 — SQK-PCS114"
NOTES = "01_cdna-pcr.html"
SOURCE = 'Commercial protocol: Oxford Nanopore <a href="https://nanoporetech.com/document/pcr-cdna-sequencing-v14-sqk-pcs114">cDNA-PCR Sequencing V14, SQK-PCS114</a>.'
SUMMARY = "Strand-switch reverse transcription adds a UMI and terminal handles, PCR selects and amplifies full-length cDNA while adding rapid-attachment tags, and a motor-loaded Rapid Adapter prepares one cDNA strand for nanopore sequencing."
CAVEAT = "SQK-PCS114 oligo and adapter bases are proprietary. Roles, UMI placement and reaction order are vendor-supported; unknown bases remain role placeholders."
FINAL_LIBRARY = Construct([
    Segment("motor-loaded Rapid Adapter", "[motor-loaded adapter]", "r2",
            placeholder=True, inferred=True),
    Segment("5-prime rapid-attachment tag", "[rapid-attachment tag]", "r1",
            placeholder=True, inferred=True),
    Segment("strand-switch UMI", "[UMI]", "umi", placeholder=True, inferred=True,
            feature=feature("umi", "umi", "random", note="exact length is proprietary")),
    Segment("full-length cDNA", "X" * 36, placeholder=True),
    Segment("3-prime rapid-attachment tag", "[rapid-attachment tag]", "r1",
            placeholder=True, inferred=True),
], name="SQK-PCS114 adapted cDNA")
AMPLIFIED_CDNA = Construct(list(FINAL_LIBRARY)[1:], name="SQK-PCS114 amplified cDNA")
FINAL_CAPTION = ("INFERRED — the vendor-supported role order after PCR and Rapid Adapter "
                 "attachment. Proprietary sequences and UMI length are not asserted.")
SEQ_PRIMERS = ()
SEQUENCING_ENDING = "There is no synthesis sequencing primer. The motor-loaded Rapid Adapter captures one strand of amplified cDNA at the pore and meters that strand through the nanopore; the adapter-bearing end enters first."

def sections():
    return [
        ("Reverse-transcribe and strand-switch", template_switch_scene(inferred=True, umi=True).rows(), "INFERRED — proprietary CRTA/RTP/SSPII bases are not drawn. The current vendor protocol states that strand switching incorporates a UMI."),
        ("Select full-length molecules by PCR", [Row(chunks=[("5′ rapid-attachment tag — UMI — full-length cDNA — rapid-attachment tag 3′", None, True)]), Row(chunks=[("                             PCR amplification", "umi", False)])], "INFERRED — cDNA Primer amplification enriches molecules carrying both terminal handles and installs the rapid-adapter attachment tags."),
        ("Attach the motor-loaded Rapid Adapter", nanopore_entry_rows("amplified cDNA", amplified=True), "The Rapid Adapter assembles at the PCR-added terminal tag; ** marks the adapter attachment boundary."),
    ]


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE)
    rna = g.input("rna", "poly(A) RNA", MaterialKind.MOLECULAR,
                  state=molecule("poly(A) RNA", (role("full-length transcript"),
                                                  role("poly(A) tail"))), primary=True)
    rt_primer = g.input("rt_primer", "proprietary poly(dT) RT primer",
                        MaterialKind.OLIGO, inferred=True,
                        source=SOURCE + " (CRTA/RTP bases proprietary)")
    rt = g.input("rt", "strand-switching reverse transcriptase", MaterialKind.CATALYST)
    first = g.product("first_strand", "first-strand cDNA at the RNA 5-prime end",
                      MaterialKind.MOLECULAR)
    g.react("reverse_transcribe", "Reverse-transcribe poly(A) RNA",
            StepKind.REVERSE_TRANSCRIPTION,
            ((rna, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
             (rt, InputRole.CATALYST)), (first,))
    tso = g.input("tso", "UMI-bearing strand-switch oligo", MaterialKind.OLIGO,
                  inferred=True, source=SOURCE + " (SSPII bases and UMI length proprietary)")
    full = g.product("full_cdna", "full-length UMI-bearing cDNA",
                     MaterialKind.MOLECULAR)
    g.react("template_switch", "Template-switch at the completed cDNA end",
            StepKind.TEMPLATE_SWITCHING,
            ((first, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL),
             (rt, InputRole.CATALYST)), (full,))
    forward = g.input("forward_primer", "cDNA PCR forward primer", MaterialKind.OLIGO,
                      inferred=True, source=SOURCE + " (primer bases proprietary)")
    reverse = g.input("reverse_primer", "cDNA PCR reverse primer", MaterialKind.OLIGO,
                      inferred=True, source=SOURCE + " (primer bases proprietary)")
    polymerase = g.input("polymerase", "cDNA PCR polymerase", MaterialKind.CATALYST)
    amplified = g.product("amplified", "PCR-amplified full-length cDNA",
                          MaterialKind.MOLECULAR, state=AMPLIFIED_CDNA, inferred=True,
                          source=SOURCE + " (terminal handle bases proprietary)")
    g.react("pcr", "Amplify full-length cDNA and install rapid-attachment tags",
            StepKind.PCR,
            ((full, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
             (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
            (amplified,), geometry=PCRGeometry.unavailable(
                AMPLIFIED_CDNA, reason="The vendor withholds both cDNA PCR-primer sequences."))
    rapid = g.input("rapid_adapter", "motor-loaded Rapid Adapter",
                    MaterialKind.ADAPTER, inferred=True,
                    source=SOURCE + " (adapter bases and attachment chemistry proprietary)")
    library = g.product("library", "motor-adapted amplified-cDNA library",
                        MaterialKind.LIBRARY, state=FINAL_LIBRARY, inferred=True,
                        source=SOURCE + " (adapter bases proprietary)")
    g.react("attach", "Attach the motor-loaded Rapid Adapter",
            StepKind.ADAPTER_ATTACHMENT,
            ((amplified, InputRole.TEMPLATE), (rapid, InputRole.STRUCTURAL)), (library,))
    limitation = ("Nanopore sequencing is motor-driven rather than primer-extended; "
                  "the current Rapid Adapter bases are proprietary.")
    g.sequence(library, motor_sequencing_plan(FINAL_LIBRARY, limitation=limitation))
    g.declare_reviewed_scope(("rna", "rt_primer", "rt", "reverse_transcribe", "tso",
                              "template_switch", "forward_primer", "reverse_primer",
                              "polymerase", "pcr", "rapid_adapter", "attach", "library"),
                             SOURCE)
    return g.require_publishable()
