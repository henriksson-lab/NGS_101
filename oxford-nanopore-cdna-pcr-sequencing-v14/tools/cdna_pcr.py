"""Oxford Nanopore SQK-PCS114 cDNA-PCR workflow."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Construct, Row, Segment, feature, revcomp
from nanopore_library import (CDNA_PRIMER_COMMON, CDNA_PRIMER_FORWARD,
                              CDNA_PRIMER_REVERSE, CRTA, RAPID_ADAPTER, RTP,
                              SSPII_AS_WRITTEN, SSPII_UMI,
                              motor_sequencing_plan)
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from protocol_model import molecule, role
from rna_special import nanopore_entry_rows, template_switch_scene

TITLE = "Oxford Nanopore cDNA-PCR Sequencing V14 — SQK-PCS114"
NOTES = "01_cdna-pcr.html"
SOURCE = 'Commercial protocol: Oxford Nanopore <a href="https://nanoporetech.com/document/pcr-cdna-sequencing-v14-sqk-pcs114">cDNA-PCR Sequencing V14, SQK-PCS114</a>.'
SUMMARY = "Strand-switch reverse transcription adds a UMI and terminal handles, PCR selects and amplifies full-length cDNA while adding rapid-attachment tags, and a motor-loaded Rapid Adapter prepares one cDNA strand for nanopore sequencing."
CAVEAT = ("ONT publishes CRTA, RTP, SSPII, both cDNA primers and the Rapid Adapter "
          "sequence. The complete motor-bearing Rapid Adapter assembly and its "
          "attachment chemistry have not been verified from the inspected vendor documents.")
PCR_TEMPLATE = Construct([
    Segment("cDNA primer forward", CDNA_PRIMER_FORWARD, "r1"),
    Segment("full-length cDNA", "X" * 36, placeholder=True),
    Segment("strand-switch UMI", SSPII_UMI, "umi", placeholder=True,
            feature=feature("umi", "umi", "random",
                            note="published partially degenerate V/T pattern")),
    Segment("cDNA primer reverse site", revcomp(CDNA_PRIMER_REVERSE), "r2"),
], name="SQK-PCS114 PCR product")
PCR_GEOMETRY = PCRGeometry.exact(
    PCR_TEMPLATE,
    PCRPrimer("cDNA Primer forward", CDNA_PRIMER_FORWARD, anneal_nt=27),
    PCRPrimer("cDNA Primer reverse", CDNA_PRIMER_REVERSE, anneal_nt=23),
    statement="Both cDNA Primer sequences are published in ONT CHTD_500 section 15.")
FINAL_LIBRARY = Construct([
    Segment("motor-loaded Rapid Adapter", "[motor-loaded adapter]", "r2",
            placeholder=True, inferred=True),
    *PCR_TEMPLATE.segments,
], name="SQK-PCS114 adapted cDNA")
AMPLIFIED_CDNA = Construct(list(FINAL_LIBRARY)[1:], name="SQK-PCS114 amplified cDNA")
FINAL_CAPTION = ("The PCR product has the two published cDNA Primer-defined ends and "
                 "published SSPII UMI pattern. The bracketed motor assembly remains unverified pending a component-level source.")
SEQ_PRIMERS = ()
SEQUENCING_MECHANISM = "nanopore_motor"
MOTOR_ADAPTER_NAME = "Rapid Adapter"
SEQUENCING_ENDING = "There is no synthesis sequencing primer. The motor-loaded Rapid Adapter captures one strand of amplified cDNA at the pore and meters that strand through the nanopore; the adapter-bearing end enters first."

def sections():
    return [
        ("Reverse-transcribe and strand-switch", template_switch_scene(inferred=False, umi=True).rows(), "CRTA is ligated at poly(A), RTP primes from the opened CRTA, and the published SSPII contributes its fixed handle, V/T UMI pattern and terminal modified GGG."),
        ("Select full-length molecules by PCR", PCR_GEOMETRY.rows(), "The two published cDNA primers face inward and define the amplified molecule's ends."),
        ("Attach the motor-loaded Rapid Adapter", nanopore_entry_rows("amplified cDNA", amplified=True), "The Rapid Adapter assembles at the PCR-added terminal tag; ** marks the adapter attachment boundary."),
    ]


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE)
    rna = g.input("rna", "poly(A) RNA", MaterialKind.MOLECULAR,
                  state=molecule("poly(A) RNA", (role("full-length transcript"),
                                                  role("poly(A) tail"))), primary=True)
    crta = g.input("crta", f"cDNA RT Adapter (CRTA): {CRTA}",
                   MaterialKind.ADAPTER, source=SOURCE)
    rt_primer = g.input("rt_primer", f"RT Primer (RTP): {RTP}",
                        MaterialKind.OLIGO, source=SOURCE)
    rt = g.input("rt", "strand-switching reverse transcriptase", MaterialKind.CATALYST)
    first = g.product("first_strand", "first-strand cDNA at the RNA 5-prime end",
                      MaterialKind.MOLECULAR)
    g.react("reverse_transcribe", "Reverse-transcribe poly(A) RNA",
            StepKind.REVERSE_TRANSCRIPTION,
            ((rna, InputRole.TEMPLATE), (crta, InputRole.STRUCTURAL),
             (rt_primer, InputRole.STRUCTURAL),
             (rt, InputRole.CATALYST)), (first,))
    tso = g.input("tso", f"Strand Switching Primer II: {SSPII_AS_WRITTEN}",
                  MaterialKind.OLIGO, source=SOURCE)
    full = g.product("full_cdna", "full-length UMI-bearing cDNA",
                     MaterialKind.MOLECULAR)
    g.react("template_switch", "Template-switch at the completed cDNA end",
            StepKind.TEMPLATE_SWITCHING,
            ((first, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL),
             (rt, InputRole.CATALYST)), (full,))
    forward = g.input("forward_primer", f"cDNA Primer forward: {CDNA_PRIMER_FORWARD}",
                      MaterialKind.OLIGO, source=SOURCE)
    reverse = g.input("reverse_primer", f"cDNA Primer reverse: {CDNA_PRIMER_REVERSE}",
                      MaterialKind.OLIGO, source=SOURCE)
    polymerase = g.input("polymerase", "cDNA PCR polymerase", MaterialKind.CATALYST)
    amplified = g.product("amplified", "PCR-amplified full-length cDNA",
                          MaterialKind.MOLECULAR, state=AMPLIFIED_CDNA, source=SOURCE)
    g.react("pcr", "Amplify full-length cDNA and install rapid-attachment tags",
            StepKind.PCR,
            ((full, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
             (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
            (amplified,), geometry=PCR_GEOMETRY)
    rapid = g.input("rapid_adapter", "motor-loaded Rapid Adapter",
                    MaterialKind.ADAPTER, inferred=True,
                    source=SOURCE + f" (published RA strand {RAPID_ADAPTER}; motor assembly unverified in inspected sources)")
    library = g.product("library", "motor-adapted amplified-cDNA library",
                        MaterialKind.LIBRARY, state=FINAL_LIBRARY, inferred=True,
                        source=SOURCE + " (complete motor-bearing assembly unverified in inspected sources)")
    g.react("attach", "Attach the motor-loaded Rapid Adapter",
            StepKind.ADAPTER_ATTACHMENT,
            ((amplified, InputRole.TEMPLATE), (rapid, InputRole.STRUCTURAL)), (library,))
    limitation = ("Nanopore sequencing is motor-driven rather than primer-extended; "
                  "the RA oligo is published, but the complete motor-bearing assembly is not.")
    g.sequence(library, motor_sequencing_plan(FINAL_LIBRARY, limitation=limitation))
    g.declare_reviewed_scope(("rna", "crta", "rt_primer", "rt", "reverse_transcribe", "tso",
                              "template_switch", "forward_primer", "reverse_primer",
                              "polymerase", "pcr", "rapid_adapter", "attach", "library"),
                             SOURCE)
    return g.require_publishable()


def oligos():
    return [
        f"CRTA: 5′-{CRTA}-3′",
        f"RTP: 5′-{RTP}-3′",
        f"SSPII (as published): 5′-{SSPII_AS_WRITTEN}-3′",
        f"cDNA Primer forward: 5′-{CDNA_PRIMER_FORWARD}-3′",
        f"cDNA Primer reverse: 5′-{CDNA_PRIMER_REVERSE}-3′",
        f"Rapid Adapter / RAP T strand: 5′-{RAPID_ADAPTER}-3′",
    ]
