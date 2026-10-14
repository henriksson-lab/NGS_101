"""Oxford Nanopore SQK-DCS109 direct-cDNA workflow."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from adapter_ligation import reported_adapter_end
from chemdraw import Construct, Row, Segment
from nanopore_library import motor_library, motor_sequencing_plan, role
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from protocol_model import molecule, role as molecular_role
from rna_special import nanopore_entry_rows, template_switch_scene

TITLE = "Oxford Nanopore Direct cDNA Sequencing — SQK-DCS109"
NOTES = "01_direct-cdna.html"
SOURCE = 'Commercial protocol: Oxford Nanopore <a href="https://nanoporetech.com/document/chemistry-technical-document">Chemistry Technical Document — Direct cDNA Sequencing Kit SQK-DCS109</a>.'
SUMMARY = "Poly(A) RNA is reverse-transcribed with strand switching, RNA is removed and second-strand DNA is synthesized. The PCR-free duplex cDNA is end-repaired, dA-tailed and ligated to a motor-loaded dT-overhang sequencing adapter."
CAVEAT = "SQK-DCS109 is a historical commercial protocol and its adapter bases are proprietary. The page distinguishes it from cDNA-PCR: no PCR or rapid-attachment tag is introduced."
FINAL_LIBRARY = motor_library([
    role("dT sequencing-adapter end", "[dT sequencing-adapter end]", "r2"),
    Segment("dT:dA ligation junction", "T", bottom="A"),
    role("PCR-free duplex cDNA", "[full-length duplex cDNA]"),
], adapter="sequencing adapter", name="SQK-DCS109 pore-entry strand")
SEQ_PRIMERS = ()
SEQUENCING_MECHANISM = "nanopore_motor"
SEQUENCING_ENDING = "There is no synthesis sequencing primer. A motor in the ligated sequencing adapter presents one strand of the PCR-free duplex cDNA to the pore. The adapter-ligated end enters first."
SEQUENCING_UNAVAILABLE = ("Nanopore sequencing is motor-driven; SQK-DCS109 adapter bases "
                          "and the motor-bearing strand are proprietary.")
FINAL_CAPTION = "The motor-bearing dT adapter is ligated to dA-tailed PCR-free duplex cDNA."
ADAPTER_END = reported_adapter_end(
    "SQK-DCS109 sequencing adapter", kind="cohesive", polarity="3-prime",
    terminal_sequence="T",
    statement="The historical kit uses a dT-overhang sequencing adapter on dA-tailed cDNA.")

def sections():
    return [
        ("Make first-strand cDNA by strand switching", template_switch_scene(inferred=True).rows(), "INFERRED — primer bases are proprietary; their published functions are poly(A)-primed RT and capture of a completed cDNA end by strand switching."),
        ("Remove RNA and synthesize the second strand", [Row(chunks=[("RNA:cDNA hybrid → RNase cocktail → single-stranded cDNA", None, False)]), Row(chunks=[("second-strand synthesis → full-length duplex cDNA", "r1", False)])], "Unlike direct RNA sequencing, the RNA is removed and DNA is the pore substrate."),
        ("End-repair and dA-tail the duplex", [Row(chunks=[("5′ — full-length duplex cDNA — A 3′", None, False)]), Row(chunks=[("3′ — full-length duplex cDNA — T 5′", None, False)])], "Repair creates ligatable ends and dA tailing makes them compatible with the dT-overhang sequencing adapter."),
        ("Ligate the motor-loaded sequencing adapter", nanopore_entry_rows("PCR-free cDNA", amplified=False), "A dT-tailed sequencing adapter ligates to a dA-tailed cDNA end; ** marks the ligation boundary."),
    ]


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE)
    rna = g.input("rna", "poly(A) RNA", MaterialKind.MOLECULAR,
                  state=molecule("poly(A) RNA", (molecular_role("transcript"),
                                                  molecular_role("poly(A) tail"))),
                  primary=True)
    primer = g.input("rt_primer", "proprietary poly(dT) RT primer",
                     MaterialKind.OLIGO, inferred=True,
                     source=SOURCE + " (primer bases proprietary)")
    rt = g.input("rt", "strand-switching reverse transcriptase", MaterialKind.CATALYST)
    hybrid = g.product("hybrid", "full-length RNA:cDNA hybrid", MaterialKind.MOLECULAR)
    g.react("reverse_transcribe", "Make first-strand cDNA by strand switching",
            StepKind.REVERSE_TRANSCRIPTION,
            ((rna, InputRole.TEMPLATE), (primer, InputRole.STRUCTURAL),
             (rt, InputRole.CATALYST)), (hybrid,))
    rnase = g.input("rnase", "RNA-removal enzyme cocktail", MaterialKind.CATALYST)
    single = g.product("single_cdna", "single-stranded full-length cDNA",
                       MaterialKind.MOLECULAR)
    g.react("remove_rna", "Digest the RNA template", StepKind.DIGESTION,
            ((hybrid, InputRole.TEMPLATE), (rnase, InputRole.CATALYST)), (single,))
    second_pol = g.input("second_strand_polymerase", "second-strand DNA polymerase",
                         MaterialKind.CATALYST)
    duplex = g.product("duplex", "full-length duplex cDNA", MaterialKind.MOLECULAR)
    g.react("second_strand", "Synthesize the second cDNA strand", StepKind.EXTENSION,
            ((single, InputRole.TEMPLATE), (second_pol, InputRole.CATALYST)), (duplex,))
    end_mix = g.input("end_prep_mix", "end-repair and dA-tailing enzymes",
                      MaterialKind.CATALYST)
    tailed = g.product("tailed", "end-repaired dA-tailed duplex cDNA",
                       MaterialKind.MOLECULAR)
    g.react("end_prep", "End-repair and dA-tail the duplex cDNA", StepKind.END_REPAIR,
            ((duplex, InputRole.TEMPLATE), (end_mix, InputRole.CATALYST)), (tailed,))
    adapter = g.input("adapter", "motor-loaded dT sequencing adapter",
                      MaterialKind.ADAPTER, inferred=True,
                      source=SOURCE + " (adapter bases proprietary)")
    ligase = g.input("ligase", "sequencing-adapter ligase", MaterialKind.CATALYST)
    library = g.product("library", "motor-adapted PCR-free cDNA library",
                        MaterialKind.LIBRARY, state=FINAL_LIBRARY, inferred=True,
                        source=SOURCE + " (adapter bases proprietary)")
    g.react("ligate", "Ligate the motor-loaded dT adapter", StepKind.LIGATION,
            ((tailed, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
             (ligase, InputRole.CATALYST)), (library,), geometry=ADAPTER_END)
    g.sequence(library, motor_sequencing_plan(
        FINAL_LIBRARY, limitation=SEQUENCING_UNAVAILABLE))
    g.declare_reviewed_scope(("rna", "rt_primer", "rt", "reverse_transcribe", "rnase",
                              "remove_rna", "second_strand_polymerase", "second_strand",
                              "end_prep_mix", "end_prep", "adapter", "ligase", "ligate",
                              "library"), SOURCE)
    return g.require_publishable()
