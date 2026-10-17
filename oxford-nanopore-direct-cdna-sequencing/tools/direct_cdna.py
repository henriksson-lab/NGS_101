"""Oxford Nanopore SQK-DCS109 direct-cDNA workflow."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from adapter_ligation import reported_adapter_end
from chemdraw import Segment
from nanopore_library import motor_library, motor_sequencing_plan, role
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from protocol_model import molecule, role as molecular_role

TITLE = "Oxford Nanopore Direct cDNA Sequencing — SQK-DCS109"
NOTES = "01_direct-cdna.html"
SOURCE = ('Commercial protocol: Oxford Nanopore <a href="https://nanoporetech.com/document/chemistry-technical-document">Chemistry Technical Document</a> '
          '(CHTD_500_v1_revAV_24July2026, historical "Direct cDNA Sequencing Kit" SQK-DCS109 section) and kit comparison in the Oxford Nanopore '
          '<a href="https://store.nanoporetech.jp/sites/default/files/s3/literature/product-brochure.pdf">product brochure</a> (BR_1009(EN)_V5_30Nov2020, "Which RNA kit?").')
SUMMARY = "Poly(A) RNA is reverse-transcribed with strand switching, RNA is removed and second-strand DNA is synthesized. The PCR-free duplex cDNA is end-repaired, dA-tailed and ligated to a motor-loaded dT-overhang sequencing adapter."
CAVEAT = ("SQK-DCS109 is a historical commercial protocol. Its kit-specific protocol, with reagent names and conditions, was not obtained, so RT-primer, strand-switching-primer and sequencing-adapter bases are unverified. "
          "The technical document's adapter-sequence table lists sequences \"for our current kits\" (LA, NA, RA, RTP, SSPII and others); none is assigned to SQK-DCS109, so none is used here. "
          "The vendor gives the 100 ng poly(A) RNA input (or 70–200 ng of prepared cDNA) and the order: RT primer annealing, strand-switching RT, RNase-cocktail RNA removal, second-strand synthesis, dA tailing and dT-tailed adapter ligation. "
          "End repair is inferred from the statement that the duplex is \"treated similarly to gDNA\". No PCR or rapid-attachment tag is introduced, unlike cDNA-PCR.")
FINAL_LIBRARY = motor_library([
    role("dT sequencing-adapter end", "[dT sequencing-adapter end]", "r2"),
    Segment("dT:dA ligation junction", "T", bottom="A"),
    role("PCR-free duplex cDNA", "[full-length duplex cDNA]"),
], adapter="sequencing adapter", name="SQK-DCS109 pore-entry strand")
SEQ_PRIMERS = ()
SEQUENCING_MECHANISM = "nanopore_motor"
SEQUENCING_ENDING = "There is no synthesis sequencing primer. A motor in the ligated sequencing adapter presents one strand of the PCR-free duplex cDNA to the pore. The adapter-ligated end enters first."
SEQUENCING_UNAVAILABLE = ("Nanopore sequencing is motor-driven; SQK-DCS109 adapter bases "
                          "and the motor-bearing strand remain unverified pending the kit-specific guide.")
FINAL_CAPTION = "The motor-bearing dT adapter is ligated to dA-tailed PCR-free duplex cDNA."
ADAPTER_END = reported_adapter_end(
    "SQK-DCS109 sequencing adapter", kind="cohesive", polarity="3-prime",
    terminal_sequence="T",
    statement="The historical kit uses a dT-overhang sequencing adapter on dA-tailed cDNA.")

def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE)
    rna = g.input("rna", "poly(A) RNA", MaterialKind.MOLECULAR,
                  state=molecule("poly(A) RNA", (molecular_role("transcript"),
                                                  molecular_role("poly(A) tail"))),
                  primary=True)
    primer = g.input("rt_primer", "annealed RT primer (poly(A) priming inferred; bases unverified)",
                     MaterialKind.OLIGO, inferred=True,
                     source=SOURCE + " (primer bases unverified)")
    ssp = g.input("strand_switch_primer", "strand-switching primer (bases unverified)",
                  MaterialKind.OLIGO, inferred=True,
                  source=SOURCE + " (strand-switching primer added before RT; bases unverified)")
    rt = g.input("rt", "strand-switching reverse transcriptase", MaterialKind.CATALYST)
    hybrid = g.product("hybrid", "full-length RNA:cDNA hybrid", MaterialKind.MOLECULAR)
    g.react("reverse_transcribe", "Make first-strand cDNA by strand switching",
            StepKind.REVERSE_TRANSCRIPTION,
            ((rna, InputRole.TEMPLATE), (primer, InputRole.STRUCTURAL),
             (ssp, InputRole.STRUCTURAL), (rt, InputRole.CATALYST)), (hybrid,))
    rnase = g.input("rnase", "RNase cocktail", MaterialKind.CATALYST)
    single = g.product("single_cdna", "single-stranded full-length cDNA",
                       MaterialKind.MOLECULAR)
    g.react("remove_rna", "Digest the RNA template", StepKind.DIGESTION,
            ((hybrid, InputRole.TEMPLATE), (rnase, InputRole.CATALYST)), (single,))
    second_pol = g.input("second_strand_polymerase", "second-strand DNA polymerase",
                         MaterialKind.CATALYST)
    duplex = g.product("duplex", "full-length duplex cDNA", MaterialKind.MOLECULAR)
    g.react("second_strand", "Synthesize the second cDNA strand", StepKind.EXTENSION,
            ((single, InputRole.TEMPLATE), (second_pol, InputRole.CATALYST)), (duplex,))
    end_mix = g.input("end_prep_mix", "dA-tailing (end-prep) enzymes; repair component inferred",
                      MaterialKind.CATALYST)
    tailed = g.product("tailed", "end-repaired dA-tailed duplex cDNA",
                       MaterialKind.MOLECULAR)
    g.react("end_prep", "End-prep and dA-tail the duplex cDNA (treated like gDNA)", StepKind.END_REPAIR,
            ((duplex, InputRole.TEMPLATE), (end_mix, InputRole.CATALYST)), (tailed,))
    adapter = g.input("adapter", "motor-loaded dT sequencing adapter",
                      MaterialKind.ADAPTER, inferred=True,
                      source=SOURCE + " (adapter bases unverified)")
    ligase = g.input("ligase", "sequencing-adapter ligase", MaterialKind.CATALYST)
    library = g.product("library", "motor-adapted PCR-free cDNA library",
                        MaterialKind.LIBRARY, state=FINAL_LIBRARY, inferred=True,
                        source=SOURCE + " (adapter bases unverified)")
    g.react("ligate", "Ligate the motor-loaded dT adapter", StepKind.LIGATION,
            ((tailed, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
             (ligase, InputRole.CATALYST)), (library,), geometry=ADAPTER_END)
    g.sequence(library, motor_sequencing_plan(
        FINAL_LIBRARY, limitation=SEQUENCING_UNAVAILABLE))
    g.declare_reviewed_scope(("rna", "rt_primer", "strand_switch_primer", "rt", "reverse_transcribe", "rnase",
                              "remove_rna", "second_strand_polymerase", "second_strand",
                              "end_prep_mix", "end_prep", "adapter", "ligase", "ligate",
                              "library"), SOURCE)
    return g.require_publishable()
