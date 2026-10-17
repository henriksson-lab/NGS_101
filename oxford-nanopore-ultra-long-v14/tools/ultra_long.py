"""Molecular model for Oxford Nanopore SQK-ULK114 ultra-long preparation."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Scene, Segment
from protocol_graph import (Disposition, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind, Terminal)
from protocol_model import molecule, role as molecular_role


def role(name: str, text: str, tag: str | None = None) -> Segment:
    return Segment(name, text, tag, placeholder=True)


def single_strand_rows(construct: Construct, label: str = "component roles"):
    sc = Scene()
    sc.strand(construct.name, list(construct), label=label)
    return sc.rows()

TITLE = "Oxford Nanopore Ultra-Long DNA Sequencing V14 — SQK-ULK114"
NOTES = "01_ultra-long-v14.html"
SOURCE = ('Commercial defining source: Oxford Nanopore '
          '<a href="https://nanoporetech.com/document/ultra-long-dna-sequencing-kit-sqk-ulk114">'
          'ULK_9177_v114</a>.')
SUMMARY = ("Diluted transposase chemistry introduces very few tags into uHMW DNA, a "
           "motor-loaded Rapid Adapter attaches, and precipitation recovers sequenceable "
           "molecules while minimizing mechanical breakage.")
CAVEAT = ("ONT publishes the Rapid Adapter strand. Fragmentation-tag bases, the complete "
          "motor assembly, attachment chemistry and exact motor-bearing strand are not "
          "assembled in this model. The endpoint shows component roles, without placing "
          "the published strand across an unverified attachment junction.")
# Explicit local assignment from the ONT Chemistry Technical Document,
# CHTD_500 v1 rev AV (24 July 2026), section 15, Rapid Adapter (RA) top strand.
# The ULK114 protocol names the RA reagent; sequence identity with another kit's
# strand does not assign that kit's bottom strand, end geometry or motor assembly.
RAPID_ADAPTER = "CCTGTACTTCGTTCAGTTACGTATTGCT"
CHEMISTRY_SOURCE = "https://oxfordnanoporedx.com/document/chemistry-technical-document"


def rapid_adapter_strand() -> Construct:
    return Construct([Segment("published RA top strand", RAPID_ADAPTER, "r2")],
                     name="Rapid Adapter: source-printed strand only")


TAGGED = Construct([
    role("fragmentation tag", "[fragmentation tag]", "me"),
    role("ultra-long genomic DNA", "[ultra-long genomic DNA]"),
], name="sparsely tagged uHMW DNA")
FINAL_LIBRARY = Construct([
    role("motor-loaded Rapid Adapter", "[motor-loaded Rapid Adapter]", "r2"),
    role("tagged uHMW DNA", "[transposase-tagged ultra-long genomic DNA]", "me"),
], name="SQK-ULK114 logical sequencing endpoint")
SEQ_PRIMERS = ()
SEQUENCING_MECHANISM = "nanopore_motor"
MOTOR_ADAPTER_NAME = "Rapid Adapter"
SEQUENCING_ENDING = ("No synthesis primer binds. The motor carried by the Rapid Adapter "
                     "docks at the nanopore and meters the attached DNA strand through the "
                     "pore, adapter-bearing end first.")
SEQUENCING_UNAVAILABLE = ("Nanopore sequencing uses motor-controlled translocation, not "
                          "primer extension; the RA strand is published but the complete "
                          "motor/attachment structure is not assembled in this model.")
FINAL_CAPTION = ("Logical component relationship only. The source-printed RA strand is "
                 "listed separately; its physical attachment to the FRA tag and the "
                 "motor-bearing assembly are not reconstructed here. The protocol "
                 "supplies no ligase for this attachment.")
FRA_DILUTION_UL = (6, 244)
TAGMENTATION_PROGRAM = (("room temperature", 10), (75, 10), ("ice", 10))
RAPID_ADAPTER_MIN = 30


def final_rows():
    # An explicitly incomplete component assembly, not a covalent strand.
    return single_strand_rows(FINAL_LIBRARY)


def sections():
    return [
        ("Extract and preserve ultra-high-molecular-weight DNA", [
            *single_strand_rows(Construct([role("uHMW genomic DNA", "[ultra-long genomic DNA]")],
                                         name="uHMW genomic DNA"))],
         "Extract DNA from six million cells and handle it with wide-bore tips and slow mixing."),
        ("Tagment with highly diluted Fragmentation Mix", single_strand_rows(TAGGED),
         "Mix 6 µl FRA with 244 µl dilution buffer, add to 750 µl DNA, then incubate 10 min at room temperature, 10 min at 75 °C and at least 10 min on ice."),
        ("Attach Rapid Adapter", final_rows(),
         "Add 5 µl Rapid Adapter and incubate 30 min at room temperature; ONT describes enzyme-free attachment to the transposase tag, without resolving the exact junction."),
        ("Precipitate and elute without shearing", final_rows(),
         "Precipitation Buffer and a Precipitation Star recover the viscous library; elute overnight before loading."),
    ]


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE)
    dna = g.input("dna", "ultra-high-molecular-weight genomic DNA",
                  MaterialKind.MOLECULAR,
                  state=molecule("uHMW DNA", (molecular_role("ultra-long genomic DNA"),)),
                  primary=True)
    tags = g.input("tag_complex", "highly diluted Fragmentation Mix tagging complex",
                   MaterialKind.ADAPTER, inferred=True,
                   source=SOURCE + " (fragmentation-tag bases unresolved)")
    transposase = g.input("transposase", "Fragmentation Mix transposase",
                          MaterialKind.CATALYST)
    tagged = g.product("tagged", "sparsely tagged uHMW DNA",
                       MaterialKind.MOLECULAR, state=TAGGED, inferred=True,
                       source=SOURCE + " (tag bases unresolved)")
    fdb = g.input("fra_dilution_buffer", "FRA Dilution Buffer (FDB)",
                  MaterialKind.CONDITION, source=SOURCE + " (6 µl FRA + 244 µl FDB)")
    off_products = g.product("off_products",
                             "candidate tagmentation products without a productive attachment end",
                             MaterialKind.MOLECULAR, inferred=True)
    g.react("tag", "Sparsely tagment ultra-long DNA", StepKind.TAGMENTATION,
            ((dna, InputRole.TEMPLATE), (tags, InputRole.STRUCTURAL),
             (transposase, InputRole.CATALYST), (fdb, InputRole.CONDITION)),
            (tagged, off_products))
    g.mark_terminal(Terminal(off_products.id, Disposition.UNRESOLVED,
        "Candidate class, not an observed population. Precipitation recovers DNA without "
        "an adapter-specific selection, so unadapted DNA is not shown as physically removed."))
    rapid = g.input("rapid_adapter", "motor-loaded Rapid Adapter",
                    MaterialKind.ADAPTER, state=rapid_adapter_strand(),
                    source=CHEMISTRY_SOURCE + "; source-printed strand only; complete reagent assembly not modelled")
    adapted = g.product("adapted", "motor-adapted ultra-long library",
                        MaterialKind.MOLECULAR, state=FINAL_LIBRARY, inferred=True,
                        source=SOURCE + " (complete motor-bearing assembly unresolved)")
    g.react("attach", "Attach the motor-loaded Rapid Adapter",
            StepKind.ADAPTER_ATTACHMENT,
            ((tagged, InputRole.TEMPLATE), (rapid, InputRole.STRUCTURAL)), (adapted,))
    precipitant = g.input("precipitant", "Precipitation Buffer (PTB) and Precipitation Star (PS)",
                          MaterialKind.SELECTION)
    elution = g.input("elution_buffer", "Elution Buffer (EB), overnight elution",
                      MaterialKind.CONDITION, source=SOURCE + " (300 µl EB, ≥12 h at room temperature)")
    library = g.product("library", "recovered ultra-long Nanopore library",
                        MaterialKind.LIBRARY, state=FINAL_LIBRARY, inferred=True,
                        source=SOURCE + " (complete motor-bearing assembly unresolved)")
    rejected = g.product("rejected", "removed supernatant (unprecipitated material)",
                         MaterialKind.WASTE)
    g.react("recover", "Precipitate on the star and elute overnight",
            StepKind.CLEANUP,
            ((adapted, InputRole.TEMPLATE), (precipitant, InputRole.SELECTION),
             (elution, InputRole.CONDITION)),
            (library, rejected), changes_molecule=False)
    g.discard(rejected, "supernatant aspirated after DNA precipitates around the Precipitation Star")
    g.sequence(library, SequencingPlan(
        library=FINAL_LIBRARY, primers=(), acquired_roles=("Nanopore strand",),
        read_uses={"Nanopore strand": "include"}, limitation=SEQUENCING_UNAVAILABLE))
    g.declare_reviewed_scope(("dna", "tag_complex", "transposase", "fra_dilution_buffer",
                              "tag", "rapid_adapter", "attach", "precipitant",
                              "elution_buffer", "recover", "library"), SOURCE)
    return g.require_publishable()


def oligos():
    return [f"Rapid Adapter strand: 5′-{RAPID_ADAPTER}-3′"]
