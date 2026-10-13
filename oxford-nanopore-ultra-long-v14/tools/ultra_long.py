"""Molecular model for Oxford Nanopore SQK-ULK114 ultra-long preparation."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct
from nanopore_library import (motor_library, motor_sequencing_plan, no_primer_ending,
                              role, single_strand_rows)
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from protocol_model import molecule, role as molecular_role

TITLE = "Oxford Nanopore Ultra-Long DNA Sequencing V14 — SQK-ULK114"
NOTES = "01_ultra-long-v14.html"
SOURCE = ('Commercial defining source: Oxford Nanopore '
          '<a href="https://nanoporetech.com/document/ultra-long-dna-sequencing-kit-sqk-ulk114">'
          'ULK_9177_v114</a>.')
SUMMARY = ("Diluted transposase chemistry introduces very few tags into uHMW DNA, a "
           "motor-loaded Rapid Adapter attaches, and precipitation recovers sequenceable "
           "molecules while minimizing mechanical breakage.")
CAVEAT = ("Fragmentation-tag and Rapid-Adapter bases, attachment chemistry and the exact "
          "motor-bearing strand are proprietary. Bracketed regions state roles only.")

TAGGED = Construct([
    role("fragmentation tag", "[fragmentation tag]", "me"),
    role("ultra-long genomic DNA", "[ultra-long genomic DNA]"),
], name="sparsely tagged uHMW DNA")
FINAL_LIBRARY = motor_library(list(TAGGED), adapter="Rapid Adapter",
                              name="SQK-ULK114 pore-entry strand")
SEQ_PRIMERS = ()
SEQUENCING_ENDING = no_primer_ending("Rapid Adapter")
SEQUENCING_UNAVAILABLE = ("Nanopore sequencing uses motor-controlled translocation, not "
                          "primer extension; current adapter bases are undisclosed.")
FINAL_CAPTION = ("One sequenceable tagged end. ** marks Rapid Adapter attachment; it is "
                 "not labelled ligation because the protocol supplies no ligase.")
FRA_DILUTION_UL = (6, 244)
TAGMENTATION_PROGRAM = (("room temperature", 10), (75, 10), ("ice", 10))
RAPID_ADAPTER_MIN = 30


def final_rows():
    return single_strand_rows(FINAL_LIBRARY, (
        ("motor-loaded sequencing adapter", "fragmentation tag", "attachment"),
    ))


def sections():
    return [
        ("Extract and preserve ultra-high-molecular-weight DNA", [
            *single_strand_rows(Construct([role("uHMW genomic DNA", "[ultra-long genomic DNA]")],
                                         name="uHMW genomic DNA"))],
         "Extract DNA from six million cells and handle it with wide-bore tips and slow mixing."),
        ("Tagment with highly diluted Fragmentation Mix", single_strand_rows(TAGGED),
         "Mix 6 µl FRA with 244 µl dilution buffer, add to 750 µl DNA, then incubate 10 min at room temperature, 10 min at 75 °C and at least 10 min on ice."),
        ("Attach Rapid Adapter", final_rows(),
         "Add 5 µl Rapid Adapter and incubate 30 min at room temperature; the proprietary attachment creates a motor-bearing pore-entry end."),
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
                   source=SOURCE + " (fragmentation-tag bases proprietary)")
    transposase = g.input("transposase", "Fragmentation Mix transposase",
                          MaterialKind.CATALYST)
    tagged = g.product("tagged", "sparsely tagged uHMW DNA",
                       MaterialKind.MOLECULAR, state=TAGGED, inferred=True,
                       source=SOURCE + " (tag bases proprietary)")
    off_products = g.product("off_products", "cut or unproductively tagged molecules",
                             MaterialKind.WASTE)
    g.react("tag", "Sparsely tagment ultra-long DNA", StepKind.TAGMENTATION,
            ((dna, InputRole.TEMPLATE), (tags, InputRole.STRUCTURAL),
             (transposase, InputRole.CATALYST)), (tagged, off_products))
    g.discard(off_products, "does not preserve a sequenceable ultra-long tagged molecule")
    rapid = g.input("rapid_adapter", "motor-loaded Rapid Adapter",
                    MaterialKind.ADAPTER, inferred=True,
                    source=SOURCE + " (adapter bases and attachment chemistry proprietary)")
    adapted = g.product("adapted", "motor-adapted ultra-long library",
                        MaterialKind.MOLECULAR, state=FINAL_LIBRARY, inferred=True,
                        source=SOURCE + " (adapter bases proprietary)")
    g.react("attach", "Attach the motor-loaded Rapid Adapter",
            StepKind.ADAPTER_ATTACHMENT,
            ((tagged, InputRole.TEMPLATE), (rapid, InputRole.STRUCTURAL)), (adapted,))
    precipitant = g.input("precipitant", "Precipitation Buffer and Star",
                          MaterialKind.SELECTION)
    library = g.product("library", "recovered ultra-long Nanopore library",
                        MaterialKind.LIBRARY, state=FINAL_LIBRARY, inferred=True,
                        source=SOURCE + " (adapter bases proprietary)")
    rejected = g.product("rejected", "unrecovered or sheared molecules",
                         MaterialKind.WASTE)
    g.react("recover", "Precipitate and recover motor-adapted molecules",
            StepKind.SELECTION,
            ((adapted, InputRole.TEMPLATE), (precipitant, InputRole.SELECTION)),
            (library, rejected), changes_molecule=False)
    g.discard(rejected, "not recovered in the ultra-long library fraction")
    g.sequence(library, motor_sequencing_plan(
        FINAL_LIBRARY, limitation=SEQUENCING_UNAVAILABLE))
    g.declare_reviewed_scope(("dna", "tag_complex", "transposase", "tag",
                              "rapid_adapter", "attach", "precipitant", "recover",
                              "library"), SOURCE)
    return g.require_publishable()
