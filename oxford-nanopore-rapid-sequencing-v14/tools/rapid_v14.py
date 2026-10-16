"""Oxford Nanopore SQK-RAD114 rapid genomic-DNA workflow."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Row, Scene, Segment, annotation_rows
from nanopore_library import RAPID_ADAPTER, motor_sequencing_plan
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from protocol_model import molecule, role as molecular_role


def role(name: str, text: str, tag: str | None = None) -> Segment:
    return Segment(name, text, tag, placeholder=True)


TITLE = "Oxford Nanopore Rapid Sequencing V14 — SQK-RAD114"
NOTES = "01_rapid-sequencing-v14.html"
SOURCE = ('Commercial defining source: Oxford Nanopore '
          '<a href="https://nanoporetech.com/document/rapid-sequencing-sqk-rad114">'
          'Rapid Sequencing Kit V14 — gDNA, RSE_9177_v114 Rev Q</a>.')
SUMMARY = ("A transposase fragments high-molecular-weight genomic DNA while installing "
           "rapid-attachment tags; the motor-loaded Rapid Adapter then attaches without "
           "ligase and presents one tagged strand to a nanopore.")
CAVEAT = ("ONT publishes the Rapid Adapter strand. The transposase tag, complete "
          "motor-bearing assembly, attachment chemistry and exact motor-bearing strand "
          "have not been verified from the checked vendor documents.")
TAGGED = Construct([role("left transposase tag", "[transposase tag]", "me"),
                    role("genomic DNA", "[genomic DNA fragment]"),
                    role("right transposase tag", "[transposase tag]", "me")],
                   name="transposase-tagged genomic DNA")
FINAL_LIBRARY = Construct([
    role("motor-loaded Rapid Adapter assembly", "[motor assembly]", "r2"),
    Segment("Rapid Adapter published strand", RAPID_ADAPTER, "r2"),
    role("attachment boundary", "[attachment]", "me"),
    role("transposase tag", "[transposase tag]", "me"),
    role("genomic DNA", "[genomic DNA fragment]"),
], name="SQK-RAD114 sequenceable strand")
SEQ_PRIMERS = ()
SEQUENCING_MECHANISM = "nanopore_motor"
MOTOR_ADAPTER_NAME = "Rapid Adapter"
SEQUENCING_ENDING = ("No synthesis primer binds. The adapter motor docks at a nanopore "
                     "and meters the attached DNA strand through the pore, adapter-bearing "
                     "end first.")
SEQUENCING_UNAVAILABLE = ("SQK-RAD114 uses motor-controlled nanopore translocation, not "
                          "primer extension; the RA strand is published but the complete "
                          "motor/attachment structure is not.")
FINAL_CAPTION = ("One sequenceable product. Rapid attachment is non-ligase chemistry, so "
                 "the boundary is labelled attachment rather than covalent ligation.")
FRAGMENTATION_PROGRAM = ((30, 2), (80, 2))
RAPID_ATTACHMENT_MIN = 5


def sections():
    sc = Scene.duplex(list(TAGGED), label="tagged fragment",
                      unpaired=("left transposase tag'", "genomic DNA'",
                                "right transposase tag'"))
    return [
        ("Fragment and tag genomic DNA", sc.rows(),
         "Fragmentation Mix cuts the input and transfers transposase tags in the same "
         "30 °C reaction; 80 °C ends the two-minute tagmentation."),
        ("Attach the Rapid Adapter", [
            Row(chunks=[("[motor assembly] ", "r2", False),
                        (RAPID_ADAPTER, "r2", False),
                        (" ** attachment ** ", None, False),
                        ("[transposase tag]", "me", False),
                        ("[genomic DNA fragment]", None, False)]),
        ], "Diluted RA attaches to a transposase-tagged end for five minutes at room "
           "temperature without a ligation enzyme."),
        ("Enter the nanopore", [
            Row(chunks=[("membrane ─ pore | motor ─ adapter ─ tagged DNA  ─────────>",
                        None, False)]),
            Row(chunks=[("                     no sequencing primer", "r1", False)]),
        ], "The motor, rather than a polymerase primer, controls strand motion through "
           "the R10.4.1 nanopore."),
    ]


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE)
    dna = g.input("dna", "high-molecular-weight genomic DNA", MaterialKind.MOLECULAR,
                  state=molecule("HMW genomic DNA", (molecular_role("genomic DNA"),)),
                  primary=True)
    tags = g.input("tag_complex", "Fragmentation Mix loaded tagging complex",
                   MaterialKind.ADAPTER, inferred=True,
                   source=SOURCE + " (tag bases and loading structure unverified)")
    transposase = g.input("transposase", "Fragmentation Mix transposase",
                          MaterialKind.CATALYST)
    tagged = g.product("tagged", "transposase-tagged genomic fragments",
                       MaterialKind.MOLECULAR, state=TAGGED, inferred=True,
                       source=SOURCE + " (tag bases unverified)")
    off_products = g.product("off_products", "non-sequenceable tagmentation products",
                             MaterialKind.WASTE)
    g.react("tag", "Fragment and tag genomic DNA", StepKind.TAGMENTATION,
            ((dna, InputRole.TEMPLATE), (tags, InputRole.STRUCTURAL),
             (transposase, InputRole.CATALYST)), (tagged, off_products))
    g.discard(off_products, "does not carry a productive Rapid Adapter attachment end")
    rapid = g.input("rapid_adapter", "motor-loaded Rapid Adapter",
                    MaterialKind.ADAPTER, inferred=True,
                    source=SOURCE + f" (RA strand {RAPID_ADAPTER}; motor/attachment unverified)")
    library = g.product("library", "motor-adapted Nanopore library",
                        MaterialKind.LIBRARY, state=FINAL_LIBRARY, inferred=True,
                        source=SOURCE + " (complete motor-bearing assembly unverified)")
    g.react("attach", "Attach the motor-loaded Rapid Adapter",
            StepKind.ADAPTER_ATTACHMENT,
            ((tagged, InputRole.TEMPLATE), (rapid, InputRole.STRUCTURAL)),
            (library,))
    g.sequence(library, motor_sequencing_plan(
        FINAL_LIBRARY, limitation=SEQUENCING_UNAVAILABLE))
    g.declare_reviewed_scope(("dna", "tag_complex", "transposase", "tag",
                              "rapid_adapter", "attach", "library"), SOURCE)
    return g.require_publishable()


def oligos():
    return [f"Rapid Adapter strand: 5′-{RAPID_ADAPTER}-3′"]
