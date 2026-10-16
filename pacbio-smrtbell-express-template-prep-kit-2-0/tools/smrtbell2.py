"""Topology model for PacBio SMRTbell Express Template Prep Kit 2.0."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Row, Segment, annotation_rows, revcomp, strand_row
from dumbbell import (Dumbbell, dumbbell_rows, smrtbell_hairpin_geometry,
                      smrtbell_sequencing_plan)
from endprep import dA_tailed_scene, repair_and_dA_tail
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind

TITLE = "PacBio SMRTbell Express Template Prep Kit 2.0"
NOTES = "01_smrtbell-express-2-0.html"
SOURCE = ('Commercial defining source: PacBio procedure <b>101-693-800 Version 01</b>, '
          '“Preparing gDNA Libraries Using the SMRTbell Express Template Preparation Kit 2.0”; '
          '<a href="https://www.pacb.com/wp-content/uploads/SMRTbell-Express-TPK-2.0-Overview-Large-Insert-gDNA-Library-Prep-Customer-Training.pdf">'
          'PacBio training PN 101-742-500 Version 01</a>, pp. 7 and 9.')
SUMMARY = ("Remove overhangs and repair sheared long DNA, dA-tail it, then ligate a "
           "hairpin adapter to both ends to make a closed template for repeated "
           "polymerase passes.")
CAVEAT = ("PacBio's February 2019 training document identifies a T overhang on Overhang "
          "Adapter v3 and A/T ligation. The hairpin body, complete adapter sequence and "
          "sequencing-primer bases remain unverified in the inspected sources; the "
          "hairpin segments encode closed topology with placeholder bases and length.")

INSERT = Construct([Segment("long genomic insert", "X" * 48, placeholder=True)],
                   name="10 kb to >30 kb insert")
END_PREP = repair_and_dA_tail(INSERT)
HAIRPIN = "N" * 18
SMRTBELL = Dumbbell(INSERT, HAIRPIN, HAIRPIN, name="SMRTbell 2.0 library",
                    insert_overhang="A", adapter_overhang="T")
FINAL_LIBRARY = Construct([
    Segment("insert strand 1", INSERT.top(), placeholder=True),
    Segment("right hairpin / primer site", HAIRPIN, "r1", placeholder=True),
    Segment("insert strand 2 reverse complement", revcomp(INSERT.top()), placeholder=True),
    Segment("left hairpin / primer site", HAIRPIN, "r1", placeholder=True),
], name="one closed polymerase circuit")
if FINAL_LIBRARY.top() != SMRTBELL.pass_sequence:
    raise AssertionError("the drawn polymerase circuit must follow the closed dumbbell")

SEQ_PRIMERS = ()
SEQUENCING_ENDING = ("A proprietary PacBio sequencing primer anneals within either "
                     "hairpin adapter; polymerase follows the derived closed circuit "
                     "through both insert strands for repeated consensus passes.")
SEQUENCING_UNAVAILABLE = ("The kit guide names its sequencing primer but does not print "
                          "the primer or Overhang Adapter v3 nucleotide sequence.")
FINAL_CAPTION = ("The ds insert plus two hairpins is one covalently closed path. ** marks "
                 "the four sealed insert–adapter junctions; hairpin bases and drawn "
                 "length are schematic placeholders, while the ligating T is source-reported.")
OVERHANG_REMOVAL = (37, 15)
DAMAGE_REPAIR = (37, 30)
END_PREP_PROGRAM = ((20, 10), (65, 30))
LIGATION = (20, 60)
HAIRPIN_END = smrtbell_hairpin_geometry(
    SMRTBELL, statement="PacBio training PN 101-742-500 v01 p. 7 states that Overhang "
    "Adapter v3 has a T overhang; p. 9 specifies A/T ligation. Procedure 101-693-800 "
    "v01 specifies dA-tailed insert ends. The 3-prime pairing follows from these ends; "
    "the proprietary hairpin body remains unspecified.")


def final_rows():
    return dumbbell_rows(SMRTBELL)


def sections():
    primer_path = [strand_row(FINAL_LIBRARY, "top",
                              prefix="polymerase ---> 5'- ", suffix=" -3'"),
                   *annotation_rows(FINAL_LIBRARY)]
    return [
        ("Remove single-stranded overhangs and repair damage",
         [Row(chunks=[("sheared long duplex  -- overhang removal -->  repaired blunt duplex", None, False)])],
         "DNA Prep Enzyme: 37 °C for 15 min; DNA Damage Repair Mix v2: 37 °C for 30 min."),
        ("End-repair and dA-tail", dA_tailed_scene(END_PREP).rows(),
         "End Prep Mix: 20 °C for 10 min, then 65 °C for 30 min."),
        ("Ligate Overhang Adapter v3 hairpins", dumbbell_rows(SMRTBELL),
         "Overhang Adapter v3 is ligated for 60 min at 20 °C. PacBio training "
         "PN 101-742-500 v01 pp. 7, 9 identifies the T overhang and A/T ligation; "
         "the full hairpin sequence remains undisclosed."),
        ("Remove molecules with free ends", dumbbell_rows(SMRTBELL),
         "Bead cleanup and nuclease treatment enrich completely closed SMRTbells."),
        ("Anneal sequencing primer and bind polymerase", primer_path,
         "The proprietary primer site is present in each hairpin; the path shown is computed from the closed template."),
    ]


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE)
    dna = g.input("dna", "sheared long genomic DNA", MaterialKind.MOLECULAR,
                  state=INSERT, primary=True)
    repair = g.input("repair", "overhang-removal, damage-repair and end-prep enzymes",
                     MaterialKind.CATALYST)
    tailed = g.product("tailed", "repaired dA-tailed long DNA",
                       MaterialKind.MOLECULAR, state=END_PREP)
    g.react("end_prep", "Remove overhangs, repair damage and dA-tail",
            StepKind.END_REPAIR,
            ((dna, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (tailed,))
    hairpins = g.input("hairpins", "Overhang Adapter v3 hairpins",
                       MaterialKind.ADAPTER, inferred=True,
                       source=SOURCE + " (T overhang known; body bases proprietary)")
    ligase = g.input("ligase", "SMRTbell ligase", MaterialKind.CATALYST)
    closed = g.product("closed", "closed SMRTbell molecules", MaterialKind.MOLECULAR,
                       state=FINAL_LIBRARY, inferred=True,
                       source=SOURCE + " (hairpin bases proprietary)")
    g.react("ligate", "Ligate a hairpin to each dA-tailed insert end",
            StepKind.LIGATION,
            ((tailed, InputRole.TEMPLATE), (hairpins, InputRole.STRUCTURAL),
             (ligase, InputRole.CATALYST)), (closed,), geometry=HAIRPIN_END)
    nuclease = g.input("nuclease", "free-end nuclease", MaterialKind.CATALYST)
    library = g.product("library", "nuclease-resistant closed SMRTbell library",
                        MaterialKind.LIBRARY, state=FINAL_LIBRARY, inferred=True,
                        source=SOURCE + " (hairpin bases proprietary)")
    linear = g.product("linear", "molecules retaining a free DNA end",
                       MaterialKind.WASTE)
    g.react("remove_linear", "Digest molecules that retain a free end",
            StepKind.DIGESTION,
            ((closed, InputRole.TEMPLATE), (nuclease, InputRole.CATALYST)),
            (library, linear), changes_molecule=False)
    g.discard(linear, "not a completely closed SMRTbell")
    g.sequence(library, smrtbell_sequencing_plan(
        SMRTBELL, FINAL_LIBRARY, limitation=SEQUENCING_UNAVAILABLE))
    g.declare_reviewed_scope(("dna", "repair", "end_prep", "hairpins", "ligase",
                              "ligate", "nuclease", "remove_linear", "library"), SOURCE)
    return g.require_publishable()
