"""Topology model for PacBio SMRTbell Express Template Prep Kit 2.0."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Row, Segment, annotation_rows, revcomp, strand_row
from dumbbell import Dumbbell, dumbbell_rows
from endprep import dA_tailed_scene, repair_and_dA_tail

TITLE = "PacBio SMRTbell Express Template Prep Kit 2.0"
NOTES = "01_smrtbell-express-2-0.html"
SOURCE = ('Commercial defining source: PacBio procedure <b>101-693-800 Version 01</b>, '
          '“Preparing gDNA Libraries Using the SMRTbell Express Template Preparation Kit 2.0.”')
SUMMARY = ("Remove overhangs and repair sheared long DNA, dA-tail it, then ligate a "
           "T-overhang hairpin to both ends to make a closed template for repeated "
           "polymerase passes.")
CAVEAT = ("Overhang Adapter v3 and sequencing-primer bases are proprietary. The hairpin "
          "segments below encode their disclosed topology, not invented sequence.")

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
                 "the four sealed insert–adapter junctions; hairpin bases are proprietary.")
OVERHANG_REMOVAL = (37, 15)
DAMAGE_REPAIR = (37, 30)
END_PREP_PROGRAM = ((20, 10), (65, 30))
LIGATION = (20, 60)


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
         "T-overhang hairpins pair with both insert dA ends and are sealed for 60 min at 20 °C."),
        ("Remove molecules with free ends", dumbbell_rows(SMRTBELL),
         "Bead cleanup and nuclease treatment enrich completely closed SMRTbells."),
        ("Anneal sequencing primer and bind polymerase", primer_path,
         "The proprietary primer site is present in each hairpin; the path shown is computed from the closed template."),
    ]
