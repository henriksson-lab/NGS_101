"""Topology model for PacBio SMRTbell prep kit 3.0 libraries."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Segment
from dumbbell import (Dumbbell, smrtbell_circuit, smrtbell_hairpin_geometry,
                      smrtbell_sequencing_plan)
from endprep import dA_tailed_scene, repair_and_dA_tail
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


INSERT = Construct([seg("genomic insert", "X" * 48, placeholder=True)],
                   name="15–20 kb genomic insert")
END_PREP = repair_and_dA_tail(INSERT)
# Full current adapter bases are proprietary. N is molecularly unknown DNA, not a guessed
# sequence; the terminal T overhang is represented separately by Dumbbell.
HAIRPIN_PLACEHOLDER = "N" * 18
SMRTBELL = Dumbbell(INSERT, HAIRPIN_PLACEHOLDER, HAIRPIN_PLACEHOLDER,
                    name="SMRTbell library", insert_overhang="A", adapter_overhang="T",
                    # Derived from the dA-tailed insert; neither inspected source prints it.
                    overhang_inferred=True)
END_PREP_TIMES = ((37, 30), (65, 5))
LIGATION = (20, 30)
NUCLEASE = (37, 15)
# Procedure 102-166-600 Rev09 steps 5.9-5.22: 95 uL (1X) SMRTbell cleanup beads after
# ligation; steps 7.1-7.14: 3.1X (155 uL) of 35% v/v AMPure PB after nuclease treatment.
POST_LIGATION_BEADS_UL = 95
SIZE_SELECTION_BEADS = (3.1, 35)
FINAL_LIBRARY = None
SEQUENCING_UNAVAILABLE = ("The kit guide names the sequencing primer and SMRTbell "
                          "adapter but does not disclose their current bases.")
HAIRPIN_END = smrtbell_hairpin_geometry(
    SMRTBELL,
    statement="Procedure 102-166-600 Rev09 section 4 (Repair and A-tailing) and training "
    "102-390-900 Rev03 p. 6 (3-prime A drawn on both insert strands; adapter design "
    "'Overhang') establish dA-tailed inserts and an overhang adapter. The adapter's "
    "complementary 3-prime T is inferred from that pairing; neither inspected source "
    "prints the adapter terminal base or body.")


def end_prep_scene():
    return dA_tailed_scene(END_PREP)


def polymerase_cycle() -> Construct:
    return smrtbell_circuit(SMRTBELL)


FINAL_LIBRARY = polymerase_cycle()


def protocol_graph():
    authority = ('Commercial protocol: PacBio WGS/metagenome procedure '
                 '<a href="https://www.pacb.com/wp-content/uploads/Procedure-checklist-Preparing-whole-genome-and-metagenome-libraries-using-SMRTbell-prep-kit-3.0.pdf">102-166-600 Rev09</a>.')
    g = ProtocolGraph("PacBio SMRTbell prep kit 3.0", authority=authority)
    dna = g.input("dna", "sheared 15–20 kb genomic DNA", MaterialKind.MOLECULAR,
                  state=INSERT, primary=True)
    repair = g.input("repair", "damage-repair, end-repair and dA-tailing enzymes",
                     MaterialKind.CATALYST)
    tailed = g.product("tailed", "repaired dA-tailed genomic inserts",
                       MaterialKind.MOLECULAR, state=END_PREP)
    g.react("end_prep", "Repair DNA damage and ends, and dA-tail",
            StepKind.END_REPAIR,
            ((dna, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (tailed,))
    hairpins = g.input("hairpins", "SMRTbell overhang hairpin adapter (or index plate 96A adapter)",
                       MaterialKind.ADAPTER, inferred=True,
                       source=authority + " (adapter bases proprietary; 3-prime T inferred)")
    ligase = g.input("ligase", "SMRTbell ligase", MaterialKind.CATALYST)
    closed = g.product("closed", "closed SMRTbell molecules", MaterialKind.MOLECULAR,
                       state=FINAL_LIBRARY, inferred=True,
                       source=authority + " (hairpin bases proprietary)")
    g.react("ligate", "Ligate a hairpin to each insert end", StepKind.LIGATION,
            ((tailed, InputRole.TEMPLATE), (hairpins, InputRole.STRUCTURAL),
             (ligase, InputRole.CATALYST)), (closed,), geometry=HAIRPIN_END)
    beads1 = g.input("cleanup_beads", "1X SMRTbell cleanup beads (95 µL)",
                     MaterialKind.CONDITION)
    ligated = g.product("ligated", "bead-purified ligation products",
                        MaterialKind.MOLECULAR, state=FINAL_LIBRARY, inferred=True,
                        source=authority + " (hairpin bases proprietary)")
    lig_sup = g.product("ligation_supernatant", "ligation supernatant",
                        MaterialKind.WASTE)
    g.react("ligation_cleanup", "Purify ligation with 1X SMRTbell cleanup beads",
            StepKind.CLEANUP,
            ((closed, InputRole.TEMPLATE), (beads1, InputRole.CONDITION)),
            (ligated, lig_sup), changes_molecule=False,
            reason="Procedure 102-166-600 Rev09 steps 5.9-5.22; supernatant discarded.")
    g.discard(lig_sup, "discarded supernatant of the 1X SMRTbell bead cleanup")
    nuclease = g.input("nuclease", "Nuclease buffer and Nuclease mix",
                       MaterialKind.CATALYST)
    treated = g.product("treated", "nuclease-resistant closed SMRTbells",
                        MaterialKind.MOLECULAR, state=FINAL_LIBRARY, inferred=True,
                        source=authority + " (hairpin bases proprietary)")
    linear = g.product("linear", "unligated DNA fragments and leftover SMRTbell adapter",
                       MaterialKind.WASTE)
    g.react("remove_linear", "Nuclease treatment, 37 °C 15 min", StepKind.DIGESTION,
            ((ligated, InputRole.TEMPLATE), (nuclease, InputRole.CATALYST)),
            (treated, linear), changes_molecule=False)
    g.discard(linear, "digested by the nuclease treatment (procedure section 6)")
    ampure = g.input("ampure", "3.1X (155 µL) 35% v/v AMPure PB beads",
                     MaterialKind.CONDITION)
    library = g.product("library", "size-selected SMRTbell library",
                        MaterialKind.LIBRARY, state=FINAL_LIBRARY, inferred=True,
                        source=authority + " (hairpin bases proprietary)")
    short = g.product("short_fraction",
                      "supernatant: nucleases, contaminants and depleted fragments <5 kb",
                      MaterialKind.WASTE)
    g.react("size_select", "Diluted AMPure PB bead cleanup and size selection",
            StepKind.CLEANUP,
            ((treated, InputRole.TEMPLATE), (ampure, InputRole.CONDITION)),
            (library, short), changes_molecule=False,
            reason="Procedure 102-166-600 Rev09 section 7 and step 6.8; libraries <10 kb "
                   "use a 1X SMRTbell bead cleanup instead (not drawn).")
    g.discard(short, "not retained by the 3.1X 35% AMPure PB beads")
    g.sequence(library, smrtbell_sequencing_plan(
        SMRTBELL, FINAL_LIBRARY, limitation=SEQUENCING_UNAVAILABLE))
    g.declare_reviewed_scope(("dna", "repair", "end_prep", "hairpins", "ligase",
                              "ligate", "cleanup_beads", "ligation_cleanup", "nuclease",
                              "remove_linear", "ampure", "size_select", "library"),
                             authority)
    return g.require_publishable()
